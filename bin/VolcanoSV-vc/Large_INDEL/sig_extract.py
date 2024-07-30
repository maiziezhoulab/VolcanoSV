#!/usr/bin/env python
# This script was largly adapted from:
# ''' 
#  * All rights Reserved, Designed By HIT-Bioinformatics   
#  * @Title: cuteSV 
#  * @author: tjiang & sqcao
#  * @date: Nov 3rd 2022
#  * @version V2.0.2
# '''

import pysam
import cigar
from multiprocessing import Pool
import os
import argparse
import logging
import sys
import time
import gc
import subprocess, signal

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm
    
def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

def exe(cmd, timeout=-1):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                            stderr=subprocess.STDOUT, close_fds=True,\
                            preexec_fn=os.setsid)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout*60))  
    try:
        stdoutVal, stderrVal =  proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d" % (timeout)))
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return 214,None,None
    
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

def load_bed(bed_file, Task_list):
    # Task_list: [[chr, start, end], ...]
    bed_regions = dict()
    if bed_file != None:
        # only consider regions in BED file
        with open(bed_file, 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[0] not in bed_regions:
                    bed_regions[seq[0]] = list()
                bed_regions[seq[0]].append((int(seq[1]) - 1000, int(seq[2]) + 1000))
        region_list = [[] for i in range(len(Task_list))]
        for chrom in bed_regions:
            bed_regions[chrom].sort(key = lambda x:(x[0], x[1]))
            for item in bed_regions[chrom]:
                for i in range(len(Task_list)):
                    if chrom == Task_list[i][0]:
                        if (Task_list[i][1] <= item[0] and Task_list[i][2] > item[0]) or item[0] <= Task_list[i][1] < item[1]:
                            region_list[i].append(item)
        assert len(region_list) == len(Task_list), "parse bed file error"
        return region_list
    else:
        return None

dic_starnd = {1: '+', 2: '-'}
flag_signal = {1 << 2: 0, \
            1 >> 1: 1, \
            1 << 4: 2, \
            1 << 11: 3, \
            1 << 4 | 1 << 11: 4}
'''
    1 >> 1 means normal_foward read
    1 << 2 means unmapped read
    1 << 4 means reverse_complement read
    1 << 11 means supplementary alignment read
    1 << 4 | 1 << 11 means supplementary alignment with reverse_complement read
'''
def detect_flag(Flag):
    back_sig = flag_signal[Flag] if Flag in flag_signal else 0
    return back_sig

def analysis_bnd(ele_1, ele_2, read_name, candidate):
    '''
    *********Description*********
    *	TYPE A:		N[chr:pos[	*
    *	TYPE B:		N]chr:pos]	*
    *	TYPE C:		[chr:pos[N	*
    *	TYPE D:		]chr:pos]N	*
    *****************************
    '''
    if ele_2[0] - ele_1[1] <= 100:
        if ele_1[5] == '+':
            if ele_2[5] == '+':
                # +&+
                if ele_1[4] < ele_2[4]:
                    candidate.append(['A', 
                                        ele_1[3], 
                                        ele_2[4], 
                                        ele_2[2], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # N[chr:pos[
                else:
                    candidate.append(['D', 
                                        ele_2[2], 
                                        ele_1[4], 
                                        ele_1[3], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # ]chr:pos]N
            else:
                # +&-
                if ele_1[4] < ele_2[4]:
                    candidate.append(['B', 
                                        ele_1[3], 
                                        ele_2[4], 
                                        ele_2[3], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # N]chr:pos]
                else:
                    candidate.append(['B', 
                                        ele_2[3], 
                                        ele_1[4], 
                                        ele_1[3], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # N]chr:pos]
        else:
            if ele_2[5] == '+':
                # -&+
                if ele_1[4] < ele_2[4]:
                    candidate.append(['C', 
                                        ele_1[2], 
                                        ele_2[4], 
                                        ele_2[2], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # [chr:pos[N
                else:
                    candidate.append(['C', 
                                        ele_2[2], 
                                        ele_1[4], 
                                        ele_1[2], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # [chr:pos[N
            else:
                # -&-
                if ele_1[4] < ele_2[4]:
                    candidate.append(['D', 
                                        ele_1[2], 
                                        ele_2[4], 
                                        ele_2[3], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # ]chr:pos]N
                else:
                    candidate.append(['A', 
                                        ele_2[3], 
                                        ele_1[4], 
                                        ele_1[2], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # N[chr:pos[

def analysis_split_read(split_read, SV_size, RLength, read_name, candidate, MaxSize, query):
    '''
    read_start	read_end	ref_start	ref_end	chr	strand
    #0			#1			#2			#3		#4	#5
    '''
    SP_list = sorted(split_read, key = lambda x:x[0])

    # detect INS involoved in a translocation
    trigger_INS_TRA = 0	

    # Store Strands of INV

    if len(SP_list) == 2:
        ele_1 = SP_list[0]
        ele_2 = SP_list[1]
        if ele_1[4] == ele_2[4]:
            if ele_1[5] == ele_2[5]:
                # ins & del 
                a = 0
                if ele_1[5] == '-':
                    ele_1 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                    ele_2 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                    query = query[::-1]

                if ele_1[3] - ele_2[2] < SV_size:
                    if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                        if ele_2[2] - ele_1[3] <= 100 and (ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1] <= MaxSize or MaxSize == -1):
                            candidate.append([(ele_2[2]+ele_1[3])/2, 
                                                ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1], 
                                                read_name,
                                                str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                                "INS",
                                                ele_2[4]])
                    if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                        if ele_2[0] - ele_1[1] <= 100 and (ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3] <= MaxSize or MaxSize == -1):
                            candidate.append([ele_1[3], 
                                                ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3], 
                                                read_name,
                                                "DEL",
                                                ele_2[4]])
        else:
            trigger_INS_TRA = 1
            analysis_bnd(ele_1, ele_2, read_name, candidate)

    else:
        # over three splits
        for a in range(len(SP_list[1:-1])):
            ele_1 = SP_list[a]
            ele_2 = SP_list[a+1]
            ele_3 = SP_list[a+2]

            if ele_1[4] == ele_2[4]:
                if ele_2[4] == ele_3[4]:

                    if ele_1[5] == ele_3[5] and ele_1[5] == ele_2[5]:
                        # ins & del 
                        if ele_1[5] == '-':
                            ele_1 = [RLength-SP_list[a+2][1], RLength-SP_list[a+2][0]]+SP_list[a+2][2:]
                            ele_2 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                            ele_3 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                            query = query[::-1]

                        if ele_1[3] - ele_2[2] < SV_size:
                            if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                                if ele_2[2] - ele_1[3] <= 100 and (ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1] <= MaxSize or MaxSize == -1):
                                    if ele_3[2] >= ele_2[3]:
                                        candidate.append([(ele_2[2]+ele_1[3])/2, 
                                                            ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1], 
                                                            read_name,
                                                            str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                                            "INS",
                                                            ele_2[4]])
                            if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                                if ele_2[0] - ele_1[1] <= 100 and (ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3] <= MaxSize or MaxSize == -1):
                                    if ele_3[2] >= ele_2[3]:
                                        candidate.append([ele_1[3], 
                                                            ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3], 
                                                            read_name,
                                                            "DEL",
                                                            ele_2[4]])
                        
                        if len(SP_list) - 3 == a:
                            ele_1 = ele_2
                            ele_2 = ele_3

                            if ele_1[3] - ele_2[2] < SV_size and ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                                if ele_2[2] - ele_1[3] <= 100 and (ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1] <= MaxSize or MaxSize == -1):
                                    candidate.append([(ele_2[2]+ele_1[3])/2, 
                                                        ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1], 
                                                        read_name,
                                                        str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                                        "INS",
                                                        ele_2[4]])
                    
                            if ele_1[3] - ele_2[2] < SV_size and ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                                if ele_2[0] - ele_1[1] <= 100 and (ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3] <= MaxSize or MaxSize == -1):
                                    candidate.append([ele_1[3], 
                                                        ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3], 
                                                        read_name,
                                                        "DEL",
                                                        ele_2[4]])

            else:
                trigger_INS_TRA = 1
                analysis_bnd(ele_1, ele_2, read_name, candidate)

                if len(SP_list) - 3 == a:
                    if ele_2[4] != ele_3[4]:
                        analysis_bnd(ele_2, ele_3, read_name, candidate)

    if len(SP_list) >= 3 and trigger_INS_TRA == 1:
        if SP_list[0][4] == SP_list[-1][4]:

            if SP_list[0][5] != SP_list[-1][5]:
                pass
            else:
                if SP_list[0][5] == '+':
                    ele_1 = SP_list[0]
                    ele_2 = SP_list[-1]
                else:
                    ele_1 = [RLength-SP_list[-1][1], RLength-SP_list[-1][0]]+SP_list[-1][2:]
                    ele_2 = [RLength-SP_list[0][1],RLength-SP_list[0][0]]+SP_list[0][2:]
                    query = query[::-1]

                dis_ref = ele_2[2] - ele_1[3]
                dis_read = ele_2[0] - ele_1[1]
                if dis_ref < 100 and dis_read - dis_ref >= SV_size and (dis_read - dis_ref <= MaxSize or MaxSize == -1):

                    candidate.append([min(ele_2[2], ele_1[3]), 
                                        dis_read - dis_ref, 
                                        read_name,
                                        str(query[ele_1[1]+int(dis_ref/2):ele_2[0]-int(dis_ref/2)]),
                                        "INS",
                                        ele_2[4]])	

def acquire_clip_pos(deal_cigar):
    seq = list(cigar.Cigar(deal_cigar).items())
    if seq[0][1] == 'S':
        first_pos = seq[0][0]
    else:
        first_pos = 0
    if seq[-1][1] == 'S':
        last_pos = seq[-1][0]
    else:
        last_pos = 0

    bias = 0
    for i in seq:
        if i[1] == 'M' or i[1] == 'D' or i[1] == '=' or i[1] == 'X':
            bias += i[0]
    return [first_pos, last_pos, bias]

def organize_split_signal(primary_info, Supplementary_info, total_L, SV_size, 
    min_mapq, max_split_parts, read_name, candidate, MaxSize, query):
    split_read = list()
    if len(primary_info) > 0:
        split_read.append(primary_info)
        min_mapq = 0
    for i in Supplementary_info:
        seq = i.split(',')
        local_chr = seq[0]
        local_start = int(seq[1])
        local_cigar = seq[3]
        local_strand = seq[2]
        local_mapq = int(seq[4])
        if local_mapq >= min_mapq:
        # if local_mapq >= 0:	
            local_set = acquire_clip_pos(local_cigar)
            if local_strand == '+':
                 split_read.append([local_set[0], total_L-local_set[1], local_start, 
                     local_start+local_set[2], local_chr, local_strand])
            else:
                try:
                    split_read.append([local_set[1], total_L-local_set[0], local_start, 
                        local_start+local_set[2], local_chr, local_strand])
                except:
                    pass
    if len(split_read) <= max_split_parts or max_split_parts == -1:
        analysis_split_read(split_read, SV_size, total_L, read_name, candidate, MaxSize, query)

def generate_combine_sigs(sigs, Chr_name, read_name, svtype, candidate, merge_dis):
    # for i in sigs:
    # 	print(svtype,i, len(sigs))
    if len(sigs) == 0:
        pass
    elif len(sigs) == 1:
        if svtype == 'INS':
            candidate.append([sigs[0][0], 
                                            sigs[0][1], 
                                            read_name,
                                            sigs[0][2],
                                            svtype,
                                            Chr_name])
        else:
            candidate.append([sigs[0][0], 
                                            sigs[0][1], 
                                            read_name,
                                            svtype,
                                            Chr_name])
    else:
        temp_sig = sigs[0]
        if svtype == "INS":
            temp_sig += [sigs[0][0]]
            for i in sigs[1:]:
                if i[0] - temp_sig[3] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] += i[2]
                    temp_sig[3] = i[0]
                else:
                    candidate.append([temp_sig[0], 
                                                        temp_sig[1], 
                                                        read_name,
                                                        temp_sig[2],
                                                        svtype,
                                                        Chr_name])
                    temp_sig = i
                    temp_sig.append(i[0])
            candidate.append([temp_sig[0], 
                                                temp_sig[1], 
                                                read_name,
                                                temp_sig[2],
                                                svtype,
                                                Chr_name])
        else:
            temp_sig += [sum(sigs[0])]
            # merge_dis_bias = max([i[1]] for i in sigs)
            for i in sigs[1:]:
                if i[0] - temp_sig[2] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] = sum(i)
                else: 
                    candidate.append([temp_sig[0], 
                                                        temp_sig[1], 
                                                        read_name,
                                                        svtype,
                                                        Chr_name])
                    temp_sig = i
                    temp_sig.append(i[0])
            candidate.append([temp_sig[0], 
                                                temp_sig[1], 
                                                read_name,
                                                svtype,
                                                Chr_name])


def parse_read(read, Chr_name, SV_size, min_mapq, max_split_parts, min_read_len, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize):
    if read.query_length < min_read_len:
        return []
    candidate = list()
    Combine_sig_in_same_read_ins = list()
    Combine_sig_in_same_read_del = list()

    process_signal = detect_flag(read.flag)
    if read.mapq >= min_mapq:
        pos_start = read.reference_start # 0-based
        pos_end = read.reference_end
        shift_del = 0
        shift_ins = 0
        softclip_left = 0
        softclip_right = 0
        hardclip_left = 0
        hardclip_right = 0
        shift_ins_read = 0
        if read.cigar[0][0] == 4:
            softclip_left = read.cigar[0][1]
        if read.cigar[0][0] == 5:
            hardclip_left = read.cigar[0][1]

        for element in read.cigar:
            if element[0] in [0, 7 ,8]:
                shift_del += element[1]
            if element[0] == 2 and element[1] < min_siglength: ## changed SV_size to min_siglength
                shift_del += element[1]
            if element[0] == 2 and element[1] >= min_siglength: ## changed SV_size to min_siglength
                Combine_sig_in_same_read_del.append([pos_start+shift_del, element[1]])
                shift_del += element[1]

            # calculate offset of an ins sig in read
            if element[0] != 2:
                shift_ins_read += element[1]

            if element[0] in [0, 2, 7, 8]:
                shift_ins += element[1]
            if element[0] == 1 and element[1] >= min_siglength: ## changed SV_size to min_siglength
                Combine_sig_in_same_read_ins.append([pos_start+shift_ins, element[1],
                    str(read.query_sequence[shift_ins_read-element[1]-hardclip_left:shift_ins_read-hardclip_left])])

        
        if read.cigar[-1][0] == 4:
            softclip_right = read.cigar[-1][1]
        if read.cigar[-1][0] == 5:
            hardclip_right = read.cigar[-1][1]

        if hardclip_left != 0:
            softclip_left = hardclip_left
        if hardclip_right != 0:
            softclip_right = hardclip_right

    # ************Combine signals in same read********************
    generate_combine_sigs(Combine_sig_in_same_read_ins, Chr_name, read.query_name, "INS", candidate, merge_ins_threshold)
    generate_combine_sigs(Combine_sig_in_same_read_del, Chr_name, read.query_name, "DEL", candidate, merge_del_threshold)

    if process_signal == 1 or process_signal == 2: # 0 / 16
        Tags = read.get_tags()
        if read.mapq >= min_mapq:
            if process_signal == 1:
                primary_info = [softclip_left, read.query_length-softclip_right, pos_start, 
                pos_end, Chr_name, dic_starnd[process_signal]]
            else:
                primary_info = [softclip_right, read.query_length-softclip_left, pos_start, 
                pos_end, Chr_name, dic_starnd[process_signal]]
        else:
            primary_info = []

        for i in Tags:
            if i[0] == 'SA':
                Supplementary_info = i[1].split(';')[:-1]
                organize_split_signal(primary_info, Supplementary_info, read.query_length, 
                    SV_size, min_mapq, max_split_parts, read.query_name, candidate, MaxSize, read.query_sequence)
    return candidate

def single_pipe(sam_path, min_length, min_mapq, max_split_parts, min_read_len, temp_dir, 
                task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, bed_regions):

    candidate = list()
    reads_info_list = list()
    Chr_name = task[0]
    samfile = pysam.AlignmentFile(sam_path)

    for read in samfile.fetch(Chr_name, task[1], task[2]):
        pos_start = read.reference_start # 0-based
        pos_end = read.reference_end
        in_bed = False
        if bed_regions != None:
            for bed_region in bed_regions:
                if pos_end <= bed_region[0] or pos_start >= bed_region[1]:
                    continue
                else:
                    in_bed = True
                    break
        else:
            in_bed = True
        if read.reference_start >= task[1] and in_bed:
            read_candidate = parse_read(read, Chr_name, min_length, min_mapq, max_split_parts, 
                                    min_read_len, min_siglength, merge_del_threshold, 
                                    merge_ins_threshold, MaxSize)
            candidate.extend(read_candidate)
            if read.mapq >= min_mapq:
                is_primary = 0
                if read.flag in [0, 16]:
                    is_primary = 1
                reads_info_list.append([pos_start, pos_end, is_primary, read.query_name])
    samfile.close()
    # print('finish %s:%d-%d in %f seconds.'%(task[0], task[1], len(reads_info_list), time.time() - start_time))
   
    if len(candidate) == 0:
        logging.info("Skip %s:%d-%d."%(Chr_name, task[1], task[2]))
        return

    output = "%ssignatures/_%s_%d_%d.bed"%(temp_dir, Chr_name, task[1], task[2])
    file = open(output, 'w')
    for ele in candidate:
        if len(ele) == 5:
            file.write("%s\t%s\t%d\t%d\t%s\n"%(ele[-2], ele[-1], ele[0], ele[1], ele[2]))
        elif len(ele) == 7:
            file.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\n"%(ele[-2], ele[-1], ele[0], 
                ele[1], ele[2], ele[3], ele[4]))
        elif len(ele) == 6:
            try:
                file.write("%s\t%s\t%s\t%d\t%d\t%s\n"%(ele[-2], ele[-1], ele[0], ele[1], ele[2], ele[3]))
                # INV chr strand pos1 pos2 read_ID
            except:
                file.write("%s\t%s\t%d\t%d\t%s\t%s\n"%(ele[-2], ele[-1], ele[0], ele[1], ele[2], ele[3]))
                # INS chr pos len read_ID seq
    file.close()
    reads_output = "%ssignatures/_%s_%d_%d.reads"%(temp_dir, Chr_name, task[1], task[2])
    reads_file = open(reads_output, 'w')
    for ele in reads_info_list:
        reads_file.write("%s\t%d\t%d\t%d\t%s\n"%(Chr_name, ele[0], ele[1], ele[2], ele[3]))
    reads_file.close()
    logging.info("Finished %s:%d-%d."%(Chr_name, task[1], task[2]))	
    gc.collect()

def multi_run_wrapper(args):
    return single_pipe(*args)

def main_ctrl(args, argv):
    if not os.path.isfile(args.reference):
        raise FileNotFoundError("[Errno 2] No such file: '%s'"%args.reference)
    if not os.path.exists(args.work_dir):
        # raise FileNotFoundError("[Errno 2] No such directory: '%s'"%args.work_dir)
        os.system("mkdir -p " + args.work_dir)

    samfile = pysam.AlignmentFile(args.input)
    contig_num = len(samfile.get_index_statistics())
    logging.info("The total number of chromsomes: %d"%(contig_num))

    Task_list = list()
    chr_name_list = list()
    contigINFO = list()
    if args.work_dir[-1] == '/':
        temporary_dir = args.work_dir
    else:
        temporary_dir = args.work_dir+'/'

    ref_ = samfile.get_index_statistics()
    for i in ref_:
        chr_name_list.append(i[0])
        local_ref_len = samfile.get_reference_length(i[0])
        contigINFO.append([i[0], local_ref_len])
        if local_ref_len < args.batches:
            Task_list.append([i[0], 0, local_ref_len])
        else:
            pos = 0
            task_round = int(local_ref_len/args.batches)
            for j in range(task_round):
                Task_list.append([i[0], pos, pos+args.batches])
                pos += args.batches
            if pos < local_ref_len:
                Task_list.append([i[0], pos, local_ref_len])
    bed_regions = load_bed(args.include_bed, Task_list)
    #'''
    analysis_pools = Pool(processes=int(args.threads))
    os.mkdir("%ssignatures"%temporary_dir)
    for i in range(len(Task_list)):
        para = [(args.input, 
                    args.min_size, 
                    args.min_mapq, 
                    args.max_split_parts, 
                    args.min_read_len, 
                    temporary_dir, 
                    Task_list[i], 
                    args.min_siglength, 
                    args.merge_del_threshold, 
                    args.merge_ins_threshold, 
                    args.max_size,
                    None if bed_regions == None else bed_regions[i])]
        analysis_pools.map_async(multi_run_wrapper, para)
    analysis_pools.close()
    analysis_pools.join()
    #'''
    #'''
    logging.info("Rebuilding signatures of structural variants.")
    analysis_pools = Pool(processes=int(args.threads))
    cmd_del = ("cat %ssignatures/*.bed | grep -w DEL | sort -u -T %s | sort -k 2,2 -k 3,3n -T %s > %sDEL.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    cmd_ins = ("cat %ssignatures/*.bed | grep -w INS | sort -u -T %s | sort -k 2,2 -k 3,3n -T %s > %sINS.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    cmd_reads = ("cat %ssignatures/*.reads > %sreads.sigs"%(temporary_dir, temporary_dir))
    # for i in [cmd_ins, cmd_del, cmd_dup, cmd_tra, cmd_inv, cmd_reads]:
    for i in [cmd_ins, cmd_del, cmd_reads]:
        analysis_pools.map_async(exe, (i,))
    analysis_pools.close()
    analysis_pools.join()

def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))


def run(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
    main_ctrl(args, argv)
    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

def parseArgs(argv):
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter)

	# **************Parameters of input******************
	parser.add_argument("input", 
		metavar="[BAM]", 
		type = str, 
		help ="Sorted .bam file from NGMLR or Minimap2.")
	parser.add_argument("reference",  
		type = str, 
		help ="The reference genome in fasta format.")
	parser.add_argument('work_dir', 
		type = str, 
		help = "Work-directory for distributed jobs")

	# ************** Other Parameters******************
	parser.add_argument('-t', '--threads', 
		help = "Number of threads to use.[%(default)s]", 
		default = 16, 
		type = int)
	parser.add_argument('-b', '--batches', 
		help = "Batch of genome segmentation interval.[%(default)s]", 
		default = 10000000, 
		type = int)
	# The description of batches needs to improve.
	parser.add_argument('-S', '--sample',
		help = "Sample name/id",
		default = "NULL",
		type = str)

	parser.add_argument('--report_readid',
		help = "Enable to report supporting read ids for each SV.",
		action="store_true")

	# **************Parameters in signatures collection******************
	GroupSignaturesCollect = parser.add_argument_group('Collection of SV signatures')
	GroupSignaturesCollect.add_argument('-p', '--max_split_parts', 
		help = "Maximum number of split segments a read may be aligned before it is ignored. All split segments are considered when using -1. \
			(Recommand -1 when applying assembly-based alignment.)[%(default)s]", 
		default = 7, 
		type = int)
	GroupSignaturesCollect.add_argument('-q', '--min_mapq', 
		help = "Minimum mapping quality value of alignment to be taken into account.[%(default)s]", 
		default = 20, 
		type = int)
	GroupSignaturesCollect.add_argument('-r', '--min_read_len', 
		help = "Ignores reads that only report alignments with not longer than bp.[%(default)s]", 
		default = 500, 
		type = int)
	GroupSignaturesCollect.add_argument('-md', '--merge_del_threshold', 
		help = "Maximum distance of deletion signals to be merged.",
		default = 0, 
		type = int)
	GroupSignaturesCollect.add_argument('-mi', '--merge_ins_threshold', 
		help = "Maximum distance of insertion signals to be merged.",
		default = 100, 
		type = int)
	GroupSignaturesCollect.add_argument('-include_bed', 
		help = "Optional given bed file. Only detect SVs in regions in the BED file. [NULL]",
		default = None,
        type = str)
	# The min_read_len in last version is 2000.
	# signatures with overlap need to be filtered

	# **************Parameters in clustering******************
	GroupSVCluster = parser.add_argument_group('Generation of SV clusters')
	GroupSVCluster.add_argument('-s', '--min_support', 
		help = "Minimum number of reads that support a SV to be reported.[%(default)s]", 
		default = 10, 
		type = int)
	GroupSVCluster.add_argument('-l', '--min_size', 
		help = "Minimum size of SV to be reported.[%(default)s]", 
		default = 30, 
		type = int)
	GroupSVCluster.add_argument('-L', '--max_size', 
		help = "Maximum size of SV to be reported. All SVs are reported when using -1. [%(default)s]", 
		default = 100000, 
		type = int)
	GroupSVCluster.add_argument('-sl', '--min_siglength', 
		help = "Minimum length of SV signal to be extracted.[%(default)s]", 
		default = 10, 
		type = int)

	# **************Parameters in genotyping******************
	GroupGenotype = parser.add_argument_group('Computing genotypes')
	GroupGenotype.add_argument('--genotype',
		help = "Enable to generate genotypes.",
		action="store_true")
	GroupGenotype.add_argument('--gt_round', 
		help = "Maximum round of iteration for alignments searching if perform genotyping.[%(default)s]", 
		default = 500, 
		type = int)

	# **************Advanced Parameters******************
	GroupAdvanced = parser.add_argument_group('Advanced')

	# ++++++INS++++++
	GroupAdvanced.add_argument('--max_cluster_bias_INS', 
		help = "Maximum distance to cluster read together for insertion.[%(default)s]", 
		default = 100, 
		type = int)
	GroupAdvanced.add_argument('--diff_ratio_merging_INS', 
		help = "Do not merge breakpoints with basepair identity more than [%(default)s] for insertion.", 
		default = 0.3, 
		type = float)

	# ++++++DEL++++++
	GroupAdvanced.add_argument('--max_cluster_bias_DEL', 
		help = "Maximum distance to cluster read together for deletion.[%(default)s]", 
		default = 200, 
		type = int)
	GroupAdvanced.add_argument('--diff_ratio_merging_DEL', 
		help = "Do not merge breakpoints with basepair identity more than [%(default)s] for deletion.", 
		default = 0.5, 
		type = float)

	GroupAdvanced.add_argument('--remain_reads_ratio', 
		help = "The ratio of reads remained in cluster. Set lower when the alignment data have high quality but recommand over 0.5.[%(default)s]", 
		default = 1.0, 
		type = float)

	args = parser.parse_args(argv)
	return args

if __name__ == '__main__':
    run(sys.argv[1:])
