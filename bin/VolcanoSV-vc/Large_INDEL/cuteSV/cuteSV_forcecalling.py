from cuteSV_genotype import cal_CIPOS, overlap_cover, assign_gt
from multiprocessing import Pool
from pysam import VariantFile
import math
import time
import logging
import numpy as np

def parse_svtype(sv_type):
    if 'DEL' in sv_type:
        return 'DEL'
    if 'INS' in sv_type:
        return 'INS'
    if 'INV' in sv_type:
        return 'INV'
    if 'DUP' in sv_type:
        return 'DUP'
    if 'BND' in sv_type or 'TRA' in sv_type:
        return 'TRA'
    return 'NA'

def parse_to_int(sth):
    if sth == None:
        return 0
    elif isinstance(sth, str):
        return int(sth)
    elif isinstance(sth, list):
        return parse_to_int(sth[0])
    elif isinstance(sth, tuple):
        return parse_to_int(sth[0])
    elif isinstance(sth, int):
        return sth
    else:
        return sth

# INDEL end->len, others:end
def parse_record(record):
    sv_type = parse_svtype(record.info['SVTYPE'])
    chrom1 = record.chrom
    start = parse_to_int(record.pos)
    chrom2 = ''
    svlen = 0
    end = 0 # len for indel, end for others
    if 'SVLEN' in record.info:
        svlen = abs(parse_to_int(record.info['SVLEN']))
    if 'END' in record.info:
        end = parse_to_int(record.info['END'])
    try:
        end = parse_to_int(record.stop)
    except:
        pass
    if sv_type == 'INS' or sv_type == 'DEL':
        end = svlen
    else:
        if start == end:
            end = start + svlen
    if sv_type == 'TRA':
        if 'CHR2' in record.info:
            chrom2 = record.info['CHR2']
        if 'END' in record.info:
            end = parse_to_int(record.info['END'])
        try:
            tra_alt = str(record.alts[0])
            if tra_alt[0] == 'N':
                if tra_alt[1] == '[':
                    tra_type = 'A'
                else:
                    tra_type = 'B'
            elif tra_alt[0] == '[':
                tra_type = 'C'
            else:
                tra_type = 'D'
            if tra_alt[0] == 'N':
                tra_alt = tra_alt[2:-1]
            else:
                tra_alt = tra_alt[1:-2]
            if ':' in tra_alt:
                chrom2 = tra_alt.split(':')[0]
                end = int(tra_alt.split(':')[1])
        except:
            pass
    strand = '.'
    if sv_type != 'TRA':
        chrom2 = record.chrom
        if 'STRAND' in record.info:
            strand = record.info['STRAND']
        elif 'STRANDS' in record.info:
            strand = record.info['STRANDS']
    if isinstance(strand, tuple) or isinstance(strand, list):
        strand = strand[0]
    svid = record.id
    ref = record.ref
    alts = record.alts[0]
    if 'SEQ' in record.info:
        if record.info['SVTYPE'] == 'INS' and record.alts[0] == '<INS>':
            alts = record.info['SEQ']
        if record.info['SVTYPE'] == 'DEL' and record.alts[0] == '<DEL>':
            ref = record.info['SEQ']
    return sv_type, chrom1, chrom2, start, end, strand, svid, ref, alts

def parse_sigs_chrom(var_type, work_dir, chrom_list):
    if var_type == 'DEL' or var_type == 'DUP':
        var_dict = dict()  #var_dict[chrom] = [chrom, start, len/end, read_id]
        with open(work_dir + var_type + '.sigs', 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[1] in chrom_list:
                    if seq[1] not in var_dict:
                        var_dict[seq[1]] = []
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
        return var_dict
    if var_type == 'INS':
        var_dict = dict()  #var_dict[chrom] = [chrom, start, len, read_id, seq]
        with open(work_dir + 'INS.sigs', 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[1] in chrom_list:
                    if len(seq) < 6:
                        cigar = '<INS>'
                    else:
                        cigar = seq[5]
                    cigar = '<INS>'
                    if seq[1] not in var_dict:
                        var_dict[seq[1]] = []
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4], cigar])
        return var_dict
    if var_type == 'INV':
        var_dict = dict() # var_dict[chrom] = [chrom, start, end, read_id]
        with open(work_dir + 'INV.sigs', 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[1] in chrom_list:
                    chrom = seq[1]
                    if chrom not in var_dict:
                        var_dict[chrom] = []
                    var_dict[chrom].append([chrom, int(seq[3]), int(seq[4]), seq[5]])
        for chrom in var_dict:
            var_dict[chrom].sort(key=lambda x:x[1])
        return var_dict
    if var_type == 'TRA':
        var_dict = dict()  #var_dict[chrom1][chrom2] = [[chrom2, pos1, pos2, read_id]]
        with open(work_dir + 'TRA.sigs', 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[1] not in chrom_list:
                    continue
                chrom1 = seq[1]
                tra_type = seq[2]
                pos1 = int(seq[3])
                chrom2 = seq[4]
                pos2 = int(seq[5])
                read_id = seq[6]
                '''
                if chrom1 in var_dict:
                    if tra_type in var_dict[chrom1]:
                        if chrom2 in var_dict[chrom1][tra_type]:
                            var_dict[chrom1][tra_type][chrom2].append([chrom2, pos1, pos2, read_id])
                        else:
                            var_dict[chrom1][tra_type][chrom2] = [[chrom2, pos1, pos2, read_id]]
                    else:
                        var_dict[chrom1][tra_type] = dict()
                        var_dict[chrom1][tra_type][chrom2] = [[chrom2, pos1, pos2, read_id]]

                else:
                    var_dict[chrom1] = dict()
                    var_dict[chrom1][tra_type] = dict()
                    var_dict[chrom1][tra_type][chrom2] = [[chrom2, pos1, pos2, read_id]]
                '''
                if chrom1 not in var_dict:
                    var_dict[chrom1] = dict()
                if chrom2 not in var_dict[chrom1]:
                    var_dict[chrom1][chrom2] = []
                var_dict[chrom1][chrom2].append([chrom2, pos1, pos2, read_id])
        for chr1 in var_dict:
            for chr2 in var_dict[chr1]:
                var_dict[chr1][chr2].sort(key=lambda x:x[1])
        return var_dict

def check_same_variant(sv_type, end1, end2):
    if sv_type == 'INS' or sv_type == 'DEL':
        return 0.7 < min(end1, end2) / max(end1, end2) <= 1
    return abs(end1 - end2) < 1000

# var_list read_id which is similar to (pos, sv_end), the reads support the variant
def find_in_list(var_type, var_list, bias, pos, sv_end):
    '''
    if interval_start == 44058795:
        print('start find in list')
        print('interval_start=%d, interval_end=%d, pos=%d, sv_end=%d'%(interval_start, interval_end, pos, sv_end))
    '''
    if len(var_list) == 0:
        return [], 0
    left = 0
    right = len(var_list) - 1
    mid = 0
    while left < right:
        mid = (left + right) >> 1
        if var_list[mid][1] < pos:
            left = mid + 1
        else:
            right = mid
    read_id_list = set()
    search_start = -1
    search_end = -1
    if right > 0 and pos - var_list[right - 1][1] <= bias: ##
        for i in range(right - 1, -1, -1):
            if check_same_variant(var_type, var_list[i][2], sv_end):
                read_id_list.add(var_list[i][3])
                search_start = var_list[i][1]
            if i > 0 and (var_list[i][1] - var_list[i - 1][1] > bias or pos - var_list[i - 1][1] > 2000):
                break
    if var_list[right][1] - pos <= bias:
        for i in range(right, len(var_list)):
            if check_same_variant(var_type, var_list[i][2], sv_end):  # if abs(var_list[i][2] - sv_end) < 1000:
                read_id_list.add(var_list[i][3])
                search_end = var_list[i][1]
            if i < len(var_list) - 1 and (var_list[i + 1][1] - var_list[i][1] > bias or var_list[i + 1][1] - pos > 2000):
                break
    if search_start == -1:
        search_start = pos
    if search_end == -1:
        search_end = pos
    search_threshold = max(abs(pos - search_start), abs(pos - search_end))
    return list(read_id_list), search_threshold

def compare_len(len1, len2): # len1 < len2
    if len2 < 100:
        if len1 / len2 > 0.6:
            return True
    else:
        if len1 / len2 > 0.8:
            return True
    return False

def find_in_indel_list(var_type, var_list, bias_origin, pos, sv_end, threshold_gloab):
    # bias = min(bias_origin * 10, 2000)
    bias = bias_origin
    debug_pos = -1
    if len(var_list) == 0:
        return [], 0, '.,.', '.,.'
    left = 0
    right = len(var_list) - 1
    mid = 0
    while left < right:
        mid = (left + right) >> 1
        if var_list[mid][1] < pos:
            left = mid + 1
        else:
            right = mid
    candidates = []
    if right > 0 and pos - var_list[right - 1][1] <= bias:
        for i in range(right - 1, -1, -1):
            candidates.append(var_list[i]) # [chrom, start, len, read_id] or [chrom, start, len, read_id, seq]
            if i > 0 and (var_list[i][1] - var_list[i - 1][1] > bias or pos - var_list[i - 1][1] > 2 * bias):
                break
    if var_list[right][1] - pos <= bias:
        for i in range(right, len(var_list)):
            candidates.append(var_list[i])
            if i < len(var_list) - 1 and (var_list[i + 1][1] - var_list[i][1] > bias or var_list[i + 1][1] - pos > 2 * bias):
                break
    if len(candidates) == 0:
        return [], 0, '.,.', '.,.'
    read_tag = dict()
    for element in candidates:
        if element[3] not in read_tag:
            read_tag[element[3]] = []
        read_tag[element[3]].append(element)

    # merge sigs on the same read
    read_tag2SortedList = []
    for read_id in read_tag:
        for i in range(len(read_tag[read_id])):
            read_tag2SortedList.append(read_tag[read_id][i])
            for j in range(i + 1, len(read_tag[read_id]), 1):
                if var_type == 'DEL':
                    read_tag2SortedList.append([read_tag[read_id][i][0], int((read_tag[read_id][i][1]+read_tag[read_id][j][1])/2), read_tag[read_id][i][2]+read_tag[read_id][j][2], read_tag[read_id][i][3] ])
                else:
                    read_tag2SortedList.append([read_tag[read_id][i][0], int((read_tag[read_id][i][1]+read_tag[read_id][j][1])/2), read_tag[read_id][i][2]+read_tag[read_id][j][2], read_tag[read_id][i][3], read_tag[read_id][i][4] ])
                for k in range(j + 1, len(read_tag[read_id]), 1):
                    if var_type == 'DEL':
                        read_tag2SortedList.append([read_tag[read_id][i][0], int((read_tag[read_id][i][1]+read_tag[read_id][j][1]+read_tag[read_id][k][1])/3), read_tag[read_id][i][2]+read_tag[read_id][j][2]+read_tag[read_id][k][2], read_tag[read_id][i][3] ])
                    else:
                        read_tag2SortedList.append([read_tag[read_id][i][0], int((read_tag[read_id][i][1]+read_tag[read_id][j][1]+read_tag[read_id][k][1])/3), read_tag[read_id][i][2]+read_tag[read_id][j][2]+read_tag[read_id][k][2], read_tag[read_id][i][3], read_tag[read_id][i][4] ])
    
    read_tag2SortedList = sorted(read_tag2SortedList, key = lambda x:x[2])
    # print(read_tag2SortedList)
    global_len = [i[2] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_TEMP = threshold_gloab * np.mean(global_len)
    if var_type == 'DEL':
        last_len = read_tag2SortedList[0][2]
        cur_bias = last_len * threshold_gloab
        allele_collect = list()
        allele_collect.append([[read_tag2SortedList[0][1]],  # start
                                [read_tag2SortedList[0][2]],  # len
                                [],  # support
                                [read_tag2SortedList[0][3]]])  # read_id
        for i in read_tag2SortedList[1:]:
            if i[2] - last_len > cur_bias:
                allele_collect[-1][2].append(len(allele_collect[-1][0]))
                allele_collect.append([[],[],[],[]])
            allele_collect[-1][0].append(i[1])
            allele_collect[-1][1].append(i[2])
            allele_collect[-1][3].append(i[3])
            last_len = (last_len * (len(allele_collect[-1][0]) - 1) + i[2]) / len(allele_collect[-1][0])
            cur_bias = last_len * threshold_gloab

        allele_collect[-1][2].append(len(allele_collect[-1][0]))
        allele_idx = -1
        nearest_gap = 0x3f3f3f3f
        for i in range(len(allele_collect)):
            allele = allele_collect[i]
            signalLen = np.mean(allele[1])
            if min(signalLen, sv_end) / max(signalLen, sv_end) > 0.7 and abs(signalLen - sv_end) < nearest_gap:
                allele_idx = i
                nearest_gap = abs(signalLen - sv_end)
        
        if allele_idx == -1:
            read_id_set = set()
            CIPOS = "-0,0"
            CILEN = "-0,0"
            seq = '<DEL>'
            search_threshold = 0
        else:
            final_alleles = [[],[],[],[]]
            for i in range(len(allele_collect[allele_idx][0])):
                if min(allele_collect[allele_idx][1][i], sv_end) / max(allele_collect[allele_idx][1][i], sv_end) > 0:
                    final_alleles[0].append(allele_collect[allele_idx][0][i])
                    final_alleles[1].append(allele_collect[allele_idx][1][i])
                    final_alleles[3].append(allele_collect[allele_idx][3][i])
            if len(final_alleles[0]) == 0:
                read_id_set = set()
                CIPOS = "-0,0"
                CILEN = "-0,0"
                seq = '<DEL>'
                search_threshold = 0
            else:
                read_id_set = set(final_alleles[3])
                CIPOS = cal_CIPOS(np.std(final_alleles[0]), len(final_alleles[0]))
                CILEN = cal_CIPOS(np.std(final_alleles[1]), len(final_alleles[1]))
                seq = '<DEL>'
                search_start = min(final_alleles[0])
                search_end = max(final_alleles[0])
                search_threshold = min(abs(pos - search_start), abs(pos - search_end))
    if var_type == 'INS':
        last_len = read_tag2SortedList[0][2]
        cur_bias = last_len * threshold_gloab
        allele_collect = list()
        allele_collect.append([[read_tag2SortedList[0][1]],  # start
                                [read_tag2SortedList[0][2]],  # len
                                [],  # support
                                [read_tag2SortedList[0][3]],  # read_id
                                [read_tag2SortedList[0][4]]])  # ins_seq
        for i in read_tag2SortedList[1:]:
            if i[2] - last_len > cur_bias:
                allele_collect[-1][2].append(len(allele_collect[-1][0]))
                allele_collect.append([[],[],[],[],[]])
            allele_collect[-1][0].append(i[1])
            allele_collect[-1][1].append(i[2])
            allele_collect[-1][3].append(i[3])
            allele_collect[-1][4].append(i[4])
            last_len = (last_len * (len(allele_collect[-1][0]) - 1) + i[2]) / len(allele_collect[-1][0])
            cur_bias = last_len * threshold_gloab
        allele_collect[-1][2].append(len(allele_collect[-1][0]))
        allele_idx = -1
        nearest_gap = 0x3f3f3f3f
        for i in range(len(allele_collect)):
            allele = allele_collect[i]
            signalLen = np.mean(allele[1])
            if min(signalLen, sv_end) / max(signalLen, sv_end) > 0.7 and abs(signalLen - sv_end) < nearest_gap:
                allele_idx = i
                nearest_gap = abs(signalLen - sv_end)
 
        if allele_idx == -1:
            read_id_set = set()
            CIPOS = "-0,0"
            CILEN = "-0,0"
            seq = '<INS>'
            search_threshold = 0
        else:
            final_alleles = [[],[],[],[],[]]
            for i in range(len(allele_collect[allele_idx][0])):
                if min(allele_collect[allele_idx][1][i], sv_end) / max(allele_collect[allele_idx][1][i], sv_end) > 0:
                    final_alleles[0].append(allele_collect[allele_idx][0][i])
                    final_alleles[1].append(allele_collect[allele_idx][1][i])
                    final_alleles[3].append(allele_collect[allele_idx][3][i])
                    final_alleles[4].append(allele_collect[allele_idx][4][i])
            if len(final_alleles[0]) == 0:
                read_id_set = set()
                CIPOS = "-0,0"
                CILEN = "-0,0"
                seq = '<DEL>'
                search_threshold = 0
            else:
                read_id_set = set(final_alleles[3])
                CIPOS = cal_CIPOS(np.std(final_alleles[0]), len(final_alleles[0]))
                CILEN = cal_CIPOS(np.std(final_alleles[1]), len(final_alleles[1]))
                seq = '<DEL>'
                search_start = min(final_alleles[0])
                search_end = max(final_alleles[0])
                search_threshold = min(abs(pos - search_start), abs(pos - search_end))
    return list(read_id_set), search_threshold, CIPOS, CILEN

def generate_dispatch(reads_count, chrom_list):
    dispatch = [[]]
    cur_count = 0
    for item in reads_count:
        if cur_count >= 10000:
            dispatch.append([])
            cur_count = 0
        cur_count = item[1]
        dispatch[-1].append(item[0])
    # chrom without reads
    for chrom in chrom_list:
        flag = 0
        for x in reads_count:
            if chrom == x[0]:
                flag = 1
        if flag == 0:
            dispatch[0].append(chrom)
    return dispatch

def force_calling_chrom(ivcf_path, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round, threads):
    logging.info('Check the parameter -Ivcf: OK.')
    logging.info('Enable to perform force calling.')

    # parse svs tobe genotyped
    vcf_reader = VariantFile(ivcf_path, 'r')
    svs_tobe_genotyped = dict()
    for record in vcf_reader.fetch():
        sv_type, chrom, sv_chr2, pos, sv_end, sv_strand, svid, ref, alts = parse_record(record)
        if sv_type not in ["DEL", "INS", "DUP", "INV", "TRA"]:
            continue
        if chrom not in svs_tobe_genotyped:
            svs_tobe_genotyped[chrom] = list()
        svs_tobe_genotyped[chrom].append([sv_type, sv_chr2, pos, sv_end, svid, ref, alts, sv_strand, chrom])
    
    # parse reads in alignment
    reads_count = dict()
    with open('%sreads.sigs'%(temporary_dir), 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[0] not in reads_count:
                reads_count[seq[0]] = 0
            reads_count[seq[0]] += 1
    reads_count = sorted(reads_count.items(), key=lambda x:x[1])
    dispatch = generate_dispatch(reads_count, svs_tobe_genotyped.keys())
    
    # force calling
    pool_result = list()
    result = list()
    process_pool = Pool(processes = threads)
    # dispatch = [['MT']]
    for chroms in dispatch:
        genotype_sv_list = dict()
        for chrom in chroms:
            if chrom in svs_tobe_genotyped:
                genotype_sv_list[chrom] = svs_tobe_genotyped[chrom]
        if len(genotype_sv_list) == 0:
            continue
        # pool_result.append(solve_fc(chroms, genotype_sv_list, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round))
        fx_para = [(chroms, genotype_sv_list, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round)]
        pool_result.append(process_pool.map_async(solve_fc_wrapper, fx_para))
    process_pool.close()
    process_pool.join()

    for x in pool_result:
        result.extend(x.get()[0])
    return result

def solve_fc_wrapper(args):
    return solve_fc(*args)
def solve_fc(chrom_list, svs_dict, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round):
    reads_info = dict() # [10000, 10468, 0, 'm54238_180901_011437/52298335/ccs']
    readsfile = open("%sreads.sigs"%(temporary_dir), 'r')
    for line in readsfile:
        seq = line.strip().split('\t')
        chr = seq[0]
        if chr not in chrom_list:
            continue
        if chr not in reads_info:
            reads_info[chr] = list()
        reads_info[chr].append([int(seq[1]), int(seq[2]), int(seq[3]), seq[4]])
    readsfile.close()
    
    sv_dict = dict()
    for sv_type in ["DEL", "DUP", "INS", "INV", "TRA"]:
        sv_dict[sv_type] = parse_sigs_chrom(sv_type, temporary_dir, chrom_list)
    
    gt_list = list()
    for chrom in svs_dict:
        read_id_dict = dict()
        ci_dict = dict()
        search_list = list()
        for i in range(len(svs_dict[chrom])):
            record = svs_dict[chrom][i]
            sv_type = record[0]
            sv_chr2 = record[1]
            sv_start = record[2]
            sv_end = record[3]
            chrom = record[8]
            search_id_list = list()
            # rewrite!
            if sv_type == 'TRA' and 'TRA' in sv_dict and chrom in sv_dict['TRA'] and sv_chr2 in sv_dict['TRA'][chrom]:
                search_id_list = sv_dict['TRA'][chrom][sv_chr2]
            elif sv_type != 'TRA' and sv_type in sv_dict and chrom in sv_dict[sv_type]:
                search_id_list = sv_dict[sv_type][chrom]
            max_cluster_bias = 0
            if sv_type == 'INS' or sv_type == 'DEL':
                sigs_bias = max_cluster_bias_dict[sv_type]
                for j in range(i+1, len(svs_dict[chrom]), 1):
                    if svs_dict[chrom][j][0] == sv_type:
                        sigs_bias = min(sigs_bias, max(100, svs_dict[chrom][j][2] - sv_start))
                        break
                read_id_list, max_cluster_bias, CIPOS, CILEN = find_in_indel_list(sv_type, search_id_list, sigs_bias, sv_start, sv_end, threshold_gloab_dict[sv_type])
            else:
                read_id_list, max_cluster_bias = find_in_list(sv_type, search_id_list, max_cluster_bias_dict[sv_type], sv_start, sv_end)
                CIPOS = '.'
                CILEN = '.'
            
            # max_cluster_bias = max(max_cluster_bias_dict[sv_type], max_cluster_bias)
            max_cluster_bias = max(1000, max_cluster_bias)
            if sv_type == 'INS' or sv_type == 'DEL' or sv_type == 'TRA':
                search_list.append((max(sv_start - max_cluster_bias, 0), sv_start + max_cluster_bias))
            elif sv_type == 'INV' or sv_type == 'DUP':
                # search_list.append((max(sv_start - max_cluster_bias/2, 0), sv_start + max_cluster_bias/2))
                search_list.append((sv_start, sv_end))

            read_id_dict[i] = read_id_list
            ci_dict[i] = (CIPOS, CILEN)

        if chrom in reads_info:
            iteration_dict, primary_num_dict, cover_dict = overlap_cover(search_list, reads_info[chrom]) # both key(sv idx), value(set(read id))
        else:
            iteration_dict = dict()
            primary_num_dict = dict()
            cover_dict = dict()
            for i in read_id_dict:
                iteration_dict[i] = 0
                primary_num_dict[i] = 0
                cover_dict[i] = set()
        assert len(iteration_dict) == len(read_id_dict), "overlap length error"
        assert len(cover_dict) == len(svs_dict[chrom]), "cover length error"
        assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict)
        for i in range(len(svs_dict[chrom])):
            assert len(assign_list[i]) == 6, "assign genotype error"
            record = svs_dict[chrom][i]
            rname = ','.join(read_id_dict[i])
            if rname == '':
                rname = 'NULL'
            if record[6] == '<TRA>' or record[6] == '<BND>':
                seq = str(record[1]) + ':' + str(record[3])
            else:
                seq = '<' + record[0] + '>'
            gt_list.append([record[8], record[2], assign_list[i][2], record[0], record[3],
                            ci_dict[i][0], ci_dict[i][1], assign_list[i], rname, record[4],
                            record[5], record[6],
                            record[7], seq])
        logging.info("Finished calling %s."%(chrom))
    return gt_list
