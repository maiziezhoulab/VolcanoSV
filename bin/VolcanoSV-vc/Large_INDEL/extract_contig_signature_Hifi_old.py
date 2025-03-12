import matplotlib.pyplot as plt
import numpy as np
import pickle
from tqdm import tqdm
import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import logging
import os
from argparse import ArgumentParser
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_number','-chr',type = int)
parser.add_argument('--bam_path','-bam')
parser.add_argument('--contig_path','-contig')
parser.add_argument('--header_path','-header')
parser.add_argument('--ref_path','-ref')
parser.add_argument('--output_dir','-o')

parser.add_argument('--max_shift',type = int, default = 100)
parser.add_argument('--max_shift_ratio',type = float, default = 0.1)
parser.add_argument('--min_reads_support',type = int, default = 1)
parser.add_argument('--min_siglen',type =int, default = 30)
parser.add_argument('--min_cigar_mapq',type = int, default = 50)
parser.add_argument('--min_split_mapq',type = int, default = 50 )

args = parser.parse_args()
chr_number = args.chr_number
bam_path = args.bam_path
header_path = args.header_path
ref_path = args.ref_path
contig_path = args.contig_path
output_dir = args.output_dir


max_shift = args.max_shift
max_shift_ratio = args.max_shift_ratio
min_reads_support = args.min_reads_support
min_siglen = args.min_siglen
min_cigar_mapq = args.min_cigar_mapq
min_split_mapq = args.min_cigar_mapq



# chr_name = "chr%d"%chr_number
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

def extract_sig_from_cigar(read,min_svlen):
    ref_name = read.reference_name
    start = read.pos
    cigar = read.cigar
    qname = read.qname
    
    offset_ref = start
    del_sig = []
    offset_contig = 0
    ins_sig = []
    hard_clip_head = 0
    if cigar[0][0]==5:
        hard_clip_head = cigar[0][1]
        
    if read.is_reverse:
        direction = '-'
    else:
        direction = '+'
        
    for tp in cigar:
        if tp[0] == 0:
            offset_ref+=tp[1]
            offset_contig +=tp[1]
        elif tp[0]==4:
            offset_contig +=tp[1]
        elif tp[0]==2:
            if tp[1]>=min_svlen:
                del_sig.append([ref_name,'DEL',offset_ref,tp[1],qname,offset_contig + hard_clip_head,offset_contig + hard_clip_head+1,direction,'cigar',read.mapq])  
            offset_ref+=tp[1]
        elif tp[0] == 1:   
            if tp[1]>=min_svlen:
                ins_sig.append([ref_name,'INS',offset_ref, tp[1],qname,offset_contig + hard_clip_head,offset_contig + hard_clip_head+tp[1],direction,'cigar',read.mapq])
            offset_contig +=tp[1]


    ### cluster ins_sig
    cluster_list = [0]*len(ins_sig)

    def merge_two_ins(sig1,sig2):
        pos = sig1[2]
        read_start = sig1[5]
        read_end = sig2[6]
        direction = sig1[7]
        merged_sig = [sig1[0],sig1[1],pos,read_end-read_start, sig1[4],read_start,read_end,direction,'cigar',read.mapq]
        return merged_sig

    def merge_two_del(sig1,sig2):
        pos = sig1[2]
        read_start = sig1[5]
        read_end = sig1[5]+1
        direction = sig1[7]
        svlen = sig2[2]+sig2[3]-sig1[2]
        merged_sig = [sig1[0],sig1[1],pos,svlen, sig1[4],read_start,read_end,direction,'cigar',read.mapq]
        return merged_sig

    def cluster_ins_one_read(ins_list):
        if len(ins_list)>=2:
            result_list = [ins_list[0]]

            for i in range(1,len(ins_list)):
                sig1 = result_list[-1]
                sig2 = ins_list[i]

                if (sig1[3]>250) and \
                (sig2[3]>250) and \
                (abs(sig2[2]-sig1[2])<250):
                    result_list[-1] = merge_two_ins(sig1,sig2)

                elif (sig1[3]>320) and \
                (sig2[3]>320) and \
                (abs(sig2[2]-sig1[2])<380):
                    result_list[-1] = merge_two_ins(sig1,sig2)


                elif (sig1[3]>100) and \
                (sig2[3]>100) and \
                (abs(sig2[2]-sig1[2])<250):
                # if ((sig1[3]<1200) and (sig2[3]<1200)) and (abs(sig2[2]-sig1[2])<200):
                #     result_list[-1] = merge_two_ins(sig1,sig2)
                # elif ((sig1[3]>=1200) and (sig2[3]>=1200)) and (abs(sig2[2]-sig1[2])<500):
                    result_list[-1] = merge_two_ins(sig1,sig2)
                else:
                    result_list.append(sig2)
            return result_list
        else:
            return ins_list

    def cluster_del_one_read(del_list):
        if len(del_list)>=2:
            result_list = [del_list[0]]

            for i in range(1,len(del_list)):
                sig1 = result_list[-1]
                sig2 = del_list[i]
                sig1_end = sig1[2]+sig1[3]
                sig2_start = sig2[2]

                if (sig1[3]>150) and \
                (sig2[3]>150) and \
                (abs(sig2[2]-sig1[2])<150):
                # if ((sig1[3]<1200) and (sig2[3]<1200)) and (abs(sig2[2]-sig1[2])<200):
                #     result_list[-1] = merge_two_del(sig1,sig2)
                # elif ((sig1[3]>=1200) and (sig2[3]>=1200)) and (abs(sig2[2]-sig1[2])<500):
                    result_list[-1] = merge_two_del(sig1,sig2)
                else:
                    result_list.append(sig2)
            return result_list
        else:
            return del_list

    ins_sig = cluster_ins_one_read(ins_sig)
    del_sig = cluster_del_one_read(del_sig)

    return (del_sig,ins_sig,offset_ref,offset_contig)



def sort_sig(sig_list):
    pos_list = []
    for sig in sig_list:
        pos_list.append(sig[2])
    idxs = np.argsort(pos_list)
    sorted_sig_list = []
    
    for idx in idxs:
        sorted_sig_list.append(sig_list[idx])
    return sorted_sig_list

def write_sig_cigar(sig_list,output_path):
#     sig_list = sig_list1.copy()
    with open(output_path,'w') as f:
        for sig1 in sig_list:
            sig = sig1.copy()
            sig[2]=str(sig[2])
            sig[3]=str(sig[3])
            sig[5]=str(sig[5])
            sig[6] =str(sig[6])
            sig[9] =str(sig[9])
            # print(sig)
            line = '\t'.join(sig)+'\n'
            f.write(line)
    return

def cluster_del(sig_list,max_shift = 100, min_overlap_ratio = 0.5, min_size_similarity = 0.5):
    
    cluster = [-1]*len(sig_list)
#     print(sig_list[:10])
    def calculate_overlap_ratio(start1,end1,start2,end2):
        len1 = end1-start1
        len2 = end2-start2
        minlen = min(len1,len2)
        return (min(end1,end2)-max(start1,start2))/minlen
    
    while (np.array(cluster)==-1).sum()!=0:
        for i in range(len(sig_list)):
            if cluster[i]==-1:
                sig1 = sig_list[i]
                cluster[i] = i
                for j in range(len(sig_list)):
                    if cluster[j]==-1:
                        sig2 = sig_list[j]
                        start1 = sig1[2]
                        start2 = sig2[2]
                        end1 = start1+sig1[3]
                        end2 = start2+sig2[3]
#                         print(type(start1),end1)
                        overlap_ratio = calculate_overlap_ratio(start1,end1,start2,end2)
                        size_similarity = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
                        shift = abs(sig1[2]-sig2[2])
                        
                        if (shift<=max_shift) and\
                        (overlap_ratio >= min_overlap_ratio) and\
                        (size_similarity >= min_size_similarity):
                            cluster[j]=cluster[i]
                            
    cluster = np.array(cluster)
    dc = Counter(cluster)
    valid_cluster = []
    for idx in dc:
        if dc[idx]>=1:
            valid_cluster.append(idx)
            
    final_sig_list = []
    for cluster_idx in valid_cluster:
        sig_idxs = np.where(cluster==cluster_idx)
#         print(sig_idxs)
        best_sig = sig_list[sig_idxs[0][0]]
        for idx in sig_idxs[0]:
            new_sig = sig_list[idx]
            best_sig_len = best_sig[3]
            new_sig_len = new_sig[3]
            # update best sig
            if new_sig_len > best_sig_len:
                best_sig = new_sig
        final_sig_list.append(best_sig)
        
    return final_sig_list
    
def cluster_ins(sig_list,max_shift = 100, min_size_similarity = 0.5):
    cluster = [-1]*len(sig_list)
    while (np.array(cluster)==-1).sum()!=0:
        for i in range(len(sig_list)):
            if cluster[i]==-1:
                sig1 = sig_list[i]
                cluster[i] = i
                for j in range(len(sig_list)):
                    if cluster[j]==-1:
                        sig2 = sig_list[j]
                        size_similarity = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
                        shift = abs(sig1[2]-sig2[2])
                        
                        if (shift<=max_shift) and\
                        (size_similarity >= min_size_similarity):
                            cluster[j]=cluster[i]                           
    cluster = np.array(cluster)
    dc = Counter(cluster)
    valid_cluster = []
    for idx in dc:
        if dc[idx]>=1:
            valid_cluster.append(idx)
            
    final_sig_list = []
    for cluster_idx in valid_cluster:
        sig_idxs = np.where(cluster==cluster_idx)
#         print(sig_idxs)
        best_sig = sig_list[sig_idxs[0][0]]
        for idx in sig_idxs[0]:
            new_sig = sig_list[idx]
            best_sig_len = best_sig[3]
            new_sig_len = new_sig[3]
            # update best sig
            if new_sig_len > best_sig_len:
                best_sig = new_sig
        final_sig_list.append(best_sig)
        
    return final_sig_list

def get_read_start_end(cigar):
    rl =0
    for tp in cigar:
        if tp[0] in {0,1,4,5}:
            rl+=tp[1]

    start , end = 0, rl

    if cigar[0][0] in {4,5}:
        start = cigar[0][1]


    if cigar[-1][0] in {4,5}:
        end = rl - cigar[-1][1]

    return start,end
    
def extract_sig_from_split(read1,read2,min_mapq,max_svlen):  
    def get_readlen(cigar):
        rl =0
        for tp in cigar:
            if tp[0] in {0,1,4,5}:
                rl+=tp[1]
        return rl
    
    assert read1.pos<=read2.pos
    assert read1.qname == read2.qname
    assert read1.reference_name == read2.reference_name
    reverse1,start1,end1,cigar1,seq1,mapq1 = \
    read1.is_reverse,read1.pos,read1.reference_end,read1.cigar,read1.seq,read1.mapq
    reverse2,start2,end2 ,cigar2,seq2,mapq2= \
    read2.is_reverse,read2.pos,read2.reference_end,read2.cigar,read2.seq,read2.mapq
    del_list = []
    ins_list = []
#     print("judge0")
    if (reverse1==reverse2) and ( (mapq1>=min_mapq) and (mapq2>=min_mapq)  )\
    and (cigar1[-1][0] in {4,5}) and (cigar2[0][0] in {4,5}):
#         print("judge")
        
        rl1 = get_readlen(cigar1)
        rl2 = get_readlen(cigar2)
        assert rl1==rl2
        
        Ref1e = end1
        Ref2s = start2
 
        Read1e = rl1-cigar1[-1][1]
        Read2s = cigar2[0][1]

            
        Diffdis = (Ref2s-Ref1e)-(Read2s-Read1e)
        Diffolp = Ref1e - Ref2s
        
        if reverse1:
            direction = '-'
        else:
            direction = '+'
            


        if abs(Diffdis)<=max_svlen:
            if (Diffolp<30)  and (Diffdis >= 30 ):
                sigdel = [read1.reference_name,'DEL', Ref1e, Diffdis,read1.qname, Read1e,Read2s,direction,'split-alignment',"%d-%d"%(read1.mapq,read2.mapq)]
                del_list.append(sigdel)
            elif (Diffolp<3000)  and (Diffdis >= 30 ):
                sigdel = [read1.reference_name,'DEL', Ref1e-Diffdis, Diffdis,read1.qname, Read1e-Diffdis,Read2s-Diffdis,direction,'split-alignment',"%d-%d"%(read1.mapq,read2.mapq)]
                del_list.append(sigdel)
            elif (Diffolp<3000) and (Diffdis<= -30) :
                svlen = abs(Read2s-Read1e+Diffolp)
                if abs(Diffolp) > 400:
                    pos_ref = int((Ref1e+Ref2s)/2)
                else:
                    pos_ref = Ref2s
                sigins = [read1.reference_name,'INS', pos_ref, svlen , read1.qname,Read1e-Diffolp, Read2s,direction,'split-alignment',"%d-%d"%(read1.mapq,read2.mapq)]

                ins_list.append(sigins)
    return del_list,ins_list

                    

def write_sig_split(sig_list,output_path):
    sig_array = np.array(sig_list)
    with open(output_path,'w') as f:
        for sig in sig_array:
            line = '\t'.join(sig)+'\n'
            
            f.write(line)
    return             

    

def extract_signature_from_cigar(bam_path,chr_name,output_dir,hp,min_cigar_mapq):
    samfile = pysam.AlignmentFile(bam_path)
    del_sig_cigar = []
    ins_sig_cigar = []
    name_list = []
    for read in samfile.fetch(chr_name):
        if hp in read.qname:
            name_list.append(read.qname)
            if read.mapq>=min_cigar_mapq:
                del_sig,ins_sig,offset_ref,offset_contig = extract_sig_from_cigar(read,min_svlen=30)
                assert offset_ref==read.reference_end
                if read.seq:
                    assert len(read.seq)==offset_contig
                del_sig_cigar.extend(del_sig)
                ins_sig_cigar.extend(ins_sig)

    del_sig_cigar_sorted = sort_sig(del_sig_cigar)
    ins_sig_cigar_sorted = sort_sig(ins_sig_cigar)
    write_sig_cigar(del_sig_cigar_sorted,"%s/%s_DEL_contig_cigar_%s.txt"%(output_dir,chr_name,hp) )
    write_sig_cigar(ins_sig_cigar_sorted,"%s/%s_INS_contig_cigar_%s.txt"%(output_dir,chr_name,hp) )

    del_sig_cigar_sorted_clustered = cluster_del(del_sig_cigar_sorted,
                      max_shift = 100, 
                      min_overlap_ratio = 0.5, 
                      min_size_similarity = 0.5
                )
    ins_sig_cigar_sorted_clustered = cluster_ins(ins_sig_cigar_sorted,
                      max_shift = 100, 
                      min_size_similarity = 0.5
                )

    print(len(del_sig_cigar_sorted),len(ins_sig_cigar_sorted))
    print(len(del_sig_cigar_sorted_clustered),len(ins_sig_cigar_sorted_clustered))
    return del_sig_cigar_sorted_clustered,ins_sig_cigar_sorted_clustered

def extract_sig_from_split_reads(bam_path,chr_name,output_dir,hp,min_split_mapq):
    samfile = pysam.AlignmentFile(bam_path)
    name_list = []
    for read in samfile.fetch(chr_name):
        if hp in read.qname:
            if read.mapq >= min_split_mapq:
                name_list.append(read.qname)
    dc = Counter(name_list)
    split_name_dc ={}
    for name in dc:
        if dc[name]>1:
            split_name_dc[name]=dc[name]


    samfile = pysam.AlignmentFile(bam_path)
    name_list = []
    read_list =[]
    for read in samfile.fetch(chr_name):
        if read.qname in split_name_dc:
            if read.mapq >= min_split_mapq:
                read_list.append(read)
                name_list.append(read.qname)
    name_array = np.array(name_list)     

    del_sig_split = []
    ins_sig_split = []
    for name in split_name_dc:
        idxs = np.where(name_array==name)[0]
    #     print(idxs)

        for i in range(len(idxs)-1):
    #         print(name)
            read1 = read_list[idxs[i]]
            read2 = read_list[idxs[i+1]]
            del_list,ins_list = extract_sig_from_split(read1,read2,min_mapq=min_split_mapq,max_svlen = 50000)
            del_sig_split.extend(del_list)
            ins_sig_split.extend(ins_list)

    del_sig_split_sorted = sort_sig(del_sig_split)
    ins_sig_split_sorted = sort_sig(ins_sig_split)

    write_sig_split(del_sig_split_sorted,output_dir+"/%s_DEL_contig_split_%s.txt"%(chr_name,hp) )
    write_sig_split(ins_sig_split_sorted,output_dir+"/%s_INS_contig_split_%s.txt"%(chr_name,hp) )

    del_sig_split_sorted_clustered = cluster_del(del_sig_split_sorted,
                      max_shift = 100, 
                      min_overlap_ratio = 0.5, 
                      min_size_similarity = 0.5
                )
    ins_sig_split_sorted_clustered = cluster_ins(ins_sig_split_sorted,
                      max_shift = 100, 
                      min_size_similarity = 0.5
                )
    print(len(del_sig_split_sorted),len(ins_sig_split_sorted))
    print(len(del_sig_split_sorted_clustered),len(ins_sig_split_sorted_clustered))
    return del_sig_split_sorted_clustered,ins_sig_split_sorted_clustered

def merge_sig_ins(ins_sig_cigar,ins_sig_split):
    sig_list = ins_sig_cigar+ins_sig_split
    sig_list_sorted = sort_sig(sig_list)
    sig_list_clustered = cluster_ins(sig_list_sorted)
    print(len(sig_list_sorted),len(sig_list_clustered))
    return sig_list_clustered

def merge_sig_del(del_sig_cigar,del_sig_split):
    sig_list = del_sig_cigar+del_sig_split
    sig_list_sorted = sort_sig(sig_list)
    sig_list_clustered = cluster_del(sig_list_sorted)
    print(len(sig_list_sorted),len(sig_list_clustered))
    return sig_list_clustered

def merge_all(del_sig_cigar,ins_sig_cigar,del_sig_split,ins_sig_split):
    ins_final = merge_sig_ins(ins_sig_cigar,ins_sig_split)
    del_final = merge_sig_del(del_sig_cigar,del_sig_split)
    final_sig_sorted = sort_sig(ins_final+del_final)
    print("final INS signal %d"%(len(ins_final)))
    print("final DEL signal %d"%(len(del_final)))
    print("final signal %d"%(len(final_sig_sorted)))
    return final_sig_sorted

def extract_signature_one_hap(bam_path,chr_name,output_dir,hp):
    ### extract from cigar
    del_sig_cigar_sorted_clustered,ins_sig_cigar_sorted_clustered = \
    extract_signature_from_cigar(bam_path,chr_name,output_dir,hp,min_cigar_mapq)
    ### extract from split reads
    del_sig_split_sorted_clustered,ins_sig_split_sorted_clustered = \
    extract_sig_from_split_reads(bam_path,chr_name,output_dir,hp,min_split_mapq)
    ### merge both information
    final_sig_sorted = merge_all(del_sig_cigar_sorted_clustered,
                                 ins_sig_cigar_sorted_clustered,
                                 del_sig_split_sorted_clustered,
                                 ins_sig_split_sorted_clustered)
    return final_sig_sorted

def pair_ins(sig1,sig2, max_shift, min_size_similarity):
    assert sig1[:2]==sig2[:2]
    shift = abs(sig1[2]-sig2[2])
    size_sim = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
    if (shift<= max_shift) and (size_sim>=min_size_similarity):

        return 1
    else:
        return 0
            
def pair_del(sig1,sig2, max_shift, min_overlap_ratio,min_size_similarity):
    def calculate_overlap_ratio(start1,end1,start2,end2):
        len1 = end1-start1
        len2 = end2-start2
        minlen = min(len1,len2)
        return (min(end1,end2)-max(start1,start2))/minlen
    assert sig1[:2]==sig2[:2]
    start1 = sig1[2]
    start2 = sig2[2]
    end1 = start1+sig1[3]
    end2 = start2+sig2[3]

    overlap_ratio = calculate_overlap_ratio(start1,end1,start2,end2)
    size_similarity = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
    shift = abs(sig1[2]-sig2[2])

    if (shift<=max_shift) and\
    (overlap_ratio >= min_overlap_ratio) and\
    (size_similarity >= min_size_similarity):
        return 1
    else:
        return 0
    
def pair_sig(sig_hp1,sig_hp2,max_compare_dist,max_shift, min_overlap_ratio,min_size_similarity):
    pair_status_hp1 = [-1]*len(sig_hp1)
    pair_status_hp2 = [-1]*len(sig_hp2)
    
    for i in range(len(sig_hp1)):
        sig1 = sig_hp1[i]
        for j in range(len(sig_hp2)):
            sig2 = sig_hp2[j]
            dist = sig2[2]-sig1[2]
            if dist>max_compare_dist:
#                 print(i,j,dist)
                break
            elif (sig1[:2]==sig2[:2]) and (pair_status_hp2[j]==-1):
                if sig1[1]=='DEL':
                    result = pair_del(sig1,sig2,200,0.5,0.5)
                else:
                    result = pair_ins(sig1,sig2,200,0.5)
#                 print(i,j,result)
                if result==1:
                    pair_status_hp1[i] = j
                    pair_status_hp2[j] = i
                    break
    paired_sig =[]        
    for i in range(len(sig_hp1)):
        sig1 = sig_hp1[i]
        sig1_contig_info = "%s:%d-%d"%(sig1[4],sig1[5],sig1[6])
        if pair_status_hp1[i]==-1:
            paired_sig.append(sig1+['0/1',sig1_contig_info,sig1[7], sig1[8], str(sig1[9]) ])
        else:
            j = pair_status_hp1[i]
            sig2 = sig_hp2[j]
            sig2_contig_info = "%s:%d-%d"%(sig2[4],sig2[5],sig2[6])
            contig_info = sig1_contig_info+','+sig2_contig_info
            sig_source = sig1[8]+','+sig2[8]
            sig_mapq = str(sig1[9])+','+str(sig2[9])
            if sig1[3]>sig2[3]:
                paired_sig.append(sig1+['1/1',contig_info,sig1[7]+','+sig2[7], sig_source, sig_mapq])
            else:
                paired_sig.append(sig2+['1/1',contig_info,sig1[7]+','+sig2[7], sig_source, sig_mapq ])
                
    for i in range(len(sig_hp2)):
        sig2 = sig_hp2[i]
        if pair_status_hp2[i]==-1:
            sig2_contig_info = "%s:%d-%d"%(sig2[4],sig2[5],sig2[6])
            paired_sig.append(sig2+['0/1', sig2_contig_info,sig2[7], sig2[8], str(sig2[9])])
            
    paired_sig_sorted = sort_sig(paired_sig)
    homo_cnt =0
    for sig in paired_sig_sorted:
        if sig[10]=='1/1':
            homo_cnt +=1
    heter_cnt = len(paired_sig_sorted)-homo_cnt
    print("homo variants %d"%homo_cnt)
    print("heter variants %d"%heter_cnt)
    print("total variant %d"%len(paired_sig_sorted))
    return paired_sig_sorted  

def load_contigs(fasta_path):
    logger.info("loading "+fasta_path)
    with open(fasta_path,'r') as f:
        s = f.readlines()
    dc = {}

    for line in tqdm(s,desc='load contig'):
        if '>' in line:
            cur_name = line[1:-1]
        else:
            if cur_name in dc:
                dc[cur_name].append(line[:-1])
            else:
                dc[cur_name] = [line[:-1]]

    for name in tqdm(dc):
        dc[name] = ''.join(dc[name])

    logger.info("finish loading")
    return dc 

def reverse_compelement(seq):
    seq1 = seq.upper()[::-1]
    map_base = {'N':'N','A':'T','T':'A','G':'C','C':'G'}
    newseq = ''
    for i in range(len(seq1)):
        newseq+=map_base[seq1[i]]
    return newseq

def load_seq(fasta_path):
    logger.info("loading "+fasta_path)
    with open(fasta_path,'r') as f:
        s = f.read().split('\n')[1:-1]
    logger.info("finish loading")
    return ''.join(s)


def add_seq_to_sig(paired_sig_sorted,ref_path,contig_path):
    # ref_seq = load_seq(ref_path)
    # dc_contig = load_contigs(contig_path)
    sig_with_seq = []
    for i in tqdm(range(len(paired_sig_sorted)),desc = "add seq"):
        sig = paired_sig_sorted[i].copy()
        if sig[4] in dc_contig:
            if sig[1]=='INS':
                seq_contig = dc_contig[sig[4]]
                # if len(sig)==9:
                #     inserted_seq = seq_contig[sig[5]:sig[6]]
                # else:
                #     inserted_seq = seq_contig[sig[5]:sig[5]+sig[3]]

                # if sig[7]=='-':
                #     if len(sig)==9:
                #         inserted_seq = reverse_compelement(seq_contig[-sig[6]:-sig[5]])
                #     else:
                #         inserted_seq = reverse_compelement(seq_contig[-sig[5]-sig[3]:-sig[5]])
                # else:
                #     if len(sig)==9:
                #         inserted_seq = seq_contig[sig[5]:sig[6]]
                #     else:
                #         inserted_seq = seq_contig[sig[5]:sig[5]+sig[3]]
                if sig[7]=='-':
                    inserted_seq = reverse_compelement(seq_contig[-sig[6]:-sig[5]])
                else:
                    inserted_seq = seq_contig[sig[5]:sig[6]]
                sig.append(inserted_seq)
            else:
                deleted_seq = ref_seq[sig[2]:sig[2]+sig[3]]
                sig.append(deleted_seq)
            sig_with_seq.append(sig)
    return sig_with_seq
            
            
def write_vcf(paired_sig_sorted,vcf_path,ref_path,contig_path):
    sig_with_seq = add_seq_to_sig(paired_sig_sorted,ref_path,
                              contig_path)
    with open(header_path,"r") as f:
        header = f.readlines()
    fv = open(vcf_path,'w')
    fv.writelines(header)
    # seq_chr = load_seq(ref_path)
    seq_chr = ref_seq
    ins_cnt = 0
    del_cnt = 0
    for sig in sig_with_seq:
        chr_name = sig[0]
        svtype = sig[1]
        pos = sig[2]-1
        pos_one_base = sig[2]

        if svtype == 'DEL':
            alt_allele = seq_chr[pos]
            ref_allele = alt_allele+sig[-1]
            del_cnt+=1
            index_cnt = del_cnt

        else:
            ref_allele = seq_chr[pos]
            alt_allele = ref_allele+sig[-1]
            ins_cnt+=1
            index_cnt = ins_cnt

        svlen = len(alt_allele)-len(ref_allele)
        svinfo = "SVLEN=%d;SVTYPE=%s;TIG_REGION=%s;QUERY_STRAND=%s;SIG_SOURCE=%s;TIG_MAPQ=%s"%(svlen,svtype,sig[11],sig[12],sig[13],sig[14])
        gt = sig[10]

        line = chr_name+'\t'+str(pos_one_base)+'\tvolcano.%s.%s.%d\t%s\t%s\t%d\tPASS\t%s\tGT\t%s\n'%(chr_name,svtype, index_cnt,ref_allele.upper(),alt_allele.upper(),20,svinfo, gt)
        fv.write(line)
    fv.close()
    return 
     

# bam_path = '/data/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/CanLuo/\
# long_reads_project/subsampled_result/12x_result/final_contigs/NA24385_k12_onehot_cov12_asm_sorted.bam'

# contig_path ='/data/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/CanLuo/long_reads_project/subsampled_result/12x_result/final_contigs/final_contig.p_ctg.fa'

# bam_path = "./Na24385_contigs.bam"
# contig_path = "./final_contig.p_ctg.fa"

signature_dir = output_dir+'/signature/'
os.system("mkdir -p "+signature_dir)


if chr_number is not None:
    start_i = chr_number
    end_i = chr_number+1 
else:
    start_i = 1 
    end_i = 23

logger.info("load asm contigs...")
dc_contig = load_contigs(contig_path)
print(dc_contig.keys())
logger.info("load reference...")
dc_ref = load_contigs(ref_path)

for i in range(start_i, end_i) :
    chr_name = 'chr'+str(i)
    chr_number = i
    print(chr_name)
    ref_seq = dc_ref[chr_name]
    fw = open(signature_dir+'/%s_cigar.txt'%chr_name,'w')
    samfile = pysam.AlignmentFile(bam_path)
    samiter = samfile.fetch(chr_name)

    for read in samiter:
        start,end = get_read_start_end(read.cigar)
        if read.is_reverse:
            direction = '-'
        else:
            direction = '+'
        print(read.qname+'\t'+str(read.mapq)+'\t'+str(read.pos)+'\t'+str(read.reference_end)+'\t'+"%d\t%d\t%s\t"%(start,end,direction)+str(read.cigar),file = fw)
    fw.close()

    final_sig_hp1 = extract_signature_one_hap(bam_path,
                                 chr_name,output_dir+'/signature/',
                                hp ='hp1')

    final_sig_hp2 = extract_signature_one_hap(bam_path,
                                 chr_name,output_dir+'/signature/',
                                hp ='hp2')

    paired_sig_sorted = pair_sig(final_sig_hp1,final_sig_hp2,1000,200,0.5,0.5)

    write_vcf(paired_sig_sorted,
              vcf_path=output_dir+"/volcano_variant_chr%d.vcf"%chr_number,
              ref_path=ref_path,
              contig_path = contig_path)