from argparse import ArgumentParser
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_number','-chr',type = int)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')

parser.add_argument('--max_shift',type = int, default = 100)
parser.add_argument('--max_shift_ratio',type = float, default = 0.1)
parser.add_argument('--min_reads_support',type = int, default = 1)
parser.add_argument('--min_siglen',type =int, default = 30)

args = parser.parse_args()
chr_number = args.chr_number
input_path = args.input_path
output_dir = args.output_dir


max_shift = args.max_shift
max_shift_ratio = args.max_shift_ratio
min_reads_support = args.min_reads_support
min_siglen = args.min_siglen



chr_name = "chr%d"%chr_number

import matplotlib.pyplot as plt
import numpy as np
import pickle
from tqdm import tqdm
import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import os
import logging
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
        if tp[0] in [0,7,8]:
            offset_ref+=tp[1]
            offset_contig +=tp[1]
        elif tp[0]==4:
            offset_contig +=tp[1]
        elif tp[0]==2:
            if tp[1]>=min_svlen:
                del_sig.append([ref_name,'DEL',offset_ref,tp[1],qname,offset_contig + hard_clip_head,direction,'cigar'])  
            offset_ref+=tp[1]
        elif tp[0] == 1:   
            if tp[1]>=min_svlen:
                ins_sig.append([ref_name,'INS',offset_ref, tp[1],qname,offset_contig + hard_clip_head,direction,'cigar'])
            offset_contig +=tp[1]
        elif tp[0] == 3:  # Skipped region
            offset_ref += tp[1]  # Advance reference without consuming contig bases

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
            line = '\t'.join(sig)+'\n'
            f.write(line)
    return

def extract_signature_from_cigar(bam_path,chr_name,output_dir):
    samfile = pysam.AlignmentFile(bam_path)
    del_sig_cigar = []
    ins_sig_cigar = []
    name_list = []
    cnt = 0
    for read in samfile.fetch(chr_name):
        if cnt%10000==0:
            logger.info("Processed %d reads"%cnt)
        cnt+=1
#         if hp in read.qname:
        name_list.append(read.qname)
    
        if read.mapq>=50:
            del_sig,ins_sig,offset_ref,offset_contig = extract_sig_from_cigar(read,min_svlen=30)
            assert offset_ref==read.reference_end
            if read.seq:
                assert len(read.seq)==offset_contig
            del_sig_cigar.extend(del_sig)
            ins_sig_cigar.extend(ins_sig)

    del_sig_cigar_sorted = sort_sig(del_sig_cigar)
    ins_sig_cigar_sorted = sort_sig(ins_sig_cigar)
    write_sig_cigar(del_sig_cigar_sorted,"%s/%s_DEL_reads_cigar.txt"%(output_dir,chr_name) )
    write_sig_cigar(ins_sig_cigar_sorted,"%s/%s_INS_reads_ciga.txt"%(output_dir,chr_name) )

#     del_sig_cigar_sorted_clustered = cluster_del(del_sig_cigar_sorted,
#                       max_shift = 100, 
#                       min_overlap_ratio = 0.5, 
#                       min_size_similarity = 0.5
#                 )
#     ins_sig_cigar_sorted_clustered = cluster_ins(ins_sig_cigar_sorted,
#                       max_shift = 100, 
#                       min_size_similarity = 0.5
#                 )

    print(len(del_sig_cigar_sorted),len(ins_sig_cigar_sorted))
#     print(len(del_sig_cigar_sorted_clustered),len(ins_sig_cigar_sorted_clustered))
    return del_sig_cigar_sorted,ins_sig_cigar_sorted

def extract_sig_from_split(read1,read2,min_mapq,max_svlen):  
    def get_readlen(cigar):
        rl =0
        for tp in cigar:
            if tp[0] in {0,7,8,1,4,5}:
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
                sigdel = [read1.reference_name,'DEL', Ref1e, Diffdis,read1.qname, Read1e,Read2s,direction,'split-alignemnt']
                del_list.append(sigdel)
            elif (Diffolp<30) and (Diffdis<= -30) :
                sigins = [read1.reference_name,'INS', int((Ref1e+Ref2s)/2),abs(Diffdis), read1.qname,Read1e, Read2s,direction,'split_alignment']
                ins_list.append(sigins)
    return del_list,ins_list

def extract_sig_from_split_reads(bam_path,chr_name,output_dir):
    samfile = pysam.AlignmentFile(bam_path)
    name_list = []
    for read in samfile.fetch(chr_name):
        name_list.append(read.qname)
    dc = Counter(name_list)
    split_name_dc ={}
    for name in dc:
        if dc[name]>1:
            split_name_dc[name]=dc[name]


    samfile = pysam.AlignmentFile(bam_path)
    name_list = []
    read_list =[]
    cnt = 0
    for read in samfile.fetch(chr_name):
        if cnt%10000==0:
            logger.info("Processed %d reads"%cnt)
        cnt+=1
        if read.qname in split_name_dc:

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
            del_list,ins_list = extract_sig_from_split(read1,read2,min_mapq=0,max_svlen = 50000)
            del_sig_split.extend(del_list)
            ins_sig_split.extend(ins_list)

    del_sig_split_sorted = sort_sig(del_sig_split)
    ins_sig_split_sorted = sort_sig(ins_sig_split)

    write_sig_split(del_sig_split_sorted,output_dir+"/%s_DEL_reads_split.txt"%(chr_name) )
    write_sig_split(ins_sig_split_sorted,output_dir+"/%s_INS_reads_split.txt"%(chr_name) )


    print(len(del_sig_split_sorted),len(ins_sig_split_sorted))

    return del_sig_split_sorted,ins_sig_split_sorted


def write_sig_split(sig_list,output_path):
    # sig_array = np.array(sig_list)
    with open(output_path,'w') as f:
        for sig in sig_list:
            sig1 =[str(p) for p in sig]
            line = '\t'.join(sig1)+'\n'
            
            f.write(line)
    return   

def merge_all(sig_list,output_dir,chr_name):
    output_path = output_dir + '/'+chr_name+'_reads_sig.txt'
    sig_list_sorted = sort_sig(sig_list)
    write_sig_split(sig_list_sorted,output_path)
    return sig_list_sorted


output_dir = output_dir+'/reads_signature'
os.system("mkdir -p "+output_dir)

# bam_path = '/data/h_constantinidis_lab/ukbb/MaizieZhouLab_backup/CanLuo/long_reads_project/subsampled_result/3x_result/%s/NA24385_aligned_by_ngmlr_%s.bam'%(chr_name,chr_name)

bam_path = input_path

del_sig_cigar_sorted,ins_sig_cigar_sorted = \
extract_signature_from_cigar(bam_path,chr_name,output_dir)

del_sig_split_sorted,ins_sig_split_sorted = \
extract_sig_from_split_reads(bam_path,chr_name,output_dir)

sig_list = del_sig_cigar_sorted+\
ins_sig_cigar_sorted+\
del_sig_split_sorted+\
ins_sig_split_sorted

merged_sig = merge_all(sig_list,output_dir, chr_name)

