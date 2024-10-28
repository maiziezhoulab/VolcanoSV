import pysam
import numpy as np 
from tqdm import tqdm 
from joblib import Parallel, delayed
import os
from collections import defaultdict

def load_indel_kmer_file(file_path):
    """
    Load the indel kmer file and return a list of dictionaries containing the information.

    Args:
        file_path (str): The path to the indel kmer file.

    Returns:
        list: A list of dictionaries containing the information from the indel kmer file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    results = []
    cnt = 0 
    for line in lines:
        line = line.strip()
        if cnt==0:
            cnt+=1
            continue
        
        if not line:
            continue
        fields = line.split('\t')
        result = {
            'svid': fields[0],
            'qname': fields[1],
            # 'start': int(fields[2]),
            # 'end': int(fields[3]),
            'contig_seq': fields[4],
            'alt': fields[5]
        }
        results.append(result)

    return results



def kmer_sequencer(seq, kmer_size):
    kmers  = []
    for i in range(len(seq) - kmer_size + 1):
        kmer = seq[i:i+kmer_size]
        kmers.append(kmer)
    return kmers


# def count_kmers_in_region(bam_file, region, kmer_size):
#     kmers = {}
#     bam = pysam.AlignmentFile(bam_file, "rb")
#     chrom, coords = region.split(":")
#     start, end = coords.split("-")
#     start, end = int(start), int(end)
#     x = 0
#     for read in bam.fetch(chrom, start, end):
#         seq = read.seq.upper()
#         for i in range(len(seq) - kmer_size + 1):
#             kmer = seq[i:i+kmer_size]
#             if kmer in kmers:
#                 kmers[kmer] += 1
#             else:
#                 kmers[kmer] = 1
#         x+=1
#     # print(f"num reads in {region}: {x}")
#     # print(len(kmers))

#     return kmers

def get_seq(bamfile, chrom, start_pos, end_pos):
    samfile = pysam.AlignmentFile(bamfile)
    dc = defaultdict(list)
    for pileupcolumn in samfile.pileup(chrom, start_pos, end_pos, min_base_quality = 0, truncate = True): 
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                # Deletion or reference skip

                continue
            elif pileupread.indel > 0:
                # Insertion
                seq = pileupread.alignment.query_sequence[pileupread.query_position:\
                                                          pileupread.query_position+pileupread.indel+1]
                dc[pileupread.alignment.query_name].append(seq)


            else:
                # No indel
                seq = pileupread.alignment.query_sequence[pileupread.query_position]
                dc[pileupread.alignment.query_name].append(seq)

    for qn,bs in dc.items():
        dc[qn] = ''.join(bs).upper()
    samfile.close()
    return dc


def count_kmers_in_region(bam_file, region, kmer_size):
    kmers = {}
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    start, end = int(start), int(end)
    x = 0
    dc = get_seq(bam_file, chrom, start, end)
    for qname, seq in dc.items():
        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i:i+kmer_size]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
        x+=1
    # print(f"num reads in {region}: {x}")
    # print(len(kmers))

    return kmers


def check_one_var(var,bamfile,kmer_size):
    svid = var['svid']
    chrom = svid.split('-')[0]
    pos = int(svid.split('-')[-1])
    region =f"{chrom}:{pos-20}-{pos+70}"
    # print("region:",region)
    kmer_counter_reads = count_kmers_in_region(bamfile, region, kmer_size)
    kmers_contig = kmer_sequencer(var['contig_seq'],kmer_size)
    # print("kmers from contig:",kmers_contig)
    # kmer_count  = [ (kmer,kmer_counter_reads[kmer]) if kmer in kmer_counter_reads else (kmer,0) for kmer in kmers_contig ]
    kmer_count  = [ kmer_counter_reads[kmer] if kmer in kmer_counter_reads else 0 for kmer in kmers_contig ]
    # print("kmer count in reads:",kmer_count)
    return kmer_count[7:-7]


def load_index(vcf_file):
    index_list = []
    with open(vcf_file,'r') as f:
        for line in f:
            if line[0]!='#':
                index_list.append(line.split()[2])

    return set(index_list) 

def write_kmer_cnt(kmer_cnt_list, svid_list, outfile):

    with open(outfile,'w') as f:
        for i in range(len(svid_list)):
            kmer_cnt = [str(k) for k in kmer_cnt_list[i]]
            kmer_str = '\t'.join(kmer_cnt)
            line = f"{svid_list[i]}\t{kmer_str}\n"
            f.write(line)

    return

def load_kmer_cnt(outfile):
    data_list = []
    
    with open(outfile,'r') as f:
        for line in f:
            data =[ int(x) for x in line[:-1].split()[1:]]
            data_list.append(data)

    return data_list


def load_summary(eval_dir):
    with open(eval_dir+'/summary.txt','r') as f:
        s = f.readlines()
    dc = {}
    for line in s:
        data = line.strip().split(', ')
        for x in data:
            [key,val] = x.split(': ')
            val = eval(val)
            dc[key]= val 
    return dc

def filter_indel( prefix,vcffile, indel_kmer_file, reads_bamfile, output_folder, 
                  eval_dir = None,
                  kmer_size = 15,
                   n_thread = 30, ratio = 0.3, min_support = 5, restart = False ):

    os.system('mkdir -p '+output_folder)
    # chrom = 'chr'+chr_number


    var_list = load_indel_kmer_file(indel_kmer_file)
    final_list = []
    data_list = []

    svid_list = np.array([var['svid'] for var in var_list])
    outfile = output_folder+'/kmer_count_bysvid.txt'

    if not restart:
        sequences = Parallel(n_jobs=n_thread)(delayed(check_one_var)(var_list[i],reads_bamfile,kmer_size) for i in tqdm(range(len(var_list))))
        write_kmer_cnt(sequences , svid_list, outfile)
    else:
        sequences = load_kmer_cnt(outfile)

    # add debugger
    for i in range(len(var_list)):
        if len(len(sequences[i])) == 0:
            print("encounter empty kmer support array!")
            print("Var information: ", f"{i}-th var ",var_list[i])
            exit()


    if eval_dir:

        for cnt_list in tqdm(sequences):
            zero_count = ((np.array(cnt_list)<= min_support )).sum()
            data_list.append(cnt_list)
            final_list.append(zero_count/len(cnt_list))

        final_list = np.array(final_list)
        # ratio = 0.3
        print(f"under {min_support} ratio: {ratio}")
        rm_svids = set(svid_list[final_list> ratio])


        #### check how many TPs and FPs are removed

        tp_index = load_index(eval_dir+'/tp-call.vcf')
        fp_index = load_index(eval_dir+'/fp.vcf')
        dc_eval = load_summary(eval_dir)


        tp_cnt = 0
        fp_cnt = 0
        for svid in rm_svids:
            if svid in tp_index:
                tp_cnt+=1
            else:
                fp_cnt +=1
        print("remove tp: ",tp_cnt, "remove fp: ",fp_cnt)
        tp = dc_eval['TP'] - tp_cnt
        fp = dc_eval['FP'] - fp_cnt 
        fn = dc_eval['FN'] + tp_cnt
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * recall * precision / ( recall + precision )
        # print("=====old result====")
        # for key in dc_eval:
        #     print(key,': ', dc_eval[key])
        # print("=====new result====")
        print("TP:",tp)
        print("FP:",fp)
        print("FN:",fn)
        print("Recall:", recall)
        print("Precision:",precision)
        print("F1:", f1)

        # print(final_list)
        # print((final_list>0).sum())
        # print((final_list==0).sum())

    if vcffile:
        cnt = 0
        with open(vcffile,'r') as f:
            with open(output_folder+f"/{prefix}_volcanosv_small_indel.vcf",'w') as fw:
                for line in f:
                    
                    if line[0]=='#':
                        ### add PS tag ================== to do !!
                        if "ID=CONTEXT" in line :
                            line += '''##INFO=<ID=PS,Number=.,Type=Integer,Description="phase block name">\n'''
                        elif "#CHROM" in line:
                            line = line.replace("syndip",prefix)
                        fw.write(line)
                    else:
                        # svid = line.split()[2]
                        # if svid not in rm_svids:

                        cnt_list = sequences[cnt]
                        cnt+=1
                        zero_rate = ((np.array(cnt_list)<= min_support )).mean()
                        if zero_rate <= ratio:
                            data = line.split()
                            gt_og = data[-1].split(':')[0]
                            gt_ele = set(gt_og.split('|'))
                            if gt_ele - {0,1,'.'}:
                                # multi-ale
                                if "hp1" in data[7]:
                                    gt = "1|0"
                                else:
                                    gt = "0|1"
                            else:
                                gt = gt_og 
                            
                            data[-1] = gt 
                            data[-2] = 'GT'
                            ps = line.split("TIG_REGION=PS")[1].split('_')[0]
                            data[7]= data[7]+";PS=" + ps 
                            line = "\t".join(data)+'\n'
                            fw.write(line)


    return 




import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count k-mer occurrences in reads overlapping indels",
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--indel_kmer_file", '-i',type=str, help="Path to indel_kmer.txt file")
    parser.add_argument("--vcffile", '-v',type=str, help="Path to original vcf file")
    parser.add_argument("--eval_dir", '-e',type=str, help="Path to evaluation directory")
    parser.add_argument("--reads_bamfile", '-b',type=str, help="Path to BAM file containing reads")
    parser.add_argument("--output_folder",'-o', type=str, help="Output folder for results")
    # parser.add_argument('--chr_number','-chr', type = str)
    parser.add_argument("--kmer_size",'-k', type=int, help="K-mer size to count", default = 15)
    parser.add_argument("--n_thread",'-t', type = int, default = 50)
    parser.add_argument("--ratio",'-r', type = float, default = 0.3, help = "kmer support ratio")
    parser.add_argument('--min_support','-ms', type = int, default = 5, help = "min kmer support per base")
    parser.add_argument('--prefix','-px', help = "file prefix in the output folder", default = "Sample")
    parser.add_argument('--restart','-rs', action='store_true', help = "restart mode; assume there is kmer support file already.")
    args = parser.parse_args()

    # Extract variables from args
    indel_kmer_file = args.indel_kmer_file
    vcffile = args.vcffile
    eval_dir = args.eval_dir
    reads_bamfile = args.reads_bamfile
    output_folder = args.output_folder
    kmer_size = args.kmer_size
    # chr_number = args.chr_number 
    n_thread = args.n_thread
    ratio = args.ratio
    min_support = args.min_support
    restart = args.restart
    prefix=args.prefix

    filter_indel( prefix,vcffile, indel_kmer_file, reads_bamfile, output_folder, 
                  eval_dir ,
                  kmer_size ,
                   n_thread , ratio , min_support , restart  )





