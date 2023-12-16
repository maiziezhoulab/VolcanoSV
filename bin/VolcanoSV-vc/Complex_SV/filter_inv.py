import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcffile','-v')
parser.add_argument('--output_dir','-o')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--min_support','-ms',type = int, default = 1, help = "minimum reads support per breakend")
parser.add_argument('--min_mapq','-mq',type = int, default = 0, help = "minimum support read mapping quality")
parser.add_argument('--flanking','-fl',type = int, default = 1000, help = "flanking region for support reads checking")
parser.add_argument('--n_thread','-t', type = int, default = 22 )
parser.add_argument('--max_dist','-d', type = int, default = 500, help = "maximum distance to merge 2 INVs")
args = parser.parse_args()
vcffile = args.vcffile
output_dir = args.output_dir
bamfile = args.bamfile
n_thread = args.n_thread
max_dist = args.max_dist
min_support = args.min_support
min_mapq = args.min_mapq
flanking = args.flanking


import os
import pysam
from collections import defaultdict
from joblib import Parallel, delayed
from tqdm import tqdm

def eval_inv(vcffile,outdir):
    cmd = f'''python3 /data/maiziezhou_lab/CanLuo/long_reads_project/Read_Simulation/bin/Evaluation_tools/eval_inv.py \
        -callvcf {vcffile} \
        -o {outdir}'''
    os.system(cmd)
    return 


def load_vcf(vcffile, output_dir):
    var_list = []
    header =  []
    with open(vcffile,'r') as f:
        for line in f:
            if line[0]!='#':
                if 'SVTYPE=INV' in line:
                    data = line.split()
                    chrom = data[0]
                    pos = int(data[1])
                    end = int(data[7].split('END=')[1].split(';')[0])
                    var_list.append((chrom,pos, end, line))
            else:
                header.append(line)
    # print("eval original INV")
    # eval_inv(vcffile,output_dir + '/eval_og')
    return var_list,header

def merge_inv(var_list, max_dist):
    clusters = [[var_list[0]]]
    for i in range(1, len(var_list)):
        old_chrom,old_start, old_end , _ = clusters[-1][-1]
        new_chrom,new_start, new_end , _ = var_list[i]
        if (old_chrom == new_chrom) & ( abs(old_start - new_start) <= max_dist) & ( abs(old_end - new_end) <= max_dist):
            clusters[-1].append(var_list[i])
        else:
            clusters.append([var_list[i]])

    
    lines = []
    for vars in clusters:
        if len(vars)==1:
            lines.append(vars[0][-1])
        else:
            gt_hp1 = get_gt_votes(vars,hp = 1)
            gt_hp2 = get_gt_votes(vars,hp = 2)
            new_gt = str(gt_hp1)+'/'+str(gt_hp2)
            data = vars[0][-1].split()
            data[-1] = new_gt
            lines.append('\t'.join(data)+'\n')


    return clusters,lines

def get_gt_votes(var_list,hp):
    '''hp = 1 or 2 '''
    hp = hp-1
    assert hp in [0,1]

    votes=[]
    for var in var_list:
        vt = int(var[-1].split()[-1].split(':')[0].split('/')[hp])
        votes.append(vt)
    
    if sum(votes):
        return 1 
    else:
        return 0 
    

def write_vcf(new_lines,header, outfile, output_dir):

    with open(outfile,'w') as f:
        f.writelines(header+ new_lines)
 
    # print("eval merged INV")
    # eval_inv(outfile,output_dir + '/eval_merged')

    cmd = f"grep PASS {outfile} > {outfile}.passonly"
    os.system(cmd)
    # print("eval merged INV (PASS only)")
    # eval_inv(outfile+'.passonly',output_dir + '/eval_merged_passonly')
    return  

def calc_support(reads_forward, reads_reverse):
    
    reads_f = set(reads_forward)
    reads_r = set(reads_reverse)

    overlap_reads = reads_f & reads_r

    return len(overlap_reads)


def extract_reads_support_one_region(bamfile,chrom,bnd,flanking, min_mapq):

    samfile = pysam.AlignmentFile(bamfile)
    # reads_forward = defaultdict(list)
    # reads_reverse = defaultdict(list)
    reads_forward = []
    reads_reverse = []
    try:
        for read in samfile.fetch(chrom, bnd-flanking, bnd+flanking):
            ref_start = read.pos
            ref_end = read.reference_end 
            read_length = read.query_length
            cg_list = read.cigar 
            if cg_list[0][0] in [4,5]:
                read_start = cg_list[0][1]
            else:
                read_start = 0
            if cg_list[-1][0] in [4,5]:
                read_end = read_length -  cg_list[-1][1]
            else:
                read_end = read_length

            if read.mapq >= min_mapq:

                if read.is_reverse :
                    reads_reverse.append(read.qname)
                    # reads_reverse[read.qname].append((ref_start, ref_end, read_length - read_start, read_length - read_end, read_length, read.mapq ))
                else:
                    reads_forward.append(read.qname)
                    # reads_forward[read.qname].append((ref_start, ref_end, read_start, read_end, read_length, read.mapq ))
    except:
        pass
        
    # return reads_forward, reads_reverse
    return calc_support(reads_forward, reads_reverse)
    



def filter_inv(new_lines,bamfile, min_mapq, flanking, min_support, n_thread):
    line_info  = {}
    bnd_list = []
    for line in new_lines:
        data = line.split()
        chrom = data[0]
        start = int(data[1])
        end = int( data[7].split('END=')[1].split(';')[0])
        line_info[line] = (chrom, start, end )
        bnd_list.append((chrom, start))
        bnd_list.append((chrom, end ))

    n_support = Parallel(n_jobs=n_thread)(delayed( extract_reads_support_one_region )
                                          (bamfile,chrom,bnd,flanking, min_mapq) 
                                          for chrom,bnd in tqdm(bnd_list))

    dc_support = dict(zip(bnd_list, n_support))
    # print(dc_support)

    filtered_lines = []
    for key,val in line_info.items():
        chrom,start, end = val 
        n_support_start = dc_support[(chrom, start)]
        n_support_end = dc_support[(chrom, end)]

        if ( n_support_start >= min_support) & ( n_support_end >= min_support):
            filtered_lines.append(key)

    return filtered_lines


def write_filtered_vcf(filtered_lines,header, outfile, output_dir):

    with open(outfile,'w') as f:
        f.writelines(header+ filtered_lines)
 
    # print("eval filtered INV")
    # eval_inv(outfile,output_dir + '/eval_filtered')

    cmd = f"grep PASS {outfile} > {outfile}.passonly"
    os.system(cmd)
    # print("eval filtered INV (PASS only)")
    # eval_inv(outfile+'.passonly',output_dir + '/eval_filtered_passonly')
    return  















if not os.path.exists(output_dir):
    os.system("mkdir -p "+ output_dir)
            
var_list, header = load_vcf(vcffile, output_dir)
clusters, new_lines = merge_inv(var_list, max_dist)
outfile = output_dir+'/merged_inv.vcf'
write_vcf(new_lines ,header, outfile, output_dir)
        

filtered_lines = filter_inv(new_lines,bamfile, min_mapq, flanking, min_support, n_thread)
filtered_vcf = output_dir+'/INV_final.vcf'
write_filtered_vcf(filtered_lines,header, filtered_vcf, output_dir)