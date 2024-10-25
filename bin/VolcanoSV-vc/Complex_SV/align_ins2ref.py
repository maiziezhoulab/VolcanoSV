import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--reference','-ref')
parser.add_argument('--n_thread','-t', type = int, default = 50 )
parser.add_argument('--datatype','-d', choices=['Hifi','CLR','ONT'])
parser.add_argument('--min_mapq','-q', type = int, default= 0, help = 'min_mapq')
parser.add_argument('--min_size_sim','-P', type = float, default= 0.7, help = 'min_size_sim')
parser.add_argument('--max_shift','-r', type = int, default = 300, help = 'max_shift')
parser.add_argument('--max_shift_ratio','-rr', type = float, default= 0.3, help = 'max_shift_ratio' )


args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
reference = args.reference
n_thread = args.n_thread
datatype = args.datatype
min_mapq = args.min_mapq
min_size_sim = args.min_size_sim 
max_shift = args.max_shift
max_shift_ratio = args.max_shift_ratio

dc = {
    "Hifi":'hifi',
    "CLR":'pb',
    "ONT":'ont'
}
datatype = dc[datatype]
# reference = "/data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa"

import os
from subprocess import Popen
import pysam
import numpy as np 
from collections import defaultdict
from tqdm import tqdm
if not os.path.exists(output_dir):
    os.system("mkdir -p "+ output_dir)

def vcf2fasta(vcffile,ofile):
    ins_dc ={}
    with open(ofile,'w') as fo:
        with open(vcffile,'r') as fin:
            for line in fin:
                if line[0]!='#':
                    if 'SVTYPE=INS' in line:
                        data = line.split()
                        ins_seq = data[4]
                        svid = data[2]
                        chrom, pos = data[0], data[1]
                        pos = int(pos)
                        gt = data[-1].split(':')[0]
                        fo.write('>'+svid+'\n')
                        fo.write(ins_seq+'\n')
                        ins_dc[svid] = (chrom, pos, len(ins_seq),gt,line)
    
    return ins_dc

def align_insseq(fastafile, reference, datatype, bamfile):
    ''' pb hifi ont '''

    cmd = f"minimap2  -a -x map-{datatype} --cs -r2k -t {n_thread } {reference} {fastafile} | samtools sort -@ 20 -m 4G > {bamfile}; samtools index {bamfile} "
    print(cmd)
    Popen(cmd, shell = True).wait()

    return 

def extract_baminfo(bamfile):
    samfile = pysam.AlignmentFile(bamfile)
    dc = defaultdict(list)

    for read in samfile.fetch():
        info = (read.reference_name, read.pos,read.reference_end, read.mapq)
        dc[read.qname].append(info)
    return dc 

def is_dup(info_aln, info_var, min_mapq, min_size_sim, max_shift, max_shift_ratio):
    ''' max_shift_ratio: shift / svlen  '''
    chrom_aln, start_aln, end_aln, mapq_aln = info_aln
    chrom_var, bnd_var, svlen_var, _, line = info_var 
    svlen_aln = end_aln - start_aln 

    size_sim = min(svlen_aln, svlen_var)/ max(svlen_var, svlen_aln)
    shift = min ( abs( start_aln - bnd_var), abs(end_aln - bnd_var))
    shift_ratio = shift / svlen_var

    if (chrom_aln == chrom_var) & ( mapq_aln >= min_mapq) & (size_sim >= min_size_sim ) & ( shift <= max_shift) & ( shift_ratio <= max_shift_ratio):
        # is dup
        return 1, size_sim, shift, shift_ratio, line
    else:
        # is not dup
        return 0, size_sim, shift, shift_ratio, line

def is_dup_one_var(aln_list, info_var, min_mapq, min_size_sim, max_shift, max_shift_ratio):

    dec_list = [ is_dup(info_aln, info_var, min_mapq, min_size_sim, max_shift , max_shift_ratio)  for info_aln in  aln_list ]

    dup_flag = sum([dec[0] for dec in dec_list])
    assert dup_flag>=0
    if dup_flag == 0:
        return 0, ()
    elif dup_flag==1:
        for i in range(len(dec_list)):
            if dec_list[i][0]:
                dup_aln = aln_list[i]
                break 
        return 1, dup_aln
    elif dup_flag>1:
        candidate_list = []
        candidate_info_list = []
        
        for i in range(len(dec_list)):
            if dec_list[i][0]:
                _, size_sim, shift, shift_ratio, line = dec_list[i]
                aln_info = aln_list[i]
                candidate_list.append(aln_info)
                candidate_info_list.append([size_sim,  - shift,  - shift_ratio])

        candidate_info_array = np.array(candidate_info_list)
        
        # Normalize each column
        candidate_info_array_norm = (candidate_info_array - np.mean(candidate_info_array, axis=0)) / (np.std(candidate_info_array, axis=0)+0.0001)
        cand_scores = candidate_info_array_norm.sum(1)

        best_dup = candidate_list[np.argmax(cand_scores)]
        return 1, best_dup
    else:
        print(f"abnormal dup flag value {dup_flag}; allowed value: >=0, integer")
        exit()


def extract_dup(dc_aln, dc_var, min_mapq, min_size_sim, max_shift, max_shift_ratio):
    dc_final = {}
    aln_cnt = 0
    for svid in tqdm(dc_aln, desc = "process DUP candidates"):
        aln_list = dc_aln[svid]
        info_var = dc_var[svid]
        aln_cnt+=len(aln_list)
        dec, dup_aln = is_dup_one_var(aln_list, info_var,min_mapq, min_size_sim, max_shift, max_shift_ratio )
        if dec:
            
            dc_final[svid] = dup_aln 

    print(f"# Original INS: {len(dc_var)}")
    print(f"# INS that are aligned: {len(dc_aln)}")
    print(f"# alignment records: {aln_cnt}")
    print(f"# INS that are DUP: {len(dc_final)} ( {np.round((len(dc_final)/ len(dc_var) * 100), 2)} % )")
    return dc_final

def write_dup(ofile,dc_var, dc_final, og_vcffile):
    with open(ofile,'w') as f:
        with open(og_vcffile,'r') as fin:
            for line in fin:
                if line[0]=='#':
                    f.write(line)
                elif 'SVTYPE=DUP' in line:
                    f.write(line)


        cnt = 0

        for svid in dc_final:
            chrom, start, end , mapq = dc_final[svid]
            _,_,_,gt,line = dc_var[svid ]
            reads = line.split('TIG_REGION=')[1].split(';')[0]
            cnt+=1
            line = '\t'.join([chrom,str(start),\
                              'volcanosv.DUP.recover.%d'%cnt,'.','<DUP>',\
                              '20', 'PASS','SVTYPE=DUP;END=%d;SVLEN=%d;READS=%s'%(end, end-start,reads),'GT',gt,])+'\n'

            f.write(line)
    








        










fastafile = output_dir+'/ins.fa'
bamfile = output_dir+'/ins_aligned.bam'
dupfile = output_dir+'/DUP_recovered_from_INS.vcf'
dc_var = vcf2fasta(input_path, fastafile)
align_insseq(fastafile, reference, datatype, bamfile)

dc_aln = extract_baminfo(bamfile)

dc_dup = extract_dup(dc_aln, dc_var, min_mapq, min_size_sim, max_shift, max_shift_ratio)

write_dup(dupfile,dc_var, dc_dup, input_path)


# cmd = f'''python3 /data/maiziezhou_lab/CanLuo/long_reads_project_lio/Read_Simulation/bin/Evaluation_tools/eval_dup.py -callvcf {dupfile} -o {output_dir}/eval '''
# os.system(cmd)


    






