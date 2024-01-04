import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--inbam','-i')
parser.add_argument('--out_dir','-o')
parser.add_argument('--reference','-r')
parser.add_argument('--n_thread','-t', type = int, default =10 )
parser.add_argument('--chrnum','-chr', type = int, choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22] )
parser.add_argument('--dtype','-d', choices = ['CCS','CLR','ONT'])
parser.add_argument('--prefix','-px', default = "Sample")
args = parser.parse_args()



inbam=args.inbam
prefix=args.prefix
reference=args.reference
chrnum = args.chrnum
out_dir= args.out_dir + "/chr"+str(chrnum)+"/"
dtype=args.dtype
n_thread = args.n_thread

tkc=n_thread
tka=n_thread
tA=n_thread
ta=2


import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

from subprocess import Popen


import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'





cmd = f"mkdir -p {out_dir}"
Popen(cmd, shell = True).wait()


logger.info("Extract by chromosome bam...")
cmd = f'''samtools view -b {inbam} chr{chrnum} >  {out_dir}/{prefix}.bam;
samtools index  {out_dir}/{prefix}.bam'''
Popen(cmd, shell = True).wait()




logger.info("Longshot phasing...")
cmd = f"mkdir -p {out_dir}/phasing_result"
Popen(cmd, shell = True).wait()

cmd = f'''longshot --bam {out_dir}/{prefix}.bam \
--ref {reference} \
--out {out_dir}/phasing_result/{prefix}_phased.vcf \
-O {out_dir}/phasing_result/{prefix}_phased.bam -F
samtools index {out_dir}/phasing_result/{prefix}_phased.bam'''
Popen(cmd, shell = True).wait()


logger.info("k-mer based reads partition...")
cmd = f'''python3 {code_dir}/unphased_reads_assignment_kmer_norm.py \
    --phased_bam {out_dir}/phasing_result/{prefix}_phased.bam \
    --output_dir {out_dir}/kmer_assign/ \
    -k 12 \
    --threads_kc {tkc} \
    --threads_asg {tka} \
    --significance_level 0.1'''
Popen(cmd, shell = True).wait()

####### remove seq file,hash file####
cmd = f'''rm {out_dir}/kmer_assign/*.seq
rm {out_dir}/kmer_assign/*.hash
rm {out_dir}/kmer_assign/*.name
rm -r {out_dir}/kmer_assign/kmer_db_by_hap
rm -r {out_dir}/kmer_assign/up_est_hash'''
Popen(cmd, shell = True).wait()
######


logger.info("By haplotype assembly...")
cmd = f'''python3 {code_dir}/assembly.py \
    -bam {out_dir}/{prefix}.bam  \
    -dc_og {out_dir}/kmer_assign/phasing_info/read_hp_og.p \
    -dc_up {out_dir}/kmer_assign/unphased_reads_assignment_dippav_norm_asg.p \
    -pinf {out_dir}/kmer_assign/phasing_info/pb_info.csv \
    -o {out_dir}/assembly/ \
    -dtype {dtype} -t {tA} -ta {ta}'''
Popen(cmd, shell = True).wait()



