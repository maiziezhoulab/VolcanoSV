import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--hp1fa','-hp1')
parser.add_argument('--hp2fa','-hp2')
parser.add_argument('--indelvcf','-vcf')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--reference','-ref')
parser.add_argument('--datatype','-d', choices=['CCS','CLR','ONT'])
parser.add_argument('--out_dir','-o')
parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
hp1fa = args.hp1fa
hp2fa = args.hp2fa
indelvcf = args.indelvcf
bamfile = args.bamfile
out_dir = args.out_dir
reference = args.reference
n_thread = args.n_thread
datatype = args.datatype


import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")




import os
from subprocess import Popen
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

logger.info("-------------------------------Extract raw complex SV")
raw_dir = out_dir+"/Raw_Detection/"
os.system("mkdir -p " + raw_dir)
cmd = f'''minimap2 -a -x asm10 --cs -r2k -t {n_thread} \
   ${reference} \
   ${hp1fa} \
   | samtools sort -@ {n_thread} -m 4G > ${raw_dir}/assembly_hp1.bam
samtools index -@ {n_thread} ${raw_dir}/assembly_hp1.bam'''
Popen(cmd, shell = True).wait()

cmd = f'''minimap2 -a -x asm10 --cs -r2k -t {n_thread} \
   ${reference} \
   ${hp2fa} \
   | samtools sort -@ {n_thread} -m 4G > ${raw_dir}/assembly_hp2.bam
samtools index -@ {n_thread} ${raw_dir}/assembly_hp2.bam'''

Popen(cmd, shell = True).wait()

cmd = f'''python3 {code_dir}/svim-asm-1.0.2/src/svim_asm/svim-asm diploid {raw_dir}/ \
    {raw_dir}/assembly_hp1.bam  {raw_dir}/assembly_hp2.bam {reference}'''
Popen(cmd, shell = True).wait()




logger.info("-------------------------------DUP detection")
cmd = f"python3 {code_dir}/align_ins2ref.py -i {indelvcf} -o {out_dir}/DUP -d {datatype}"
Popen(cmd, shell = True).wait()
cmd = f"cat {raw_dir}/variants.vcf |grep SVTYPE=DUP > {raw_dir}/variants_dup.vcf"
Popen(cmd, shell = True).wait()
cmd = f"cat {out_dir}/DUP/DUP_recovered_from_INS.vcf {raw_dir}/variants_dup.vcf > {out_dir}/DUP/DUP_final.vcf"
Popen(cmd, shell = True).wait()


logger.info("-------------------------------TRA detection")

cmd = f"python3 {code_dir}/filter_tra.py -vcf {raw_dir}/variants.vcf -o {out_dir}/TRA/ -bam {bamfile}"
Popen(cmd, shell = True).wait()

logger.info("-------------------------------INV detection")
cmd = f"python3 {code_dir}/filter_inv.py -vcf {raw_dir}/variants.vcf -o {out_dir}/INV/ -bam {bamfile}"
Popen(cmd, shell = True).wait()


logger.info("-------------------------------INV detection")
cmd = "cat {out_dir}/DUP/DUP_final.vcf {out_dir}/TRA/TRA_final.vcf {out_dir}/INV/INV_final.vcf > {out_dir}/complex_SV.vcf "
Popen(cmd, shell = True).wait()

 








