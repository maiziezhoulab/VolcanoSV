import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_dir','-i')
parser.add_argument('--indelvcf','-vcf')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--reference','-ref')
parser.add_argument('--datatype','-d', choices=['CCS','CLR','ONT'])
parser.add_argument('--out_dir','-o')
parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_dir = args.input_dir
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


def merge_fasta(input_dir,outdir):
   # if datatype == "CCS":
   #    fasta_list =  [input_dir+"/chr"+str(i+1)+"/assembly/final_contigs/final_contig.p_ctg.fa" for i in range(22)]
   # else:
   
   fasta_list =  [input_dir+"/chr"+str(i+1)+"/assembly/final_contigs/final_contigs.fa" for i in range(22)]

   if not os.path.exists(outdir):
      os.system("mkdir -p " + outdir)
   hp1file = outdir+"/hp1.fa"
   hp2file = outdir+"/hp2.fa"
   fhp1 = open(hp1file,'w')
   fhp2 = open(hp2file,'w')
   for fasta in fasta_list:
      with open(fasta,'r') as f:
         for line in f:
            if line[0] == '>':
               if "hp1" in line:
                  fw = fhp1 
               else:
                  fw = fhp2 
            fw.write(line)
   fhp1.close()
   fhp2.close()
   return 

hp1fa = out_dir+"/hp1.fa"
hp2fa = out_dir+"/hp2.fa"
raw_dir = out_dir+"/Raw_Detection/"

logger.info("-------------------------------Merge contigs")
merge_fasta(input_dir,out_dir)

logger.info("-------------------------------Extract raw complex SV")

os.system("mkdir -p " + raw_dir)
cmd = f'''minimap2 -a -x asm10 --cs -r2k -t {n_thread} \
   {reference} \
   {hp1fa} \
   | samtools sort -@ {n_thread} -m 4G > {raw_dir}/assembly_hp1.bam
samtools index -@ {n_thread} {raw_dir}/assembly_hp1.bam'''
Popen(cmd, shell = True).wait()

cmd = f'''minimap2 -a -x asm10 --cs -r2k -t {n_thread} \
   {reference} \
   {hp2fa} \
   | samtools sort -@ {n_thread} -m 4G > {raw_dir}/assembly_hp2.bam
samtools index -@ {n_thread} {raw_dir}/assembly_hp2.bam'''
Popen(cmd, shell = True).wait()

cmd = f'''python3 {code_dir}/svim-asm-1.0.2/src/svim_asm/svim-asm diploid {raw_dir}/ \
    {raw_dir}/assembly_hp1.bam  {raw_dir}/assembly_hp2.bam {reference}'''
Popen(cmd, shell = True).wait()




logger.info("-------------------------------DUP detection")
cmd = f"python3 {code_dir}/align_ins2ref.py -i {indelvcf} -o {out_dir}/DUP -d {datatype} \
   -ref {reference}"
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


logger.info("-------------------------------Merge DUP TRA INV")
cmd = f"cat {out_dir}/DUP/DUP_final.vcf {out_dir}/TRA/TRA_final.vcf {out_dir}/INV/INV_final.vcf > {out_dir}/volcanosv_complex_SV.vcf "
Popen(cmd, shell = True).wait()

 








