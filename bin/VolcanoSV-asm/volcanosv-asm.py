import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam_file','-bam')
parser.add_argument('--output_dir','-o')
parser.add_argument('--reference','-ref')
parser.add_argument('--chrnum','-chr', type = int, choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22] )
parser.add_argument('--assembler','-asm', choices = ['wtdbg2','canu','miniasm','shasta','nextdenovo','hifiasm','hicanu','flye'])
parser.add_argument('--data_type','-dtype', choices = ['CLR','ONT','Hifi'])
parser.add_argument('--pacbio_subtype','-pb', choices = ['CLR-rs','CLR-sq'], help = "must provide when using wtdbg2 on CLR data")
parser.add_argument('--shasta_ont_config','-shacon', choices = ['Nanopore-OldGuppy-Sep2020'], help = "must provide when using shasta")
parser.add_argument('--n_thread','-t', required = True,type = int, help = "number of threads", default =10 )
parser.add_argument('--prefix','-px', help = "file prefix in the output folder", default = "Sample")
parser.add_argument('--clean','-cl', action='store_true')
args = parser.parse_args()



inbam=args.bam_file
prefix=args.prefix
reference=args.reference
chrnum = args.chrnum
out_dir= args.output_dir + "/chr"+str(chrnum)+"/"
dtype=args.data_type
pacbio_subtype = args.pacbio_subtype
n_thread = args.n_thread
assembler = args.assembler
shasta_ont_config = args.shasta_ont_config
tkc=n_thread
tka=n_thread
tA=n_thread
ta=8
# canu and hicanu has minimum threads requirement
clean = args.clean

# prefix = "sample"
import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

from subprocess import Popen
from General_Assembly_Workflow import run_assembly
from write_fastq_asm_general import write_fastqs_general
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



# write fastq

logger.info("write fastq...")


write_fastqs_general(f"{out_dir}/{prefix}.bam", f"{out_dir}/assembly/", f"{out_dir}/kmer_assign/", None, None)



logger.info("By haplotype assembly...")
fq_dir =  f"{out_dir}/assembly/fastq_by_hap/"
run_assembly(None, [fq_dir], out_dir+"/assembly", [assembler], 
                 n_thread, ta, clean, dtype, shasta_ont_config,
                  pacbio_subtype, prefix )



####### remove bam, fastq file####
if clean:
    cmd = f'''rm -r {out_dir}/assembly/fastq_by_hap
    rm {out_dir}/phasing_result/{prefix}_phased.bam*
    rm {out_dir}/{prefix}.bam*'''
    Popen(cmd, shell = True).wait()


