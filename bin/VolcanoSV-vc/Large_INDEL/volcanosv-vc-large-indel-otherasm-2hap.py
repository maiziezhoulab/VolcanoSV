from argparse import ArgumentParser
import argparse
parser = ArgumentParser(description="",
						usage='use "python3 %(prog)s --help" for more information',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# parser.add_argument('--fasta_pattern','-fa', help = "by chromosome fasta, replace chr number with *")
# parser.add_argument('--ref_pattern','-chrref', help = "by chromosome reference, replace chr number with *")

parser.add_argument('--hp1_fa','-hp1', help = "whole genome scale haplotype-resolved fasta file that includes both haplotype. Should include 'hp1' or 'hp2' in contig name to differentiate haplotype.")
parser.add_argument('--hp2_fa','-hp2', help = "whole genome scale haplotype-resolved fasta file that includes both haplotype. Should include 'hp1' or 'hp2' in contig name to differentiate haplotype.")

parser.add_argument('--output_dir','-o')
parser.add_argument('--data_type','-dtype',help='Hifi;CLR;ONT')
parser.add_argument('--bam_file','-bam', help = "reads bam file for reads signature extraction")
# parser.add_argument('--bam_file','-bam', help = "reads bam file for reads signature extraction; if both read_signature_dir and pre_cutesig are provided, you do not need to provide bam file")
parser.add_argument('--reference','-ref', help ="only needed when presig is not provided")


parser.add_argument('--read_signature_dir','-rdsig', help = "pre-extracted reads signatures; optional; if not provided, will generate one using read bamfile")
parser.add_argument('--pre_cutesig','-presig', help = "pre-extracted cutesv signature directory;optional; if not provided, will generate a new one")


parser.add_argument('--n_thread','-t',type = int, help = "number of threads",
					default = 11)
parser.add_argument('--n_thread_align','-ta',help = "number of threads for contig alignment",
					type = int, default = 10)
parser.add_argument('--mem_per_thread','-mempt', default = '768M',help = "Set maximum memory per thread for alignment; suffix K/M/G recognized; default = 768M")
parser.add_argument('--prefix','-px', help = "file prefix in the output folder", default = "Sample")
args = parser.parse_args()


global prefix
read_signature_dir = args.read_signature_dir
rbam_file = args.bam_file
reference=args.reference
pre_cutesig= args.pre_cutesig
# input_path = args.input_path
hp1_fa = args.hp1_fa
hp2_fa = args.hp2_fa
output_dir = args.output_dir
data_type = args.data_type
## optional
# chr_num = args.chr_num
chr_num = None
n_thread = args.n_thread
n_thread_align = args.n_thread_align
mem_per_thread = args.mem_per_thread
prefix=args.prefix

import os
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd
from subprocess import Popen
import subprocess
import pysam

global code_dir
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
os.system("mkdir -p "+output_dir)

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")




def create_fai(reference_fa):
	"""Create an FAI file for the reference genome if it doesn't exist."""
	fai_file = reference_fa + ".fai"

	# Check if FAI file exists
	if not os.path.exists(fai_file):
		print(f"FAI file not found. Creating FAI file for {reference_fa}.")
		try:
			# Create the .fai index file using samtools
			subprocess.run(["samtools", "faidx", reference_fa], check=True)
		except subprocess.CalledProcessError as e:
			print(f"Error in creating FAI file: {e}")
			return None

	return fai_file

def extract_contigs_from_fai(fai_file):
	"""Extract contigs and lengths from FAI file."""
	contigs = []
	with open(fai_file, 'r') as file:
		for line in file:
			parts = line.split('\t')
			contig_name = parts[0]
			contig_length = parts[1]
			contigs.append((contig_name, contig_length))
	return contigs

def generate_vcf_header(reference_fa, header_file, chr_num):
	# Step 1: Check and create .fai file if necessary
	fai_file = create_fai(reference_fa)

	# Step 2: Extract contig lengths and names
	contigs = extract_contigs_from_fai(fai_file)

	# Step 3: Generate VCF header
	"""Generate the VCF header from contig lengths."""
	vcf_header = "##fileformat=VCFv4.2\n"
	if chr_num is None:
		for contig_name, contig_length in contigs:
			vcf_header += f"##contig=<ID={contig_name},length={contig_length}>\n"
	else:
		target_contig_name = "chr"+str(chr_num)
		for contig_name, contig_length in contigs:
			if contig_name ==  target_contig_name:
				vcf_header += f"##contig=<ID={contig_name},length={contig_length}>\n"



	with open(code_dir +  "/header_info",'r') as f:
		header_info = f.read().replace("Sample", prefix)

	with open(header_file,'w') as f:
		f.write(vcf_header + header_info)

	return vcf_header



	

def split_reference(input_path,output_dir, chr_num):
	logger.info("load reference...")
	with open(input_path,'r') as f:
		s = f.read()
	chr_list = s.split('>')[1:]
	if not os.path.exists(output_dir):
		os.system("mkdir -p "+output_dir)
	for i in tqdm(range(22),desc = "write by chromosome ref"):
		if (chr_num is None) or ( int(i+1) == chr_num):
			path = output_dir+'/chr%d.fa'%(i+1)
			with open(path, 'w') as f:
				f.write('>'+chr_list[i])
			cmd = "samtools faidx "+path
			os.system(cmd)
	return 


def run_one_chrom(i,):
	ref_one_chr = ref_list[i]
	contig_one_chr = fasta_list[i]
	logger.info("process chromosome %d..."%(i+1))

	if read_signature_dir is None:
		cmd = "python3 "+code_dir+"/Raw_variant_call-otherasm.py \
		-contig %s  \
		-cbam %s -ref %s \
		-rbam %s  -o %s -dtype %s -t %d -chr %d"%(
			contig_one_chr,
			contig_bam,
			ref_one_chr,
			rbam_file,
			output_dir+"/chr%d/"%(i+1),
			data_type,
			n_thread_align,
			i+1
			)
	else:
		cmd = "python3 "+code_dir+"/Raw_variant_call-otherasm.py \
		-contig %s  \
		-cbam %s -ref %s \
		-sigd %s -o %s -dtype %s -t %d -chr %d"%(
			contig_one_chr,
			contig_bam,
			ref_one_chr,
			read_signature_dir,
			output_dir+"/chr%d/"%(i+1),
			data_type,
			n_thread_align,
			i+1
			)

	logger.info(cmd)
	Popen(cmd,shell = True).wait()
	logger.info("finish chromosome %d"%(i+1))
	return 


def read_vcf(vcf_path):
	with open(vcf_path,'r') as f:
		header = []
		body = []
		for line in f:
			if line[0]=='#':
				header.append(line)
			else:
				body.append(line)

	return header,body 


def phase_vcf(infile, outfile, header_file):

	body = []
	with open(infile,'r') as fin:
		for line in fin:
			if line[0]!='#':
				body.append(line)

	### process header
	with open(header_file,'r') as f:
		header = f.readlines()

	with open(outfile,'w') as fout:
		fout.writelines(header)
		for line in body:
			data = line.split()
			info = data[7].split('TIG_REGION=')[1].split(';')[0]
			hps = info.split(',')
			ps = hps[0].split('_')[0][2:]
			if data[-1] == '0/1':
				if 'hp1' in hps[0]:
					phased_gt = "1|0"
				else:
					phased_gt = "0|1"
			else:
				phased_gt = "1|1"
			data[7] = data[7]+';PS='+ps 
			data[-1] = phased_gt 
			line = '\t'.join(data)+'\n'
			fout.write(line)

# bam_path = output_dir+'/'+prefix+'.sorted.bam'

def align_contig_one_hap(contig_path, reference_path, bam_path, n_thread, mem_per_thread,):
	cmd = "minimap2 -a -x asm20 --cs -r2k -t %d \
		%s \
		%s \
			| samtools sort -@ %d -m %s > %s"%(n_thread,reference_path,contig_path,n_thread, mem_per_thread,bam_path )
	logger.info(cmd)
	Popen(cmd,shell=True).wait()

	cmd = "samtools index "+bam_path
	logger.info(cmd)
	Popen(cmd,shell=True).wait()


def align_contig(hp1_fa,hp2_fa, reference_path, bam_path, out_dir, n_thread, mem_per_thread,):

	hp1_bam = out_dir+"/hp1.bam"
	hp2_bam = out_dir+"/hp2.bam"
	align_contig_one_hap(hp1_fa, reference_path, hp1_bam, n_thread, mem_per_thread,)
	align_contig_one_hap(hp2_fa, reference_path, hp2_bam, n_thread, mem_per_thread,)
	cmd = f"samtools merge -u - {hp1_bam} {hp2_bam} | samtools sort -o {bam_path};samtools index {bam_path}"
	logger.info(cmd)
	Popen(cmd,shell=True).wait()


def split_contigs_by_chromosome(bam_file, hp1_fa, hp2_fa, output_folder):
	"""
	Splits contigs from a contig file by chromosomes (chr1-22) based on BAM file alignments.

	Args:
		bam_file (str): Path to the BAM file with contig alignments.
		contig_file (str): Path to the contig FASTA file.
		output_folder (str): Folder to store chromosome-specific FASTA files.
	"""
	# Ensure output folder exists
	os.makedirs(output_folder, exist_ok=True)

	# merge fasta
	contig_file = output_folder+"/merged.fa"

	cmd =f"cat {hp1_fa} {hp2_fa} > {contig_file}"
	logger.info(cmd)
	Popen(cmd,shell=True).wait()


	# Load BAM file
	bam = pysam.AlignmentFile(bam_file, "rb")

	# Filter chromosome names (chr1-22)
	valid_chromosomes = {f"chr{i}" for i in range(1, 23)}

	# Dictionary to store chromosome-specific contigs
	chromosome_contigs = {chrom: [] for chrom in valid_chromosomes}
	logger.info("Iterate over BAM alignments to group contigs by chromosome...")
	# Iterate over BAM alignments to group contigs by chromosome
	for read in bam.fetch(until_eof=True):
		if not read.is_unmapped and read.reference_name in valid_chromosomes:
			chromosome_contigs[read.reference_name].append(read.query_name)

	bam.close()

	# Load contig sequences
	logger.info("Load contig sequences...")
	contigs = {}
	with pysam.FastaFile(contig_file) as fasta:
		for contig_name in fasta.references:
			contigs[contig_name] = fasta.fetch(contig_name)

	# Write contigs grouped by chromosome into separate files
	for chrom, contig_names in chromosome_contigs.items():
		if contig_names:
			output_path = os.path.join(output_folder, f"{chrom}.fasta")
			with open(output_path, "w") as f:
				for contig_name in set(contig_names):
					if contig_name in contigs:
						f.write(f">{contig_name}\n{contigs[contig_name]}\n")

	print(f"Contigs split by chromosomes saved to {output_folder}.")







# set variables
header_file = output_dir+"/VCF_header"
ref_dir = output_dir + "/ref_by_chr/"
ref_list = [ref_dir + "/chr"+str(i+1)+".fa" for i in range(22)]
contig_bam = output_dir+'/'+prefix+'.sorted.bam'
contig_by_chr_folder = output_dir + "/contig_by_chr/"
fasta_list =  [contig_by_chr_folder+"/chr"+str(i+1)+".fasta" for i in range(22)]
raw_vcf = output_dir+"/raw_variants_wgs.vcf"
chr_para = " "
final_gonotyped_vcf = f"{output_dir}/variants_filtered_GT_corrected.vcf"
final_phased_vcf = f"{output_dir}/{prefix}_volcanosv_large_indel.vcf"

# create output folder
if not os.path.exists(output_dir):
	os.system("mkdir -p "+output_dir)

# align contig
logger.info(f"Aligning contig {hp1_fa} {hp2_fa} to reference {reference}...")
align_contig(hp1_fa, hp2_fa, reference , contig_bam, output_dir, n_thread, mem_per_thread,)

# split contig
logger.info("split assemblies by chromosome...")
split_contigs_by_chromosome(contig_bam, hp1_fa,hp2_fa, contig_by_chr_folder)

# generate VCF header
logger.info("generate VCF header...")
vcf_header = generate_vcf_header(reference, header_file, chr_num)

# split reference
logger.info("split reference by chromosome...")
split_reference(reference, ref_dir, chr_num)




# loop each chr to call SV
results = Parallel(n_jobs=n_thread)(delayed(run_one_chrom)(i) for i in tqdm(range(22)))

# collect all result 
logger.info("collect all chromosome result")
body_wgs = []
for i in range(22):
	vcf_path = output_dir+'/chr%d/final_vcf/volcano_variant_no_redundancy.vcf'%(i+1)
	header,body = read_vcf(vcf_path)
	body_wgs.extend(body)
with open(raw_vcf,'w') as f:
	f.writelines(header + body_wgs)

# filter and GT correct
logger.info("Filter VCF and correct genotype...")
if pre_cutesig is not None:
	cmd = f'''python3 {code_dir}/filter_GT_correction.py \
		-vcf {raw_vcf} \
		-presig {pre_cutesig} \
		-dtype {data_type} \
		-t {n_thread} {chr_para}'''
else:
	cmd = f'''python3 {code_dir}/filter_GT_correction.py \
		-vcf {raw_vcf} \
		-bam {rbam_file} \
		-ref {reference} \
		-dtype {data_type} \
		-t {n_thread} {chr_para}'''
print(cmd)
Popen(cmd, shell = True).wait()

# phase VCF
phase_vcf(final_gonotyped_vcf, final_phased_vcf,header_file)

