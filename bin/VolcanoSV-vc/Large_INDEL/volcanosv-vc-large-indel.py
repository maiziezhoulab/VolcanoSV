from argparse import ArgumentParser

parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')

# parser.add_argument('--fasta_pattern','-fa', help = "by chromosome fasta, replace chr number with *")
# parser.add_argument('--ref_pattern','-chrref', help = "by chromosome reference, replace chr number with *")

parser.add_argument('--input_dir','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--data_type','-dtype',help='CCS;CLR;ONT')

parser.add_argument('--rbam_file','-rbam', help = "reads bam file for reads signature extraction; if both read_sognature_dir and pre_cutesig are provided, you do not need to provide bam file")
parser.add_argument('--reference','-ref', help ="only needed when presig is not provided")


parser.add_argument('--read_signature_dir','-rdsig', help = "pre-extracted reads signatures; optional; if not provided, will generate one using read bamfile")
parser.add_argument('--pre_cutesig','-presig', help = "pre-extracted cutesv signature directory;optional; if not provided, will generate a new one")

## optional
parser.add_argument('--n_thread','-t',type = int, default = 11)
parser.add_argument('--n_thread_align','-ta',type = int, default = 10)
parser.add_argument('--mem_per_thread','-mempt', default = '768M',help = "Set maximum memory per thread for alignment; suffix K/M/G recognized; default = 768M")

args = parser.parse_args()

# fasta_pattern = args.fasta_pattern
# ref_pattern = args.ref_pattern
read_signature_dir = args.read_signature_dir
rbam_file = args.rbam_file
reference=args.reference
pre_cutesig= args.pre_cutesig
input_dir = args.input_dir
output_dir = args.output_dir
data_type = args.data_type
## optional
n_thread = args.n_thread
n_thread_align = args.n_thread_align
mem_per_thread = args.mem_per_thread



import os
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd
from subprocess import Popen


code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
os.system("mkdir -p "+output_dir)

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


def split_reference(input_path,output_dir):
    logger.info("load reference...")
    with open(input_path,'r') as f:
        s = f.read()
    chr_list = s.split('>')[1:]
    if not os.path.exists(output_dir):
        os.system("mkdir -p "+output_dir)
    for i in tqdm(range(22),desc = "write by chromosome ref"):
        path = output_dir+'/chr%d.fa'%(i+1)
        with open(path, 'w') as f:
            f.write('>'+chr_list[i])
        cmd = "samtools faidx "+path
        os.system(cmd)
    return 


logger.info("split reference by chromosome...")
ref_dir = output_dir+"/ref_by_chr/"
split_reference(reference, ref_dir)

if data_type == "CCS":
	fasta_list =  [input_dir+"/chr"+str(i+1)+"/assembly/final_contigs/final_contig.p_ctg.fa" for i in range(22)]
else:
	fasta_list =  [input_dir+"/chr"+str(i+1)+"/assembly/final_contigs/final_contigs.fa" for i in range(22)]
ref_list = [ref_dir + "/chr"+str(i+1)+"/" for i in range(22)]

def run_one_chrom(i):
	ref_one_chr = ref_list[i]
	logger.info("process chromosome %d..."%(i+1))

	if read_signature_dir is None:
		cmd = "python3 "+code_dir+"/Raw_variant_call.py \
		-contig %s -ref %s \
		-rbam %s  -o %s -dtype %s -t %d -chr %d"%(
			fasta_list[i],
			ref_one_chr,
			rbam_file,
			output_dir+'/chr%d'%(i+1),
			data_type,
			n_thread_align,
			i+1
			)
	else:
		cmd = "python3 "+code_dir+"/Raw_variant_call.py \
		-contig %s -ref %s \
		-sigd %s, -o %s -dtype %s -t %d -chr %d"%(
			fasta_list[i],
			ref_one_chr,
			read_signature_dir,
			output_dir+'/chr%d'%(i+1),
			data_type,
			n_thread_align,
			i+1
			)

	logger.info(cmd)
	Popen(cmd,shell = True).wait()
	logger.info("finish chromosome %d"%(i+1))
	return 


results = Parallel(n_jobs=n_thread)(delayed(run_one_chrom)(i) for i in tqdm(range(22)))



### collect all result to be one 

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

logger.info("collect all chromosome result")
body_wgs = []
for i in range(22):
	vcf_path = output_dir+'/chr%d/final_vcf/volcano_variant_no_redundancy.vcf'%(i+1)
	header,body = read_vcf(vcf_path)
	body_wgs.extend(body)

with open(output_dir+"/variants.vcf",'w') as f:
	f.writelines(header + body_wgs)


logger.info("Filter VCF and correct genotype...")
if pre_cutesig is not None:
	cmd = f'''python3 {code_dir}/filter_GT_correction.py \
		-vcf {output_dir}/variants.vcf \
		-presig {pre_cutesig} \
		-dtype {data_type} \
		-t {n_thread}'''
else:
	cmd = f'''python3 {code_dir}/filter_GT_correction.py \
		-vcf {output_dir}/variants.vcf \
		-bam {rbam_file} \
		-ref {reference} \
		-dtype {data_type} \
		-t {n_thread}'''
Popen(cmd, shell = True).wait()







