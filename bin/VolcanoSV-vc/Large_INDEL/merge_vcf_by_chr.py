from argparse import ArgumentParser
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_dir','-i')
parser.add_argument('--prefix','-px',default='DipPAV_variants')
parser.add_argument('--chrom_list','-cl', required=True, help = "text file with one chromosome/contig name per line")
args = parser.parse_args()
input_dir = args.input_dir


def get_chrom_list(chrom_list_file):
	with open(chrom_list_file, 'r') as f:
		return [line.strip() for line in f if line.strip()]


chrom_list = get_chrom_list(args.chrom_list)

body = []
for chrom in chrom_list:
	vcf_path = input_dir+'/'+chrom+'/variant_call/final_vcf/volcano_variant_no_redundancy.vcf'
	header = []
	with open(vcf_path,'r') as f:
		for line in f:
			if line[0]=='#':
				header.append(line)
			else:
				body.append(line)

with open(input_dir+'/'+prefix+'.vcf','w') as f:
	f.writelines(header+body)






