from argparse import ArgumentParser
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_dir','-i')
parser.add_argument('--prefix','-px',default='DipPAV_variants')
args = parser.parse_args()
input_dir = args.input_dir



body = []
for i in range(22):
	chrom = "chr%d"%(i+1)
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






