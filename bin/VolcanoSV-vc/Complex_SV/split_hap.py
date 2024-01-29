from argparse import ArgumentParser
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_file','-i')
parser.add_argument('--output_dir','-o')
args = parser.parse_args()
input_file = args.input_file
output_dir = args.output_dir

fhp1 = open(output_dir +"/wgs_hp1.fa",'w')
fhp2 = open(output_dir +"/wgs_hp2.fa",'w')
# fasta_list = [ input_dir + "/chr%d/assembly/final_contigs/final_contigs.fa"%i for i in range(1,23)]
# for fasta_file in fasta_list:
with open(input_file ,'r') as f0:
	for line in f0:
		if 'hp1' in line:
			f = fhp1
		if 'hp2' in line:
			f = fhp2 
		f.write(line)

fhp1.close()
fhp2.close()














