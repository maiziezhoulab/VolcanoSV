from argparse import ArgumentParser
import multiprocessing as mp
import pandas as pd
from subprocess import Popen
import os

parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--working_dir','-w')
parser.add_argument('--hap_name','-hap')
parser.add_argument('--final_contig_dir','-o')
args = parser.parse_args()
input_path =  args.input_path
working_dir =  args.working_dir
hap_name = args.hap_name
final_contig_dir = args.final_contig_dir
prefix = hap_name+'.asm'


def reformat_fasta(raw_fasta,hap_name):
	with open(raw_fasta,'r') as f:
		s = f.readlines()
	fw = open(raw_fasta,'w')
	idx = 0
	for line in s:
		if line[0]=='>':
			line = '>'+hap_name+'_'+str(idx)+'\n'
			idx+=1
		fw.write(line)
	fw.close()
	return raw_fasta




if not os.path.exists(working_dir):
	os.makedirs(working_dir)
if not os.path.exists(final_contig_dir ):
	os.makedirs(final_contig_dir)



global hifiasm_dir
hifiasm_dir=os.path.dirname(os.path.realpath(__file__))+'/hifiasm-0.14/'
# hifiasm_dir = "/data/maiziezhou_lab/CanLuo/Software/hifiasm-0.14/"
# run assembly
cmd = hifiasm_dir + "/hifiasm -o %s/%s -t 2 \
%s"%(working_dir,prefix,input_path)
Popen(cmd, shell= True).wait()


# convert gfa to fa
cmd = "awk '/^S/{print \">\"$2;print $3}' %s/%s.p_ctg.gfa > %s/%s.p_ctg.fa"%(working_dir,prefix,working_dir,prefix)
Popen(cmd, shell= True).wait()

cmd = "awk '/^S/{print \">\"$2;print $3}' %s/%s.p_utg.gfa > %s/%s.p_utg.fa"%(working_dir,prefix,working_dir,prefix)
Popen(cmd, shell= True).wait()

cmd = "awk '/^S/{print \">\"$2;print $3}' %s/%s.r_utg.gfa > %s/%s.r_utg.fa"%(working_dir,prefix,working_dir,prefix)
Popen(cmd, shell= True).wait()


# reformat fa
raw_fasta = "%s/%s.p_ctg.fa"%(working_dir,prefix)
reformat_fasta(raw_fasta,hap_name)

raw_fasta = "%s/%s.p_utg.fa"%(working_dir,prefix)
reformat_fasta(raw_fasta,hap_name)

raw_fasta = "%s/%s.r_utg.fa"%(working_dir,prefix)
reformat_fasta(raw_fasta,hap_name)

# concat fasta file to an overall fasta file
cmd = "cat %s/%s.p_ctg.fa >> %s/final_contigs.fa"%(working_dir,prefix,final_contig_dir)
Popen(cmd, shell= True).wait()
cmd = "cat %s/%s.p_utg.fa >> %s/final_contig.p_utg.fa"%(working_dir,prefix,final_contig_dir)
Popen(cmd, shell= True).wait()
cmd = "cat %s/%s.r_utg.fa >> %s/final_contig.r_utg.fa"%(working_dir,prefix,final_contig_dir)
Popen(cmd, shell= True).wait()


# remove intermediate files

cmd = "rm %s/%s*"%(working_dir,prefix)
Popen(cmd, shell= True).wait()









