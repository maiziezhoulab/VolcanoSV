from joblib import Parallel, delayed
import os
from argparse import ArgumentParser
import pandas as pd
from subprocess import Popen
import time
from tqdm import tqdm
global tool_name, input_dir, final_contig_dir, fasta_by_hap_dir,working_dir, code_dir,threads_asm,log_file_path
parser = ArgumentParser(description="Author: can.luo@vanderbilt.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_dir','-i')
parser.add_argument('--working_dir','-w')
parser.add_argument('--final_contig_dir','-o')
# parser.add_argument('--pb_info_path','-p')
parser.add_argument('--tool_name','-tool',help = "shasta;flye",default = 'flye')
parser.add_argument('--data_type','-dtype',help="CLR;ONT")
parser.add_argument('--n_thread','-t',type=int )
parser.add_argument('--threads_asm','-ta',type = int)
args = parser.parse_args()
input_dir =  args.input_dir
working_dir =  args.working_dir
final_contig_dir = args.final_contig_dir
# pb_info_path = args.pb_info_path
tool_name = args.tool_name
data_type = args.data_type
n_thread = args.n_thread
threads_asm = args.threads_asm
import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

code_dir=os.path.dirname(os.path.realpath(__file__))+'/'
fasta_by_hap_dir = final_contig_dir+'/fasta_by_hap/'

# df = pd.read_csv(pb_info_path)

os.system("mkdir -p " + working_dir)
os.system("mkdir -p " + final_contig_dir)
os.system("mkdir -p " + fasta_by_hap_dir)

log_file_path = final_contig_dir+'/log.txt'
os.system("rm "+log_file_path)
os.system("touch "+log_file_path)


def assemble_one_hp(hap_name):
	temp_dir = working_dir+'/'+hap_name+'/'
	fastq_path = input_dir+'/'+hap_name+'.fastq'
	output_path = fasta_by_hap_dir + '/' + hap_name+'.fa'
	cmd = "python3 %s/%s_single_task.py -i %s \
	-o %s \
	-temp %s -t %d -dtype %s"%(code_dir,tool_name,fastq_path,output_path, temp_dir,threads_asm, data_type)
	logger.info(cmd)
	Popen(cmd, shell= True).wait()

	# write log file
	# log_file_path = final_contig_dir+'/log.txt'
	cmd = "echo '%s' >> %s"%(hap_name,log_file_path)
	Popen(cmd, shell= True).wait()
	return output_path


hap_name_list = [ fq.split('.')[0]  for fq in os.listdir(input_dir)]
# hap_name_list = []
# for i in range(df.shape[0]):
# 	pb_name = 'PS%d_%d_%d'%(df['name'][i],df['start'][i],df['end'][i])
# 	# hap_name = pb_name+'_hp1'
# 	# fastq_path = input_dir+'/'+hap_name+'.fastq'
# 	hap_name_list.append(pb_name+'_hp1')
# 	hap_name_list.append(pb_name+'_hp2')

order_log_path = final_contig_dir+'/order.txt'
with open(order_log_path,'w') as f:
	f.write('\n'.join(hap_name_list)+'\n')

if n_thread==1:
	contig_paths = [assemble_one_hp(hap_name) for hap_name in tqdm(hap_name_list)]
else:


	contig_paths = Parallel(n_jobs=n_thread)(delayed(assemble_one_hp)(hap_name) for hap_name in tqdm(hap_name_list))

## merge files

def reformat_and_write_fasta(fin,fout,hap_name):
	cnt = -1
	for line in fin:
		if line[0]=='>':
			cnt+=1
			fout.write(">%s_%d\n"%(hap_name, cnt))
		else:
			fout.write(line)
	return


final_contig_path = final_contig_dir+'/final_contigs.fa'
fout = open(final_contig_path,'w')

for hap_name in hap_name_list:
	contig_path = fasta_by_hap_dir+'/%s.fa'%hap_name
	if os.path.exists(contig_path):
		fin = open(contig_path)
		reformat_and_write_fasta(fin,fout,hap_name)
	else:
		logger.info(hap_name + " fail to assemble")

fout.close()
