from argparse import ArgumentParser
from subprocess import Popen
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
parser.add_argument('--data_type','-dtype',help = "CLR;ONT")
parser.add_argument('--temp_dir','-temp')

parser.add_argument('--n_thread','-t', type = int)


args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
temp_dir = args.temp_dir
n_thread = args.n_thread
data_type = args.data_type

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

### remove duplicate id in fastq

def remove_duplicate(input_path,output_path):

	fw = open(output_path,'w')
	name_dc = {}
	write_or_not = []
	cnt = 0
	for line in open(input_path):
		if cnt%4==0:
			name = line[1:-1]
			if name not in name_dc:
				name_dc[name] = 1 
				write_or_not.append(1)
			else:
				write_or_not.append(0)
		read_num = cnt//4
		decision = write_or_not[read_num]
		if decision==1:
			fw.write(line)
		cnt+=1
	fw.close()
	return 


## run asm
def assemble(input_path,output_dir, data_type, n_thread):
	if data_type =='ONT':

		cmd = "flye --nano-raw %s --out-dir %s --threads %d"%(input_path,output_dir,n_thread)
	elif data_type =='CLR':

		cmd = "flye --pacbio-raw %s --out-dir %s --threads %d"%(input_path,output_dir,n_thread)
	else:
		logger.info("only support ONT and CLR")
		return 

	Popen(cmd, shell = True).wait()
	return 



os.system("mkdir -p " + temp_dir )
new_fq_path = temp_dir+'/hap.fq'
remove_duplicate(input_path, new_fq_path)
assemble(new_fq_path,temp_dir, data_type, n_thread)

cmd = "cp %s/assembly.fasta %s"%(temp_dir,output_path)
Popen(cmd, shell = True).wait()

cmd = "rm -r "+temp_dir
Popen(cmd, shell = True).wait()







