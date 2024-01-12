from argparse import ArgumentParser
import multiprocessing as mp
import pandas as pd
from subprocess import Popen
import os
import pandas as pd
import os
import time
from multiprocessing import Pool,cpu_count,active_children,Manager

parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_dir','-i')
parser.add_argument('--working_dir','-w')
parser.add_argument('--final_contig_dir','-o')
parser.add_argument('--pb_info_path','-p')
parser.add_argument('--num_of_threads','-t',type=int,default = 10)
parser.add_argument('--preserve_intermediate_file','-pr', action='store_true')
args = parser.parse_args()
input_dir =  args.input_dir
working_dir =  args.working_dir
final_contig_dir = args.final_contig_dir
pb_info_path = args.pb_info_path
num_of_threads = args.num_of_threads
preserve_intermediate_file = args.preserve_intermediate_file
if preserve_intermediate_file:
	preserve_option = '-pr'
else:
	preserve_option = ''

import os
global code_dir
code_dir=os.path.dirname(os.path.realpath(__file__))+'/'
df = pd.read_csv(pb_info_path)



# need to delete old contig files
if not os.path.exists(final_contig_dir ):
	os.makedirs(final_contig_dir)
else:
	cmd = "rm -r "+final_contig_dir
	Popen(cmd, shell= True).wait()
	os.makedirs(final_contig_dir)


os.system("touch "+final_contig_dir+'/log.txt')
# Step 1: Init multiprocessing.Pool()
# pool = mp.Pool(mp.cpu_count())
# pool = mp.Pool(3)
print("cpu_count: ",mp.cpu_count())
with open("cpu.log",'w') as f:
	f.write("cpu_count: %d"%mp.cpu_count())
# Step 2: `pool.apply` the `howmany_within_range()`
def assembly_one_hp(fastq_path,working_dir,final_contig_dir,hap_name,code_dir):
	cmd = "python3 %s/hifiasm_single_task.py -i %s \
	-w %s \
	-o %s \
	-hap %s"%(code_dir,fastq_path,working_dir,final_contig_dir,hap_name)
	Popen(cmd, shell= True).wait()

	# write log file
	log_file_path = final_contig_dir+'/log.txt'
	cmd = "echo '%s' >> %s"%(hap_name,log_file_path)
	Popen(cmd, shell= True).wait()
	return 

def get_lines(file_path):
	with open(file_path,'r') as f:
		num_lines = len(f.readlines())
	return num_lines

## assembly hp1
pool = mp.Pool(num_of_threads)
count = 1
total_num = df.shape[0]
print('total_num',total_num)
for i in range(total_num):
	count+=1
	pb_name = 'PS%d_%d_%d'%(df['name'][i],df['start'][i],df['end'][i])
	hap_name = pb_name+'_hp1'
	fastq_path = input_dir+'/'+hap_name+'.fastq'
	pool.apply_async(assembly_one_hp, args=(fastq_path,working_dir,final_contig_dir,hap_name,code_dir))
	if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
		pool.close()
		# finished_tasks_num = get_lines(final_contig_dir+'/log.txt')
		while len(active_children()) > 1:
			# print(active_children())
			time.sleep(0.5)
		# while (finished_tasks_num == i+1 ) & (len(active_children()) > 1):
		# 	print(active_children())
		# 	time.sleep(0.5)
		# 	finished_tasks_num = get_lines(final_contig_dir+'/log.txt')
		# 	print("sleep 0.5 s")
		pool.join()

		if (count - 1) == total_num:
			print("finished all phase block assembly" )
		else:
			pool = Pool(num_of_threads)

## assembly hp2
pool = mp.Pool(num_of_threads)
count = 1
for i in range(total_num):
	# print("hp2",count)
	count+=1
	pb_name = 'PS%d_%d_%d'%(df['name'][i],df['start'][i],df['end'][i])
	hap_name = pb_name+'_hp2'
	fastq_path = input_dir+'/'+hap_name+'.fastq'
	pool.apply_async(assembly_one_hp, args=(fastq_path,working_dir,final_contig_dir,hap_name,code_dir))
	if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
		pool.close()
		
		# finished_tasks_num = get_lines(final_contig_dir+'/log.txt')

		# while finished_tasks_num < i+1 :
		# 	time.sleep(0.5)
		# 	finished_tasks_num = get_lines(final_contig_dir+'/log.txt')

		# 	# pass

		while len(active_children()) > 1:
			time.sleep(0.5)
		pool.join()

		if (count - 1) == total_num:
			print("finished all phase block assembly" )
		else:
			pool = Pool(num_of_threads)
