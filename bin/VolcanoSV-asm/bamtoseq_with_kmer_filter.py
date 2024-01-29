import pysam
import pickle
import os
import pandas as pd
from argparse import ArgumentParser
import logging
## set logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
parser.add_argument('--kmer','-k',type = int)
# parser.add_argument('--pb_info_path','-p')
args = parser.parse_args()
output_path = args.output_path
bam_path = args.input_path
kmer = args.kmer
# pb_info_path = args.pb_info_path


# load samfile
samfile = pysam.AlignmentFile(bam_path)
samiter = samfile.fetch()

# load pb information
# df = pd.read_csv(pb_info_path)
# pb_dict = {}
# for i in range(df.shape[0]):
# 	pb_name = df['name'][i]
# 	pb_start = df['start'][i]
# 	pb_end = df['end'][i]
# 	pb_str = 'PS%d_%d_%d'%(pb_name,pb_start,pb_end)
# 	pb_dict[pb_name]=pb_str

os.system("rm  %s"%output_path)
write_name = []
i = 0
j = 0
line_list = []
filter_list = []
cnt = 0
for read in samiter:
	read_name = read.query_name
	seq = read.seq
	if cnt%10000 ==0:
		logger.info("processed %d reads"%cnt)
	cnt +=1
	if (len(seq)>=(kmer)) and ('N' not in seq):
		i+=1
		line_list.append("%d\t"%i+read_name+'\t'+seq+'\n')
		write_name.append(read_name+'\n')
	else:
		j+=1
		filter_list.append("%d\t"%j+read_name+'\t'+seq+'\n')

with open(output_path,'w') as f:
	f.writelines(line_list)

with open(output_path.replace('.seq','.name'),'w') as f:
	f.writelines(write_name)

with open(output_path.replace('.seq','_filterout.seq'),'w') as f:
	f.writelines(filter_list)