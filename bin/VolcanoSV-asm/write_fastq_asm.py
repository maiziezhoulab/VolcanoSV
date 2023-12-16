import pickle
from subprocess import Popen
from argparse import ArgumentParser
import pysam
import logging
import os
from tqdm import tqdm
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--bam_path','-bam')
parser.add_argument('--read_hp_og_dc','-dc_og')
parser.add_argument('--nearest_neighbor_dc','-dc_up')
parser.add_argument('--output_dir','-o')
parser.add_argument('--delete_intermediate_file','-d', action='store_true')
args = parser.parse_args()
bam_path = args.bam_path
read_hp_og_dc = args.read_hp_og_dc
nearest_neighbor_dc = args.nearest_neighbor_dc
output_dir = args.output_dir
## set logger

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("write fastq")


##clean output dir
logger.info("Deleting old output files...")
os.system("rm %s/*fastq "%output_dir)
os.system("mkdir -p "+output_dir)

dc_up = pickle.load(open(nearest_neighbor_dc,'rb'))
dc_og = pickle.load(open(read_hp_og_dc,'rb'))

# merge dc

# for read in dc_up:
# 	if read not in dc_og:
# 		dc_og[read] = dc_up[read]

logger.info("Reads overlap; %d"%len(set(dc_up.keys())&set(dc_og.keys())))


logger.info("Counting reads...")
samfile = pysam.AlignmentFile(bam_path)
samiter = samfile.fetch()
num_reads = 0
for read in samiter:
	num_reads+=1
logger.info("number of reads: %d"%num_reads)


samiter = samfile.fetch()
for i in tqdm(range(num_reads), desc = "writing progress",miniters=2000):
	read = next(samiter)
	if read.qname in dc_og:
		hp = dc_og[read.qname]
		out_path = output_dir+'/%s.fastq'%hp 
		seq = read.seq
		qual = ''.join(['!']*len(seq))
		with open(out_path,'a+') as f:
			f.writelines(['@'+read.qname+'\n',
				seq+'\n',
				'+\n',
				qual+'\n'])
	if read.qname in dc_up:
		for hp in dc_up[read.qname]:
		# hp = dc_up[read.qname]
			out_path = output_dir+'/%s.fastq'%hp 
			seq = read.seq
			qual = ''.join(['!']*len(seq))
			with open(out_path,'a+') as f:
				f.writelines(['@'+read.qname+'\n',
					seq+'\n',
					'+\n',
					qual+'\n'])


				
