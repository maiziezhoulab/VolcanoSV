from joblib import Parallel, delayed
from tqdm import tqdm
from argparse import ArgumentParser
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--kmer_size','-k',type = int)
parser.add_argument('--n_thread','-t',type = int, default = 50)
parser.add_argument('--n_thread_w','-tw',type = int, default = 10)
parser.add_argument('--batch_size','-b', type = int, default = 10000)
args = parser.parse_args()
in_file = args.input_path
output_dir = args.output_dir
kmer_size = args.kmer_size
n_thread = args.n_thread 
batch_size = args.batch_size
n_thread_w = args.n_thread_w

temp_dir = output_dir+'/temp/'
out_file = output_dir+'/reads.hash'
import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

ONE_HOT_MAP = {'A':'00', 'C':'01', 'G':'10', 'T':'11'}

def one_hot_hash(kmer):
    kmer = "".join([ONE_HOT_MAP[u] for u in kmer])
    return str(int(kmer, 2))

def generate_kmer(seq,k):
	kmers = [  seq[i:i+k] for i in range(len(seq)-k)]
	return kmers

def gen_for_line_onehot(line, k):
    seq = line.split("\t")[2].strip()
    return " ".join([one_hot_hash(u) for u in generate_kmer(seq, k)])+'\n'

def fast_write(lines, fout, temp_dir,n_thread_w):
	os.system("mkdir -p "+temp_dir)
	def write_block(block,id,temp_dir):
		with open(temp_dir+'/temp%d.hash'%id,'w') as ftemp:
			ftemp.writelines(lines)
		return 

	n_per_block = len(lines)//n_thread_w+1

	## map
	Parallel(n_jobs=n_thread_w)(delayed(write_block)(lines[i*n_per_block: (i+1)*n_per_block],
	 i, 
	 temp_dir) for i in tqdm(range(n_thread_w), desc = "map"))

	## reduce

	for i in tqdm(range(n_thread_w),desc = "reduce"):
		with open(temp_dir+'/temp%d.hash'%i,'r') as ftemp:
			x = ftemp.readlines()
		fout.writelines(x)

	## delete temp
	os.sytem("rm %s/* "%temp_dir)
	return






def convert(in_file, out_file,  kmer_size, n_thread,  batch_size ):
    logger.info("start converting...")

    def f(lines,fout):

        sequences = Parallel(n_jobs=n_thread)(delayed(gen_for_line_onehot)(line, kmer_size) for line in tqdm(lines))
        fout.writelines(sequences)

        # hash_lines = [" ".join([str(v) for v in u])+'\n' for u in sequences]

        # fast_write(hash_lines, fout, temp_dir,n_thread_w)
        # fout.writelines(hash_lines)

        # for u in sequences:
        #     fout.write(" ".join([str(v) for v in u]))
        #     fout.write("\n")
        return 



    logger.info ("parameters: " + str(locals()))


    logger.info("creating hash ...")        

          
    cnt = 0
    with open(out_file, 'wt') as fout:
        with open(in_file) as fp:
            lines = []
            for line in fp:
                lines.append(line)
                if len(lines) >= batch_size:
                    f(lines, fout)
                    cnt += len(lines)
                    lines = []
                    logger.info ("written {} lines".format(cnt))

        if len(lines) > 0:
            f(lines, fout)
            cnt += len(lines)
            lines = []
            logger.info ("written {} lines".format(cnt))            
    
    logger.info("finish converting...")
    return 



convert(in_file, 
	out_file, 
	kmer_size,
	n_thread=n_thread,
	batch_size=batch_size)





