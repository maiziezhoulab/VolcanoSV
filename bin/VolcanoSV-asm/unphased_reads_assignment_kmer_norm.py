from subprocess import Popen
from argparse import ArgumentParser
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--phased_bam','-pbam')
parser.add_argument('--output_dir','-o')
parser.add_argument('--kmer_size','-k',type = int)
parser.add_argument('--batch_size','-b',type = int, default = 10000)
parser.add_argument('--threads_hash','-th', type = int, default = 100)
parser.add_argument('--threads_kc','-tkc', type = int, default = 50)
parser.add_argument('--threads_asg','-tas', type = int, default = 40)
parser.add_argument('--significance_level','-sigl',type = float, default = 0.1)
args = parser.parse_args()
phased_bam = args.phased_bam
output_dir = args.output_dir
k = args.kmer_size
batch_size = args.batch_size
threads_hash = args.threads_hash
threads_kc = args.threads_kc
threads_asg = args.threads_asg
significance_level = args.significance_level

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

import os
global code_dir
code_dir=os.path.dirname(os.path.realpath(__file__))+'/'

## make output dir
cmd = "mkdir -p "+output_dir
Popen(cmd, shell = True).wait()

# #######################################################
#
#
#                  hash kmer
#
#
# ######################################################


os.system("mkdir -p "+output_dir)
## bam to seq
logger.info("bam to seq...")
cmd = "python3 "+code_dir+"/bamtoseq_with_kmer_filter.py \
  -i %s  \
  -o %s/reads.seq \
  -k %d"%(
  	phased_bam,
  	output_dir,
  	k)
Popen(cmd, shell = True).wait()

## hash seq
logger.info("hash seq...")
cmd = "python3 "+code_dir+"/HashSeq.py \
  -i %s/reads.seq \
  -o %s/ \
  -k %d --n_thread %d --batch_size %s"%(
  	output_dir,
  	output_dir,
  	k,
  	threads_hash,
  	batch_size)
Popen(cmd, shell = True).wait()

# #######################################################
#
#
#                  Collect Phasing Information
#
#
# ######################################################

## create dir
logger.info("collect phasing info...")
phasing_info_dir = output_dir+'/phasing_info/'
read_hp_og_path = phasing_info_dir+"/read_hp_og.p"
reads_pb_path = phasing_info_dir+'/reads_pb.p'
os.system("mkdir -p "+phasing_info_dir)
cmd = "python3 "+code_dir+"/prepare_info_v1.py \
-i %s \
-o %s/pb_info.p \
--unphased_reads_assign_dict %s/reads_pb.p \
--read_hp_og_dict %s/read_hp_og.p -t 50"%(phased_bam, phasing_info_dir, phasing_info_dir,  phasing_info_dir)

Popen(cmd,shell=True).wait()



# #######################################################
#
#
#                  construct kmer database
#
#
# ######################################################

## construct kmer database
hash_path = output_dir+'/reads.hash'
name_path = output_dir+'/reads.name'
kmer_db_dir = output_dir+'/kmer_db_by_hap/'
logger.info("construct kmer database...")
cmd = "python3 "+code_dir+"/count_kmer_v1.py \
--read_hp_og_path %s --hash_path %s --name_path %s --output_dir %s -t %d"%(read_hp_og_path,
hash_path,
name_path,
 kmer_db_dir, 
 threads_kc)
Popen(cmd, shell = True).wait()



# #######################################################
#
#
#                  unphased reads assignment
#
#
# ######################################################


## split hash by unphased block
logger.info("split hash by up block...")
hash_by_up_dir = output_dir+'/up_est_hash/'
cmd = "python3 "+code_dir+"/split_hash_by_hp.py \
-i %s  -name %s  \
-dc %s -o %s"%(hash_path,
	name_path,
	reads_pb_path,
	hash_by_up_dir)
Popen(cmd, shell = True).wait()


## assign unphased reads
logger.info("assign unphased reads...")
cmd = "python3 "+code_dir+"/get_raw_kmer_overlap_count_unphased_est_pbs_v1.py \
--unphased_dir %s \
--kmer_db %s \
 -o %s \
  --prefix %s \
  -t %d -sigl %.4f"%(
  	hash_by_up_dir,
  	kmer_db_dir,
  	output_dir,
  	"unphased_reads_assignment_dippav",
  	threads_asg,
  	significance_level  	)
Popen(cmd, shell = True).wait()

# remove up_est_hash
logger.info("delete unphased hash files...")
cmd = "rm -r "+hash_by_up_dir
Popen(cmd, shell = True).wait()



# remove hash and seq files
logger.info("delete hash files...")
cmd = "rm  "+hash_path + " "+ output_dir+"/reads.seq"
Popen(cmd, shell = True).wait()

# remove kmer database
logger.info("delete kmer db...")
cmd = "rm -r "+ kmer_db_dir
Popen(cmd, shell = True).wait()










