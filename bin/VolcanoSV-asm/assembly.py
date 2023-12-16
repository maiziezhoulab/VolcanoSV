from argparse import ArgumentParser
from subprocess import Popen
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--bam_path','-bam')
parser.add_argument('--read_hp_og_dc','-dc_og')
parser.add_argument('--nearest_neighbor_dc','-dc_up')
parser.add_argument('--pb_info_csv','-pinf')
parser.add_argument('--output_dir','-o')
parser.add_argument('--assembly_threads','-t',default = 5, type=int, help = "optional argument. default= 5")
parser.add_argument('--single_task_threads','-ta',default = 5, type=int, help = "optional argument. default= 5")
# parser.add_argument('--assembler','-asm',help = "hifiasm;flye", default = "flye")
parser.add_argument('--data_type','-dtype',help="CCS;CLR;ONT", choices= ['CCS','CLR','ONT'])
parser.add_argument('--delete_intermediate_file','-d', action='store_true')

args = parser.parse_args()
bam_path = args.bam_path
read_hp_og_dc = args.read_hp_og_dc
nearest_neighbor_dc = args.nearest_neighbor_dc
pb_info_csv = args.pb_info_csv
output_dir = args.output_dir
assembly_threads = args.assembly_threads
single_task_threads = args.single_task_threads
# assembler = args.assembler
data_type = args.data_type

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

os.system("mkdir -p "+output_dir)
fastq_dir = output_dir+"/fastq_by_hap/"

import os
global code_dir
code_dir=os.path.dirname(os.path.realpath(__file__))+'/'


# write fastq
logger.info("write fastq...")
cmd = "python3 "+code_dir+"/write_fastq_asm.py \
-bam %s \
-dc_og %s \
-dc_up %s \
-o %s "%(
	bam_path,
	read_hp_og_dc,
	nearest_neighbor_dc,
	fastq_dir)
Popen(cmd,shell=True).wait()

if data_type == 'CCS':
	assembler = 'hifiasm'
else:
	assembler = 'flye'
# assembly
logger.info("start assembly...")
if assembler == 'hifiasm':
	cmd = "python3 "+code_dir+"/run_assembly.py  \
	-i "+fastq_dir+" \
	-w "+output_dir+"/assembly_files/ \
	-o "+output_dir+"/final_contigs/ \
	-p "+pb_info_csv+" \
	-t %d"%(assembly_threads)
	Popen(cmd,shell=True).wait()
elif assembler == 'flye':
	cmd = "python3 /data/maiziezhou_lab/CanLuo/long_reads_project/DipPAV_pipeline/side_bin/bin_v2/run_assembly_diff_assembler.py \
	-i "+fastq_dir+" \
	-w "+output_dir+"/assembly_files/ \
	-o "+output_dir+"/final_contigs/ \
	-t %d -ta %d -tool flye -dtype %s "%(assembly_threads,single_task_threads,data_type)
	Popen(cmd,shell=True).wait()
else:
	logger.info("only allow hifiasm or flye as assembler")















