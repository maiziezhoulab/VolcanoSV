from argparse import ArgumentParser
from subprocess import Popen
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--contig_path','-contig')
parser.add_argument('--reference_path','-ref')
parser.add_argument('--signature_dir','-sigd', help = "pre-extracted reads signatures; can only specify either signature_dir or bamfile")
parser.add_argument('--rbam_file','-rbam', help = "reads bam file for reads signature extraction; can only specify either signature_dir or bamfile")
parser.add_argument('--output_dir','-o')
parser.add_argument('--data_type','-dtype',help="Hifi;CLR;ONT")
parser.add_argument('--chr_num','-chr')
###optional
parser.add_argument('--header_file','-header',help='optional;if not set, will use the default header')
parser.add_argument('--n_thread','-t', type = int, default = 10,help = 'default = 10')
parser.add_argument('--mem_per_thread','-mempt', default = '1G',help = "Set maximum memory per thread; suffix K/M/G recognized; default = 1G")

args = parser.parse_args()
contig_path = args.contig_path
reference_path = args.reference_path
signature_dir = args.signature_dir
rbam_file = args.rbam_file
output_dir = args.output_dir
data_type = args.data_type
chr_num = args.chr_num
###optional
header_file = args.header_file
n_thread = args.n_thread
mem_per_thread = args.mem_per_thread

print("\n\nraw ", signature_dir)
import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

code_dir=os.path.dirname(os.path.realpath(__file__))+'/'
if not header_file:
    header_file = code_dir + 'header'


logger.info("align contig to reference...")
prefix = contig_path.split('/')[-1].split('.')[0]
bam_path = output_dir+'/'+prefix+'.sorted.bam'
cmd = "minimap2 -a -x asm5 --cs -r2k -t %d \
	   %s \
	   %s \
	    | samtools sort -@ %d -m %s > %s"%(n_thread,reference_path,contig_path,n_thread, mem_per_thread,bam_path )
logger.info(cmd)
Popen(cmd,shell=True).wait()

cmd = "samtools index "+bam_path
logger.info(cmd)
Popen(cmd,shell=True).wait()





logger.info("Raw variant call by chromosome...")
cmd = "python3 "+code_dir+"/extract_contig_signature_%s.py \
-bam %s  -contig %s -header %s -ref %s -o %s -chr %s"%(data_type,bam_path,
 contig_path,
 header_file,
 reference_path,
 output_dir,
 chr_num)

Popen(cmd,shell=True).wait()

### combine all vcf together

cmd = "cat %s/volcano_variant_chr*vcf | grep -v \"#\" > %s/body "%(output_dir,output_dir)
Popen(cmd,shell=True).wait()
cmd = "cat %s %s/body > %s/volcano_raw_variant.vcf; rm %s/body "%(header_file,output_dir,output_dir,output_dir)
Popen(cmd,shell=True).wait()


if signature_dir is None:
	cmd = f'''python3 {code_dir}/extract_reads_signature.py \
		-i {rbam_file} \
		-o {output_dir} -chr {chr_num}'''
	Popen(cmd,shell=True).wait()
	signature_dir = output_dir+"/reads_signature/"


logger.info("Filter out false positive...")
vcf_path = output_dir+"/volcano_raw_variant.vcf"
vcf_path_filtered = output_dir+"/volcano_variant_filtered.vcf"
cmd = "python3 "+code_dir+"/FP_filter_v1.py \
-i %s -sigd %s -o %s"%(vcf_path,signature_dir, vcf_path_filtered)
Popen(cmd,shell=True).wait()



logger.info("Remove redundancy...")
# vcf_path_noredun = output_dir+"/volcano_variant_chr%d_no_redundancy.vcf"%chr_num
final_vcf_dir = output_dir+'/final_vcf/'
cmd = "python3 "+code_dir+"/remove_redundancy.py -i %s -o %s "%(
	vcf_path_filtered, final_vcf_dir)
Popen(cmd,shell=True).wait()




































