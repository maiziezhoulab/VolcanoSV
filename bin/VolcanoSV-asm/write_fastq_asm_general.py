import pickle
from subprocess import Popen
import pysam
import logging
import os
from tqdm import tqdm

def check_overlap(bedfile, pb_csv, chrom, out_dir):
	# print(out_dir)
	if not os.path.exists(out_dir):
		os.system("mkdir -p " + out_dir)

	code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

	bedtools = code_dir + "bedtools"

	# convert csv to bed file

	pb_bed = out_dir + "/pb_info.bed"
	cnt = 0
	with open(pb_csv,'r') as fin:
		with open(pb_bed, 'w') as fw:
			for line in fin:
				if cnt >0 :
					data = line[:-1].split(',')
					pb_name = 'PS'+line[:-1].replace(',','_')
					data[0] = chrom 
					data.append(pb_name)
					fw.write('\t'.join(data)+'\n')
				cnt +=1
	outfile = out_dir + "/pb_overlap.bed"
	cmd = f"{bedtools} intersect -a {pb_bed} -b {bedfile} -wo > {outfile}"
	Popen(cmd, shell= True).wait()
	inbed_pbs = []
	with open(outfile,'r') as f:
		for line in f:
			inbed_pbs.append(line.split()[3])

	return set(inbed_pbs)




def write_fastqs_general(bam_path, work_dir, kmer_dir, bed_file, chrom):

	inbed_output_dir = work_dir + "/fastq_by_hap_inbed/"
	outbed_output_dir =  work_dir + "/fastq_by_hap_outbed/"
	nobed_output_dir = work_dir + "/fastq_by_hap"
	read_hp_og_dc = kmer_dir + "/phasing_info/read_hp_og.p"
	nearest_neighbor_dc = kmer_dir + "unphased_reads_assignment_dippav_norm_asg.p"
	pb_csv = kmer_dir + "/phasing_info/pb_info.csv"

	## set logger

	logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
		level=logging.INFO,
		datefmt='%Y-%m-%d %H:%M:%S')
	logger = logging.getLogger("write fastq")


	##------- check overlap

	if bed_file is not None:
		olp_dir = work_dir + "/overlap_bed/"
		inbed_pbs = check_overlap(bed_file, pb_csv, chrom, olp_dir )
		##clean output dir
		logger.info("Deleting old output files...")
		os.system("rm %s/*fastq "%inbed_output_dir)
		os.system("mkdir -p "+inbed_output_dir)

		os.system("rm %s/*fastq "%outbed_output_dir)
		os.system("mkdir -p "+outbed_output_dir)
	else:
		logger.info("Deleting old output files...")
		os.system("rm %s/*fastq "%nobed_output_dir)
		os.system("mkdir -p "+nobed_output_dir)



	dc_up = pickle.load(open(nearest_neighbor_dc,'rb'))
	dc_og = pickle.load(open(read_hp_og_dc,'rb'))



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
			pb = '_'.join(hp.split('_')[:-1])
			if bed_file is not None:
				if pb in inbed_pbs:
					output_dir = inbed_output_dir
				else:
					output_dir = outbed_output_dir
			else:
				output_dir = nobed_output_dir
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
			# hp = dc_up[read.qname]		pb = '_'.join(hp.split('_')[:-1])
			# if bed_file is not None:
			# 	if pb in inbed_pbs:
				pb = '_'.join(hp.split('_')[:-1])
				if bed_file is not None:
					if pb in inbed_pbs:
						output_dir = inbed_output_dir
					else:
						output_dir = outbed_output_dir
				else:
					output_dir = nobed_output_dir
				out_path = output_dir+'/%s.fastq'%hp 
				seq = read.seq
				qual = ''.join(['!']*len(seq))
				with open(out_path,'a+') as f:
					f.writelines(['@'+read.qname+'\n',
						seq+'\n',
						'+\n',
						qual+'\n'])


if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
	parser.add_argument('--bam_path','-bam')
	parser.add_argument('--bed_file','-bed')
	parser.add_argument('--chrom','-chrom')
	parser.add_argument('--work_dir','-wd')
	parser.add_argument('--kmer_dir','-kd')

	args = parser.parse_args()
	bam_path = args.bam_path
	bed_file = args.bed_file
	chrom = args.chrom
	work_dir = args.work_dir
	kmer_dir = args.kmer_dir
	write_fastqs_general(bam_path, work_dir, kmer_dir, bed_file, chrom)
