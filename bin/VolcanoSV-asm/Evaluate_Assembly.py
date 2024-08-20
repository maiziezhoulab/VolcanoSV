from subprocess import Popen
import os


def merge_fasta(input_dir,outdir):

    fasta_list =  [input_dir+"/chr"+str(i+1)+f"/assembly/final_contigs/{prefix}_final_contigs.fa" for i in range(22)]

    if not os.path.exists(outdir):
        os.system("mkdir -p " + outdir)
    outfile = outdir+"/assemblies.fa"

    fw = open(outfile,'w')

    for fasta in fasta_list:
        with open(fasta,'r') as f:
            for line in f:
                fw.write(line)
    fw.close()
    return 


def preprocess_flagger(input_dir, outdir, fq_file, dtype,bam_prefix, n_thread,mem_per_t, ):
    merge_fasta(input_dir, outdir)
    cmd = f'''sh {code_dir}/preprocess_flagger.sh {fq_file} {outdir} {bam_prefix} {dtype} {n_thread} {mem_per_t}'''
    Popen(cmd, shell = True).wait()


def run_flagger(out_dir,  suffix, bam_prefix, sample_name, dtype, n_thread):
    bam_file = out_dir + "/" + bam_prefix + "_MD.bam"
    fasta_file = out_dir + "/assemblies.fa.gz"
    fai_file = out_dir + "/assemblies.fa.gz.fai"

    with open(code_dir+"/flagger_config_template.json",'r') as f:
        s = f.read()

    if dtype == 'hifi':
        max_rd = '0.02'
    else:
        max_rd = "0.09"

    s = s.replace("<fai_file>", fai_file)\
    .replace("<suffix>",suffix)\
    .replace("<fasta_file>",fasta_file)\
    .replace("<bam_file>", bam_file)\
    .replace("sample_name",sample_name)\
    .replace("<max_rd>", max_rd)\
    .replace("<t>",str(n_thread))

    config_file = out_dir+"/inputs.json"
    with open(config_file, 'w') as f:
        f.write(s)

    flaggerdir =  code_dir+"/Flagger/"
    cmd = f'''java -jar {flaggerdir}/cromwell-85.jar run {flaggerdir}/flagger-0.3.3/wdls/workflows/flagger_end_to_end_no_variant_calling_no_ref_no_secphase.wdl \
        -i {config_file} -m {out_dir}/outputs.json
        '''
    Popen(cmd, shell = True).wait()

    out_file = f"{out_dir}/outputs.json"
    with open(out_file,'r') as f:
        for line in f:
            if "FlaggerEndToEndNoVariantCallingNoRefNoSecphase.finalBed" in line:
                break 

    bed_file = line[:-1].split(':')[1].strip().replace('"', '')

    cmd = '''grep Col %s|cut -f 1|awk -F'_' '{print $1"_"$2"_"$3"_"$4}' |sort|uniq > %s/%s_%s_collapse_hp_names.txt'''%(bed_file,sample_name, suffix, out_dir)
    Popen(cmd, shell = True).wait()


import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_dir','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--fastq_file','-fq')
parser.add_argument('--data_type','-dtype', chpices = ['hifi', 'pb','ont'])
parser.add_argument('--n_thread','-t', type = int, default = 22 )
parser.add_argument('--mem_per_thread','-m', default = '5G' )
parser.add_argument('--sample_name','-samp', default = 'Sample' )

parser.add_argument('--lib_name','-lib',  )

args = parser.parse_args()
input_dir = args.input_dir
out_dir = args.output_dir
fq_file = args.fastq_file
dtype = args.data_type
n_thread = args.n_thread
mem_per_t = args.m_per_thread
sample_name = args.sample_name
suffix = args.lib_name




code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

bam_prefix = "reads2contig"

preprocess_flagger(input_dir, out_dir, fq_file, dtype, bam_prefix , n_thread, mem_per_t)
run_flagger(out_dir,  suffix, bam_prefix, sample_name, dtype, n_thread)