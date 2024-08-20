import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcffile','-vcf')
parser.add_argument('--bamfile','-bam', help ="only needed when presig is not provided")
parser.add_argument('--reference','-ref', help ="only needed when presig is not provided")
parser.add_argument('--pre_cutesig','-presig', help = "pre-extracted cutesv signature directory;optional; if not provided, will generate a new one")
parser.add_argument('--dtype','-dtype', choices = ['Hifi','CLR','ONT'])
parser.add_argument('--chr_num','-chr',type = int, 
					choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],
					default= None,
					help = "chrmosome number;Optional; if not provided, will assume input_dir contain chr1-chr22 results")
parser.add_argument('--n_thread','-t', type = int, default = 22 )
args = parser.parse_args()


vcffile=args.vcffile
bamfile=args.bamfile
reference=args.reference
pre_cutesig= args.pre_cutesig
dtype= args.dtype
chr_num = args.chr_num
t = args.n_thread

ft_vtype='DEL'


import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")








import os 
from subprocess import Popen

filtered_vcf = os.path.join(os.path.dirname(vcffile), f"{os.path.basename(vcffile)[:-4]}_filter_{ft_vtype}.vcf")
dtype_lc = dtype.lower()
gtdir = os.path.join(os.path.dirname(vcffile), "GT_Correction")
final_vcf = os.path.join(os.path.dirname(vcffile), f"variants_filtered_GT_corrected.vcf")
workdir = os.path.dirname(vcffile)

script_dir = os.path.dirname(__file__)  # Assuming this script is in the same directory as your Python script



################################# filter out false positive #####################################


# cuteSv extract reads signature

if pre_cutesig is None:
    new_cutesig = os.path.join(os.path.dirname(vcffile), "cute_sig")
    if os.path.exists(new_cutesig):
        os.system("rm -r " + new_cutesig)

    os.system("mkdir -p "+ new_cutesig)
    

    if chr_num is None:
        bed_para = " "
    else:
        # get chromosome length
        fai_file = reference + ".fai"
        with open(fai_file,'r') as f:
            chr_len = int(f.readlines()[chr_num -1 ].split()[1])
        bed_file = new_cutesig + "/sample.bed"
        with open(bed_file, 'w') as f:
            f.write("chr"+str(chr_num)+"\t1\t"+str(chr_len)+'\n')
        
        bed_para = " -include_bed " + bed_file


    if  dtype == "Hifi":
        para ='''--max_cluster_bias_INS      1000 \
        --diff_ratio_merging_INS    0.9 \
        --max_cluster_bias_DEL  1000 \
        --diff_ratio_merging_DEL    0.5 '''
    elif dtype == 'CLR':
        para = '''--max_cluster_bias_INS      100 \
        --diff_ratio_merging_INS    0.3 \
        --max_cluster_bias_DEL  200 \
        --diff_ratio_merging_DEL    0.5 '''
    else:
        para = '''--max_cluster_bias_INS      100 \
        --diff_ratio_merging_INS    0.3 \
        --max_cluster_bias_DEL  100 \
        --diff_ratio_merging_DEL    0.3 '''
    


    #cmd = f'''python3 {script_dir}/cuteSV/cuteSV \
    #{bamfile} \
    #{reference} \
    #{new_cutesig}/test.vcf \
    #{new_cutesig} \
    #{para} \
    #{bed_para} --retain_work_dir -t {t}'''

    cmd = f'''python3 {script_dir}/sig_extract.py \
    {bamfile} \
    {reference} \
    {new_cutesig} \
    {para} \
    {bed_para} -t {t}'''
    print(cmd)
    Popen(cmd, shell = True).wait()

    sigdir = new_cutesig 
else:
    sigdir = pre_cutesig



if chr_num is None:
    chr_para = " "
else:
    chr_para = " -chr " + str(chr_num)

# calculate signature support for each variant based on cuteSV signature
logger.info("------------------Calculate signature support for each variant")

cmd = f'''python3 {script_dir}/calculate_signature_support.py  \
    -v {vcffile} \
    -ct {sigdir} {chr_para}'''
Popen(cmd, shell = True).wait()

# Use empirical threshold to filter SV
logger.info("------------------Use empirical threshold to filter SV")
cmd = f'''python3 {script_dir}/filter_vcf_by_sig_cov_insdel.py \
    -i {vcffile} \
    -d {dtype_lc} \
    -a volcano \
    -v {ft_vtype}'''
Popen(cmd, shell = True).wait()



############################################## GT correction ###################################
logger.info("------------------Correct DEL GT")
os.system(f'mkdir -p {gtdir}')

cmd = f'''python3 {script_dir}/correct_gt_del_real_data.py \
    -i {filtered_vcf}  \
    -o {gtdir}/bnd_del_real.tsv \
    -bam {bamfile} \
    -sig {sigdir}/DEL.sigs \
    -t {t} \
    -d {dtype} -v DEL'''
Popen(cmd, shell = True).wait()


logger.info("------------------Correct INS GT")
cmd = f'''python3 {script_dir}/correct_gt_ins_real_data.py \
    -i {filtered_vcf}  \
    -o {gtdir}/bnd_ins_real.tsv \
    -bam {bamfile} \
    -sig {sigdir}/INS.sigs \
    -t {t} \
    -d {dtype} -v INS'''
Popen(cmd, shell = True).wait()

cmd = f"grep '#' {filtered_vcf} > {workdir}/header"
Popen(cmd, shell = True).wait()

cmd = f"cat {workdir}/header {filtered_vcf}.newgt.DEL {filtered_vcf}.newgt.INS | vcf-sort > {final_vcf}"
Popen(cmd, shell = True).wait()

