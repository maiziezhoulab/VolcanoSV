
'''ml GCC/10.2.0 BCFtools/1.16'''
def filter_vcf_by_bed(input_vcf_path,  output_folder, prefix,bed_path = None):
    # create output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    # filter vcf by size
    if 'vcf.gz' in input_vcf_path:

        cmd = f"zcat {input_vcf_path} | awk '($1 ~ /^#/ ||( (length($5) - length($(4)) >=2 ) && (length($5) - length($(4)) <= 50)) || ( (length($4) - length($(5)) >=2 ) && (length($4) - length($(5)) <= 50)))' > {output_folder}/raw.vcf"
        # print(cmd)
        os.system(cmd)
    else:
        cmd = f"cat {input_vcf_path} | awk '($1 ~ /^#/ ||( (length($5) - length($(4)) >=2 ) && (length($5) - length($(4)) <= 50)) || ( (length($4) - length($(5)) >=2 ) && (length($4) - length($(5)) <= 50)))'>  {output_folder}/raw.vcf"
        # print(cmd)
        os.system(cmd)

    # sort and index input VCF
    sorted_vcf_path = os.path.join(output_folder, "sorted.vcf.gz")
    subprocess.run(f"bcftools sort  {output_folder}/raw.vcf -Oz -o {sorted_vcf_path}", shell=True, check=True)
    subprocess.run(f"bcftools index -f {sorted_vcf_path}", shell=True, check=True)

    if bed_path is not None:
        # filter VCF based on bed file
        filtered_vcf_path = os.path.join(output_folder, prefix +'.vcf' )
        subprocess.run(f"bcftools view -R {bed_path} {sorted_vcf_path} -Ov -o {filtered_vcf_path}", shell=True, check=True)

        os.system(f"rm {output_folder}/sorted.vcf.gz* {output_folder}/raw.vcf")
        return filtered_vcf_path
    else:
        return sorted_vcf_path


def extract_indel_context(input_path, output_path):
    with open(output_path,'w') as fw:
        fw.write("svid\tqname\tstart\tend\tread_seq_start_end\talt\n")
        with open(input_path,'r') as f:
            for line in f:
                if line[0]!='#':
                    if 'CONTEXT=' in line:
                        data = line.split()
                        kmer = data[7].split('CONTEXT=')[1]
                        svid = data[2]
                        line = f"{svid}\t.\t.\t.\t{kmer}\t.\n"
                        fw.write(line)
    return 

def align_fa(reference, contig, hp, outdir,n_thread):
    cmd = f'''
    minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 --secondary=no -a -t {n_thread} --eqx -Y -O 5,56 -E 4,1 -B 5 \
    {reference} \
    {contig} | \
    samtools sort -o {outdir}/{hp}.bam;
    samtools index {outdir}/{hp}.bam
    '''
    logger.info(cmd)
    Popen(cmd, shell = True).wait()
    return 

def call_var(hp1_bam ,
    hp2_bam ,
    hp1_fa,
    hp2_fa,
    read_bam ,
    region ,
    output_dir ,
    n_thread ,
    kmer_size,
    ratio ,
    min_support ,
    restart ):

    ## set logger
    global logger, reference, bedfile
    logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(" ")

    stars0 = '\n'+'*'*40+'\n\n    '
    stars1 = '\n\n' + '*'*40
    os.makedirs(output_dir, exist_ok = True)

    if hp1_bam is None:
        # align contig
        logger.info(stars0 + "align contigs to reference" + stars1)
        align_fa(reference, hp1_fa, 'hp1', output_dir,n_thread)
        align_fa(reference, hp2_fa, 'hp2', output_dir,n_thread)
        hp1_bam = f'{output_dir}/hp1.bam'
        hp2_bam = f'{output_dir}/hp2.bam'

    # detect snp indel
    logger.info(stars0 + "detect SNP indel from both haplotypes" + stars1)
    cmd = f'''{code_dir}/htsbox/htsbox pileup -q5 -ecf \
        {reference} \
            {hp1_bam} {hp2_bam} -w 20 -r {region}  | {code_dir}/htsbox/htsbox bgzip > {output_dir}/var_raw.pair.vcf.gz'''
    logger.info(cmd)
    Popen(cmd, shell = True).wait()


    # pair
    logger.info(stars0 + "pair SNP, indel"+ stars1)
    cmd = f'''
    {code_dir}/k8 {code_dir}/dipcall//dipcall-aux.js vcfpair  -p /data/maiziezhou_lab/CanLuo/Software/dipcall/hs37d5.PAR.bed \
        {output_dir}/var_raw.pair.vcf.gz | {code_dir}/htsbox/htsbox bgzip > {output_dir}/var_raw.dip.vcf.gz
    '''
    logger.info(cmd)
    Popen(cmd, shell = True).wait()


    # reformat vcf; seperate multi-alt 
    logger.info(stars0 + "seperate multi-alt-allele"+ stars1)
    reformat_raw_vcf(f'{output_dir}/var_raw.dip.vcf.gz', f'{output_dir}/var_raw_rf.dip.vcf' )


    # filter vcf by size and bed
    logger.info(stars0 + f"filter 2-49 bp indel within {bedfile} region"+ stars1)
    filter_vcf_by_bed(f'{output_dir}/var_raw_rf.dip.vcf',  output_dir, 'indel_2_49', bedfile)


    # extract indel context from contig
    logger.info(stars0 + "extract indel context"+ stars1)
    extract_indel_context(f'{output_dir}/indel_2_49.vcf', f'{output_dir}/indel_2_49_context.txt'  )


    # check kmer support from reads bam file
    logger.info(stars0 + "check kmer support and remove false call"+ stars1)
    filter_indel( f'{output_dir}/indel_2_49.vcf', f'{output_dir}/indel_2_49_context.txt' , read_bam, 
                output_dir, 
                    None,
                    kmer_size ,
                    n_thread , ratio , min_support ,  restart )
    
    return 

if __name__ == '__main__':
    import argparse
    from argparse import ArgumentParser

    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--hp1_bam','-hp1', help = "either hp1 bam or hp1 fa need to be provided")
    parser.add_argument('--hp2_bam','-hp2', help = "either hp2 bam or hp2 fa need to be provided")
    parser.add_argument('--hp1_fa','-fhp1', help = "either hp1 bam or hp1 fa need to be provided")
    parser.add_argument('--hp2_fa','-fhp2', help = "either hp2 bam or hp2 fa need to be provided")
    parser.add_argument('--read_bam','-rbam')
    parser.add_argument('--output_dir','-o')
    parser.add_argument('--reference','-ref')

    parser.add_argument('--bedfile','-bed', help = "optional; a high confidence bed file")
    parser.add_argument('--region','-r', help = "optional; exmaple: chr21:2000000-2100000")
    parser.add_argument("--n_thread",'-t', type = int, default = 50)
    parser.add_argument('--kmer_size','-k', type = int, default = 15)
    parser.add_argument("--ratio",'-rt', type = float, default = 0.3, help = "maximum bad kmer ratio")
    parser.add_argument('--min_support','-ms', type = int, default = 5, help = "maximum support for bad kmer")
    parser.add_argument('--restart','-rs', action='store_true', help = "restart mode; assume there is kmer support file already.")
    parser.add_argument('--eval','-e', action='store_true')
    args = parser.parse_args()
    hp1_bam = args.hp1_bam
    hp2_bam = args.hp2_bam
    hp1_fa = args.hp1_fa
    hp2_fa = args.hp2_fa
    reference = args.reference
    bedfile = args.bedfile
    read_bam = args.read_bam
    region = args.region
    output_dir = args.output_dir
    n_thread = args.n_thread
    kmer_size = args.kmer_size
    ratio = args.ratio
    min_support = args.min_support
    restart = args.restart
    eval = args.eval

    import logging
    import subprocess
    from subprocess import Popen
    from reformat_dipcall import *
    from check_reads_kmer_support import * 
    import os
    code_dir = os.path.dirname(os.path.realpath(__file__))+'/'


    call_var(hp1_bam ,
    hp2_bam ,
    hp1_fa,
    hp2_fa,
    read_bam ,
    region ,
    output_dir ,
    n_thread ,
    kmer_size,
    ratio ,
    min_support ,
    restart )




    
