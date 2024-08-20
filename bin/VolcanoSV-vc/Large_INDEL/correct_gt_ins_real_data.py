import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--eval_dir','-eval')
parser.add_argument('--output_path','-o')
parser.add_argument('--sig_file','-sig')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
parser.add_argument('--dtype','-d',choices=['ONT','Hifi','CLR'])
parser.add_argument('--vtype','-v',choices=['INS','DEL'])
args = parser.parse_args()
input_path = args.input_path
eval_dir = args.eval_dir
output_path = args.output_path
sig_file = args.sig_file
bamfile = args.bamfile
n_thread = args.n_thread
dtype = args.dtype
vtype = args.vtype


# from collections import defaultdict
from tqdm import tqdm
import pandas as pd
import numpy as np
import gzip
from subprocess import Popen
import subprocess
from joblib import Parallel, delayed
import pysam

def load_vcf(vcffile):
    sv_list = []
    if '.gz' in vcffile:
        f = gzip.open(vcffile,'rt')
    else:
        f = open(vcffile,'r')

    for line in f:
        if line[0]!='#':
            if 'SVTYPE=INS' in line:
                data = line.split()
                chrom = data[0]
                chrom = int(chrom[3:])
                pos = int(data[1])
                gt = data[-1].split(':')[0]
                svid = data[2]
                svlen = int(line.split('SVLEN=')[1].split(';')[0])
                sv_list.append([chrom,pos,svlen,gt,svid])

    f.close()
    return sv_list



def load_sumfile(sum_file):
    with open(sum_file,'r') as f:
        s = eval(f.read())
    print(s)
    return 


def load_sig_file(sig_file):
    print(f"loading sigfile {sig_file}...")


    # sig_dc = defaultdict(list)
    sig_dc = dict()

    with open(sig_file,'r') as f:
        for line in f:
            data = line.split()
            # print(data)
            _,chrom,pos,svlen,rname = data[:5]
            pos = int(pos)
            svlen = int(svlen)
            if svlen >= 30:
                try:
                    chrom = int(chrom[3:])
                    var = (chrom,pos,svlen)
                    if var in sig_dc:
                        sig_dc[var]+=1
                    else:
                        sig_dc[var] = 1
                except:
                    continue

            # sig_dc[(chrom,pos)].append(svlen)
    print('num different signature: ', len(sig_dc))

    with open(sig_file+'.gte30auto','w') as f:
        for var in sig_dc:
            chrom,pos,svlen = var
            line = f'{chrom}\t{pos}\t{svlen}\t{sig_dc[var]}\n'
            f.write(line)

    return sig_dc





def extract_sig_support(sv_list_comp, sigfile,  max_dist_ratio, min_size_sim):
    sig_dc = load_sig_file(sigfile)
    sig_list = list(sig_dc.keys())

    last_pos_match = 0
    cnt_list =[ ]
    match_list = []
    for i in tqdm(range(len(sv_list_comp))):
        var = sv_list_comp[i]
        chrom,pos,svlen,gt,_ = var 
        max_shift = max(500,svlen * max_dist_ratio )
        min_pos = pos - max_shift
        max_pos = pos + max_shift
        min_size = svlen * min_size_sim
        max_size = svlen / min_size_sim
        cnt_support = 0
        match_index_list =[]
        # from left to right
        match_list.append(last_pos_match)
        for j in range(last_pos_match, len(sig_list)):
            sig = sig_list[j]
            sig_chrom,sig_pos,sig_svlen = sig
            if sig_chrom == chrom:
                if (min_pos<= sig_pos <= max_pos) :
                    match_index_list.append(j)
                    if (min_size<= sig_svlen <= max_size):
                        cnt_support += sig_dc[sig]

                elif sig_pos > max_pos:
                    break 

            # elif sig_chrom > chrom:
            #     break 

        # from right to left
        for j in range(last_pos_match, -1, -1):
            sig = sig_list[j]
            sig_chrom,sig_pos,sig_svlen = sig
            if sig_chrom == chrom:
                if (min_pos<= sig_pos <= max_pos):
                    match_index_list.append(j)
                    if  (min_size<= sig_svlen <= max_size):
                        cnt_support += sig_dc[sig]

                elif sig_pos < min_pos:
                    break 

            # elif sig_chrom < chrom:
            #     break 

        if match_index_list:
            # update next searching center
            last_pos_match = min(match_index_list)

        
        
        cnt_list.append(cnt_support)

        # info_list.append(var+[sv_list_base[i][-1]]+[cnt_support])
    return cnt_list,match_list





def bam_cov_one_region(bamfile, region ):
    # outbam = outdir + '/' + prefix + '.bam'
    cmd = f'''samtools depth -J -aa {bamfile} -r {region} '''+"|  awk '{sum+=$3} END { print sum/NR}' "
    out = Popen([cmd], shell = True,
           stdout=subprocess.PIPE, 
           stderr=subprocess.STDOUT)
    stdout,stderr = out.communicate()
    x = stdout.decode('utf8')[:-1]
    ### just in case the index is older than bam file
    if 'index' in x:
        cov = eval(x.split('\n')[1])
    else:
        cov = eval(x)
    return cov


def check_full_cover_reads(bamfile,chrom,pos,flanking):
    samfile = pysam.AlignmentFile(bamfile)
    max_start = pos - flanking
    min_end = pos + flanking
    rnames = []
    for read in samfile.fetch(chrom,max_start, min_end):
        if ( read.reference_start < max_start) & ( read.reference_end > min_end):
            rnames.append(read.qname)
    return len(rnames)

def extract_local_read_cov(bamfile,flanking,sv_list_comp, n_thread):
    region_list = []
    pos_list = []
    chrom_list = []
    for sv in sv_list_comp:
        pos = sv[1]
        pos_list.append(pos)
        chrom = 'chr'+str(sv[0])
        chrom_list.append(chrom)
        rg = f'{chrom}:{pos-flanking}-{pos+flanking}'
        region_list.append(rg)
    # cov_list = Parallel(n_jobs=n_thread)(delayed(bam_cov_one_region)\
    #     (bamfile, region_list[i]) for i in tqdm(range(len(region_list))))

    cov_list = Parallel(n_jobs=n_thread)(delayed(check_full_cover_reads )\
        (bamfile, chrom_list[i], pos_list[i], flanking) for i in tqdm(range(len(pos_list)),desc = "check reads coverage"))


    return cov_list

def write_new_df(sv_list_comp, sig_file, bamfile,gt_base_list, output_path,max_dist_ratio, min_size_sim,flanking, n_thread):
    cnt_sp_list,match_list = extract_sig_support(sv_list_comp, sig_file, max_dist_ratio, min_size_sim)
    cov_list = extract_local_read_cov(bamfile,flanking,sv_list_comp, n_thread)

    df = pd.DataFrame(sv_list_comp,columns = ['chrom','pos','svlen','call_gt','svid'])
    df['match_id'] = match_list
    if gt_base_list is not None:
        df['base_gt'] = gt_base_list
    df['n_support'] = cnt_sp_list
    df['n_cov'] = cov_list
    df['n_ratio'] = df['n_support']/ df['n_cov']
    df.to_csv(output_path, index = False, sep = '\t')

    return 


def read_para(para_file):
    dc = {}
    with open(para_file,'r') as f:
        for line in f:
            key,val = line[:-1].split(':')
            if val=='nan':
                val = np.nan
            else:
                val = eval(val)
            dc[key] = val 
    return dc 







def correct_gt_eval(df,t_large_11,t_small_11,t_large_01,t_small_01):
    # new_gt = np.array(['01']*df.shape[0])

    new_gt = df['call_gt'].values.copy()

    cond_large_11 = ((df['svlen']>1000) & (df['call_gt']=='1/1'))
    cond_small_11 = ((df['svlen']<=1000) & (df['call_gt']=='1/1'))
    cond_large_01 = ((df['svlen']>1000) & (df['call_gt']=='0/1'))
    cond_small_01 = ((df['svlen']<=1000) & (df['call_gt']=='0/1'))

    if not np.isnan(t_large_11):

        new_gt [ (cond_large_11 & (df['n_ratio']>t_large_11)) ]  = '1/1'
        new_gt [ (cond_large_11 & (df['n_ratio']<=t_large_11)) ]  = '0/1'
    if not np.isnan(t_small_11):
        new_gt [ cond_small_11 & (df['n_ratio']>t_small_11) ]  = '1/1'
        new_gt [ cond_small_11 & (df['n_ratio']<=t_small_11) ]  = '0/1'

    if not np.isnan(t_large_01):
        new_gt [ cond_large_01 & (df['n_ratio']>t_large_01) ]  = '1/1'
        new_gt [ cond_large_01 & (df['n_ratio']<=t_large_01) ]  = '0/1'
    if not np.isnan(t_small_01):
        # print(t_small_01)
        new_gt [ cond_small_01 & (df['n_ratio']>t_small_01) ]  = '1/1'
        new_gt [ cond_small_01 & (df['n_ratio']<=t_small_01) ]  = '0/1'
    
    
    return new_gt


def write_new_gt_vcf(vcffile,outfile,df):
    # print(df.head())
    svid_list = df['svid'].values
    new_gt_list = df['new_gt'].values
    dc = dict(zip(svid_list, new_gt_list))

    with open(vcffile,'r') as fin:
        with open(outfile,'w') as fout:
            for line in fin:
                if line[0]=='#':
                    # fout.write(line)
                    pass
                elif f'SVTYPE={vtype}' in line:
                    data = line.split()
                    svid = data[2]
                    new_gt = dc[svid]
                    data[-1] = new_gt
                    line = '\t'.join(data)+'\n'
                    fout.write(line)
    return  




sv_list_real = load_vcf(input_path)
if eval_dir is not None:
    tp_base_file = eval_dir+'/tp-base.vcf.gz'
    tp_comp_file = eval_dir+'/tp-comp.vcf.gz'
    sum_file = eval_dir+'/summary.json'
    dc = load_sumfile(sum_file)
    vars_base = load_vcf(tp_base_file)
    vars_comp = load_vcf(tp_comp_file)
    gt_list_base = [ var[3] for var in vars_base]
    svid_list_comp = [ var[4] for var in vars_comp]
    svid_list_real = [ var[4] for var in sv_list_real]
    dc_gt = dict(zip(svid_list_comp,gt_list_base))
    gt_list_base_real = [ dc_gt[svid] if svid in dc_gt else '0/0' for svid in svid_list_real]
else:
    gt_list_base_real = None




max_dist_ratio = 2.3
min_size_sim = 0.6
flanking = 100

write_new_df(sv_list_real, sig_file, bamfile,gt_list_base_real, output_path,max_dist_ratio, min_size_sim,flanking, n_thread)

### GT correction

# para_dir = "/data/maiziezhou_lab/CanLuo/long_reads_project/VCF_Collections/vcf_2023_02_03/vcfs/bin"
import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
para_dir = code_dir +"/para/"

para_path = para_dir + f"/GT_correction_para_{dtype}_{vtype}.txt"
dc_para = read_para(para_path)

df = pd.read_csv(output_path,sep ='\t')


# new_n_cov= df['n_cov'].values.copy()
# mean_cov = df['n_cov'].mean()
# new_n_cov[df['svlen']>1000] = mean_cov
# df['n_cov'] = new_n_cov
# df['n_ratio'] = df['n_support']/df['n_cov']


new_gt = correct_gt_eval(df,dc_para['t_large_11'],dc_para['t_small_11'],dc_para['t_large_01'],dc_para['t_small_01'],)
df['new_gt'] = new_gt
df.to_csv(output_path+'.newgt',sep = '\t', index = False)

write_new_gt_vcf(input_path,input_path+f'.newgt.{vtype}',df)


if eval_dir is not None:
    df = pd.read_csv(output_path+'.newgt',sep ='\t')
    df = df[df['base_gt']!='0/0'].reset_index(drop = True)

    # df['base_gt'] = df['base_gt'].apply(lambda x: x.replace('/',''))
    
    # print(df.head())
    
    tp_gt_old = (df['base_gt'] == df['call_gt']).sum()
    tp_gt = (df['base_gt'] == df['new_gt']).sum()
    gt_cor = tp_gt/df.shape[0]
    delta = tp_gt - tp_gt_old
    print("TP_GT (old):",tp_gt_old,
          "\nTP_GT (new):", tp_gt,
          "\nGT concordance (old):",tp_gt_old/df.shape[0], 
          "\nGT concordance (new):",gt_cor, 
          "\nsuccessfully corrected GT:",delta, 
          "\nnew GT_FP:",df.shape[0] - tp_gt)







