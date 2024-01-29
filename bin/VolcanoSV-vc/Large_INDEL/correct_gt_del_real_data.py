import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--eval_dir','-eval')
parser.add_argument('--output_path','-o')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--sigfile','-sig')
parser.add_argument('--n_thread','-t', type = int, default = 50 )
parser.add_argument('--dtype','-d',choices=['ONT','CCS','CLR'])
parser.add_argument('--vtype','-v',choices=['INS','DEL'])
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
eval_dir = args.eval_dir
bamfile = args.bamfile
output_path = args.output_path
n_thread = args.n_thread
sigfile = args.sigfile
dtype = args.dtype
vtype = args.vtype

# bamfile = '/data/maiziezhou_lab/Datasets/Pacbio_data/NA24385/PacBio_CCS_15kb_20kb_chemistry2_reads_NA24385/NGMLR_aligned_bam/NA24385_aligned_by_ngmlr.sorted.bam'
# sigfile = '/data/maiziezhou_lab/CanLuo/long_reads_project/Variant_Caller_Result/cuteSV/hifi/Hifi_L1/out/DEL_gte30_cluster.txt'

from subprocess import Popen
import subprocess
from joblib import Parallel, delayed
from tqdm import tqdm
import os
import gzip
import pandas as pd
import numpy as np
import pysam
from collections import defaultdict

def load_vcf(vcffile):
    vars = []
    if '.gz' in vcffile:
        f = gzip.open(vcffile,'rt') 
    else:
        f = open(vcffile,'r')
    # with gzip.open(vcffile,'rt') as f:
    for line in f:
        if line[0]!='#':
            if 'SVTYPE=DEL' in line:
                gt = line.split()[-1].split(':')[0]
                svlen = abs(int(line.split('SVLEN=')[1].split(';')[0]))
                vars.append((gt,svlen,line))
    f.close()
    return vars

def load_sum(sum_file):
    with open(sum_file,'r') as f:
        s = f.read()
    dc = eval(s)
    print(dc)
    return dc  




def load_sig(sigfile):
    dc = defaultdict(list)
    sig_list = []

    with open(sigfile,'r') as f:
        for line in f:
            data = line.split()
            # print(data[0])
            _, chrom,pos,svlen,rname = data
            # rnames = [var.split(',')[-1] for var in data]
            # num_reads = len(data)
            pos = int(pos)
            svlen = int(svlen)
            final_sig = (chrom, pos,svlen)
            dc[final_sig].append(rname)
    
    for var, rnames in dc.items():
        chrom, pos,svlen = var
        new_var = (chrom, pos,svlen, len(rnames))
        sig_list.append(new_var)
    print("num of sig: ", len(sig_list))
    return sig_list





def match_varlist_siglist(sig_list, var_list, min_size_sim, max_shift_ratio):

    last_match_idx = 0
    num_sup_list =[]
    for var in tqdm(var_list, desc = 'extract sig support'):
        gt, svlen_var, line = var
        
        data = line.split()
        chrom_var = data[0]
        pos_var = int(data[1])
        # print(type(svlen_var), type(min_size_sim))
        max_size = svlen_var / min_size_sim
        min_size = svlen_var * min_size_sim
        max_shift = max(svlen_var * max_shift_ratio, 500)

        min_pos = pos_var - max_shift
        max_pos = pos_var + max_shift

        match_id_list = []
        num_support = 0
        for i in range(last_match_idx, len(sig_list)):
            chrom_sig, pos_sig, svlen_sig, num_reads = sig_list[i]
            if chrom_sig == chrom_var:
                if min_pos <= pos_sig <= max_pos:
                    match_id_list.append(i)
                    if min_size <= svlen_sig <= max_size:
                        num_support+= num_reads
                elif pos_sig > max_pos:
                    break

        for i in range(last_match_idx, -1, -1):
            chrom_sig, pos_sig, svlen_sig, num_reads = sig_list[i]
            if chrom_sig == chrom_var:
                if min_pos <= pos_sig <= max_pos:
                    match_id_list.append(i)
                    if min_size <= svlen_sig <= max_size:
                        num_support+= num_reads
                elif pos_sig < min_pos:
                    break

        if match_id_list:
            last_match_idx = min(match_id_list)
        
        num_sup_list.append(num_support)

    return num_sup_list

        
def count_reads_span_region(chrom, max_start, min_end, bamfile):
    samfile = pysam.AlignmentFile(bamfile)
    rnames = []
    for read in samfile.fetch(chrom,max_start, min_end, ):
        if ( read.reference_start < max_start) & ( read.reference_end > min_end):
            rnames.append(read.qname)
    return len(rnames)
    

def check_full_cover_reads(bamfile,chrom,pos,svlen):
    if svlen <= 1000:
        max_start = pos 
        min_end = pos + abs(svlen )

        depth = count_reads_span_region(chrom, max_start, min_end, bamfile)


    else:
        flanking = 150
        span = 100
        left_start = pos - flanking
        left_end = left_start + span 
        left_depth = count_reads_span_region(chrom, left_start, left_end, bamfile)

        right_start = pos + svlen + flanking
        right_end = right_start + span 
        right_depth = count_reads_span_region(chrom, right_start, right_end, bamfile)

        depth = (left_depth + right_depth )/2 

    return depth


def extract_sig_support(sigfile, var_list_call,   bamfile, outfile,min_size_sim, max_shift_ratio,dc_gt = None):

    # df = pd.read_csv(bndfile,sep = '\t')
    # chrom_list =  df['prefix'].apply(lambda x: x.split('_')[1])
    # pos_list = df['prefix'].apply(lambda x: int(x.split('_')[2]))
    # svlen_list = df['prefix'].apply(lambda x: int(x.split('_')[3][3:]))
    chrom_list =[ ]
    pos_list =[]
    svlen_list =[]
    gt_list_call = []
    svid_list = []
    for var in var_list_call:
        gt, svlen, line = var
        data = line.split()
        chrom = data[0]
        pos = int(data[1])
        svid = data[2]
        svid_list.append(svid)
        chrom_list.append(chrom)
        pos_list.append(pos)
        svlen_list.append(abs(svlen))
        gt_list_call.append(gt)

    

    sig_list = load_sig(sigfile)

    num_sup_list = match_varlist_siglist(sig_list, var_list_call, min_size_sim, max_shift_ratio)
    depth_list = Parallel(n_jobs=n_thread)(delayed(check_full_cover_reads)\
                                           (bamfile,chrom_list[i],pos_list[i],svlen_list[i])
                                             for i in tqdm(range(len(chrom_list)), 
                                                           desc = "check around breakpoint read depth"))
    
    num_sup_list = np.array(num_sup_list)
    depth_list = np.array(depth_list)
    n_ratio_list = num_sup_list / depth_list

    dc = {'svlen': svlen_list,
          'svid':svid_list,
          'call_gt': gt_list_call,
          'n_support': num_sup_list,
          'n_cov':depth_list,
          'n_ratio': n_ratio_list}
    
    if dc_gt is not None:
        svid_list_comp = [ var[2].split()[2] for var in vars_comp]
        gt_list_base = [ dc_gt[svid] if svid in dc_gt else '0/0' for svid in svid_list_comp]
        dc['base_gt']= gt_list_base
    df = pd.DataFrame(dc)
    df.to_csv(outfile, sep = '\t', index = False)
    

    return df



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


    
# def correct_gt_eval(df,t_large_11,t_small_11,t_large_01,t_small_01):
#     # new_gt = np.array(['01']*df.shape[0])

#     new_gt = df['call_gt'].values.copy()

#     cond_large_11 = ((df['svlen']>1000) & (df['call_gt']=='11'))
#     cond_small_11 = ((df['svlen']<=1000) & (df['call_gt']=='11'))
#     cond_large_01 = ((df['svlen']>1000) & (df['call_gt']=='01'))
#     cond_small_01 = ((df['svlen']<=1000) & (df['call_gt']=='01'))

#     if not np.isnan(t_large_11):

#         new_gt [ (cond_large_11 & (df['n_ratio']>t_large_11)) ]  = '11'
#         new_gt [ (cond_large_11 & (df['n_ratio']<=t_large_11)) ]  = '01'
#     if not np.isnan(t_small_11):
#         new_gt [ cond_small_11 & (df['n_ratio']>t_small_11) ]  = '11'
#         new_gt [ cond_small_11 & (df['n_ratio']<=t_small_11) ]  = '01'

#     if not np.isnan(t_large_01):
#         new_gt [ cond_large_01 & (df['n_ratio']>t_large_01) ]  = '11'
#         new_gt [ cond_large_01 & (df['n_ratio']<=t_large_01) ]  = '01'
#     if not np.isnan(t_small_01):
#         # print(t_small_01)
#         new_gt [ cond_small_01 & (df['n_ratio']>t_small_01) ]  = '11'
#         new_gt [ cond_small_01 & (df['n_ratio']<=t_small_01) ]  = '01'
    
#     tp_gt_old = (df['base_gt'] == df['call_gt']).sum()
#     tp_gt = (df['base_gt'] == new_gt).sum()
#     gt_cor = tp_gt/df.shape[0]
#     delta = tp_gt - tp_gt_old
#     print("TP_GT (old):",tp_gt_old,
#           "\nTP_GT (new):", tp_gt,
#           "\nGT concordance (old):",tp_gt_old/df.shape[0], 
#           "\nGT concordance (new):",gt_cor, 
#           "\nsuccessfully corrected GT:",delta, 
#           "\nnew GT_FP:",df.shape[0] - tp_gt)
#     return new_gt


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

  


            


if eval_dir is not None:
    tp_base_file = eval_dir+'/tp-base.vcf.gz'
    tp_comp_file = eval_dir+'/tp-comp.vcf.gz'
    sum_file = eval_dir+'/summary.json'
    dc = load_sum(sum_file)
    vars_base = load_vcf(tp_base_file)
    vars_comp = load_vcf(tp_comp_file)
    gt_list_base = [ var[0] for var in vars_base]
    svid_list_comp = [ var[2].split()[2] for var in vars_comp]
    dc_gt = dict(zip(svid_list_comp,gt_list_base))
else:
    dc_gt = None


vars_comp = load_vcf(input_path)
min_size_sim = 0.6
max_shift_ratio = 2.3
df = extract_sig_support(sigfile, vars_comp,  bamfile, output_path, min_size_sim, max_shift_ratio, dc_gt)


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










# if eval_dir is not None:
#     para_dir = "/data/maiziezhou_lab/CanLuo/long_reads_project/VCF_Collections/vcf_2023_02_03/vcfs/bin"
#     para_path = para_dir + f"/GT_correction_para_{dtype}_{vtype}.txt"
#     dc_para = read_para(para_path)

#     df = pd.read_csv(output_path,sep ='\t')
#     df = df[df['base_gt']!='0/0'].reset_index(drop = True)
#     df['call_gt'] = df['call_gt'].apply(lambda x: x.replace('/',''))
#     df['base_gt'] = df['base_gt'].apply(lambda x: x.replace('/',''))
#     new_n_cov= df['n_cov'].values.copy()
#     # mean_cov = 56
#     mean_cov = df['n_cov'].mean()
#     new_n_cov[df['svlen']>1000] = mean_cov
#     df['n_cov'] = new_n_cov
#     # df['n_ratio'] = df['n_support']/df['n_cov']
#     # print(df.head())
#     new_gt = correct_gt_eval(df,dc_para['t_large_11'],dc_para['t_small_11'],dc_para['t_large_01'],dc_para['t_small_01'],)





    




























