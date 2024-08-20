from argparse import ArgumentParser
import numpy as np
import pandas as pd
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--dtype','-d',help = 'hifi/ont/clr')
parser.add_argument('--asm','-a',help ='other/volcano')
parser.add_argument('--vtype','-v',help ='apply filter to which variants, INS/DEL/INSDEL, default = INSDEL',choices=['INS','DEL','INSDEL'], default = 'INSDEL')

args = parser.parse_args()
input_path = args.input_path
dtype = args.dtype
asm = args.asm
vtype = args.vtype

assert vtype in ['INS','DEL','INSDEL']

import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
df_para = pd.read_csv(code_dir+'/filter_para.csv')
dfx = df_para[(df_para['asm']==asm)&(df_para['dtype']==dtype)].reset_index(drop = True)
# print(df_para)
# print(dfx)
lb_ins_r = dfx['lb_ins'][0]
rb_ins_r = dfx['rb_ins'][0]
lb_del_r = dfx['lb_del'][0]
rb_del_r = dfx['rb_del'][0]


csvfile = input_path.replace('.vcf','_cutesv_sig_support_mins30_fl1000.csv')
df = pd.read_csv(csvfile)
df['re_cov'] = df['cov']/df['svlen']


df_ins = df[df['svtype']=='INS']

if df_ins.shape[0]:
	if vtype !='DEL':
		med_cov_ins = np.quantile(df_ins['re_cov'],0.5)
		lb_ins = med_cov_ins * lb_ins_r
		rb_ins = med_cov_ins * rb_ins_r
		df_ins_filtered = df_ins[(df_ins['re_cov']>= lb_ins)&(df_ins['re_cov']<= rb_ins)]
	else:
		df_ins_filtered = df_ins.copy()
	
	ins_ids = set(df_ins_filtered['svid'].values)
else:
	ins_ids = set()



df_del = df[df['svtype']=='DEL']
if df_del.shape[0]:
	if vtype !='INS':
		med_cov_del = np.quantile(df_del['re_cov'],0.5)
		lb_del = med_cov_del * lb_del_r
		rb_del = med_cov_del * rb_del_r
		df_del_filtered = df_del[(df_del['re_cov']>= lb_del)&(df_del['re_cov']<= rb_del)]
	else:
		df_del_filtered = df_del.copy()

	del_ids = set(df_del_filtered['svid'].values)
else:
	del_ids = set()



pass_ids = ins_ids|del_ids

def calc_f1(tp,fp,svtype):
	if svtype == 'INS':
		total_base = 5077+ 204
	else:
		total_base = 3901 + 215

	recall = tp/total_base
	precision = tp/(tp+fp)
	f1 = 2* recall * precision / (recall + precision)
	return f1,recall,precision,tp,fp




# tp_ins_og = (df_ins['call_type']=='TP').sum()
# fp_ins_og = (df_ins['call_type']=='FP').sum()
# tp_del_og = (df_del['call_type']=='TP').sum()
# fp_del_og = (df_del['call_type']=='FP').sum()

# tp_ins_fl = (df_ins_filtered['call_type']=='TP').sum()
# fp_ins_fl = (df_ins_filtered['call_type']=='FP').sum()
# tp_del_fl = (df_del_filtered['call_type']=='TP').sum()
# fp_del_fl = (df_del_filtered['call_type']=='FP').sum()

# print("INS")
# print('before: ',calc_f1(tp_ins_og,fp_ins_og,'INS'))
# print('after: ',calc_f1(tp_ins_fl,fp_ins_fl,'INS'))

# print('DEL')
# print('before: ',calc_f1(tp_del_og,fp_del_og,'DEL'))
# print('after: ',calc_f1(tp_del_fl,fp_del_fl,'DEL'))



outvcf = input_path.replace('.vcf',f'_filter_{vtype}.vcf')


with open(outvcf,'w') as fw:
	with open(input_path,'r') as f:
		for line in f:
			if line[0]=='#':
				fw.write(line)
			else:
				svid = line.split()[2]
				if svid in pass_ids:
					fw.write(line)






