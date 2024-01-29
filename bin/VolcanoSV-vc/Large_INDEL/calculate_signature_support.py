from subprocess import Popen
import subprocess
import os
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import rankdata
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcffile','-v')
parser.add_argument('--cutesv_dir','-ct')
parser.add_argument('--n_thread','-t',type = int,default = 22, help = "number of threads")
parser.add_argument('--flanking','-f',type = int, default = 1000, help = 'flanking region around breakpoint')
parser.add_argument('--min_size','-s', type = int, default = 30, help = "min signature size")
parser.add_argument('--chr_num','-chr',type = int, 
					choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],
					default= None,
					help = "chrmosome number;Optional; if not provided, will assume input_dir contain chr1-chr22 results")

args = parser.parse_args()
vcffile = args.vcffile
cutesv_dir = args.cutesv_dir
n_thread = args.n_thread
flanking = args.flanking
min_size = args.min_size
chr_num = args.chr_num

def load_vcf(vcffile,svtype):
	info_list = []
	dc = defaultdict(list)
	with open(vcffile,'r') as f:
		for line in f:
			if line[0]!='#':
				if 'SVTYPE=%s'%svtype in line:


					data = line.split()
					chrom = data[0]
					svlen = int(data[7].split('SVLEN=')[1].split(';')[0])
					if abs(svlen)>= min_size:
						start = int(data[1])
						svid = data[2]
						gt = data[-1].split(':')[0]
						if svtype == 'INS':
							end = start+1
						else:
							end = start-svlen
						info_list.append((chrom,start,end, svlen,svid,gt))
						dc[chrom].append((start,end,svlen,svid,gt,svtype))

	return dc

def load_sig(sig_path,svtype):
	info_dc =defaultdict(list)
	cnt = 0
	with open(sig_path,'r')as f:
		for line in f:
			cnt+=1
			# if cnt > 10000:
			# 	break
			data = line.split()
			# print(data)
			chrom = data[1]
			start = int(data[2])
			svlen = int(data[3])
		
			if svtype == 'INS':

				end = start+1
			else:
				svlen = -int(svlen)
				end = start - svlen
			
			info_dc[chrom].append((start,end,svlen))
	return info_dc

def calc_ins_call_cov(call_list,sig_list):
	pos_call_list = [call[0] for call in call_list]
	pos_sig_list = [call[0] for call in sig_list]
	len_sig_list = [call[2] for call in sig_list]
	uniq_pos = sorted(list(set(pos_sig_list)))
	dc_map = dict(zip(uniq_pos,list(range(len(uniq_pos)))))
	rank_list = np.vectorize(dc_map.get)(pos_sig_list)

	weights_list = np.bincount(rank_list, weights=len_sig_list)

	sorted_call_pos = sorted(list(set(pos_call_list)))
	bins_list = [ (pos-flanking, pos + flanking) for pos in sorted_call_pos]
	start_i = 0
	bin_w_list = np.zeros(len(bins_list))
	for j in range(len(uniq_pos)):
		pos = uniq_pos[j]
		w = weights_list[j]
		cnt = 0
		for i in range(start_i, len(bins_list)):
			lb,rb = bins_list[i]
			if lb > pos:
				break
			if (pos<=rb):
				cnt +=1 
				if cnt==1:
					real_i = i
				bin_w_list[i]+=w 
		if cnt:
			start_i = real_i

	# bin_w_list1 = np.zeros(len(bins_list))
	# for j in tqdm(range(len(uniq_pos))):
	# 	pos = uniq_pos[j]
	# 	w = weights_list[j]
	# 	for i in range(len(bins_list)):
	# 		lb,rb = bins_list[i]
	# 		if lb > pos:
	# 			break
	# 		if (pos<=rb):
	# 			bin_w_list1[i]+=w 
	# print(bin_w_list)
	# print(bin_w_list1)
	# print((bin_w_list==bin_w_list1).mean())
	dc = dict(zip(sorted_call_pos,bin_w_list))
	return dc

def sort_region(bed_list,sort_by):
	if sort_by=='end':
		x = 1
	else:
		x = 0
		key_list = [bed[x] for bed in bed_list]

	ids = np.argsort(key_list)
	sorted_bed = [bed_list[i] for i in ids]
	return sorted_bed

def calc_del_call_cov(call_list,sig_list):
	#### intialize data for sig
	dc_sig_start = defaultdict(list)
	dc_sig_end = defaultdict(list)
	siglen_list = [sig[2] for sig in sig_list]
	bed_sig_list = [(sig[0],sig[1]) for sig in sig_list]
	for i in range(len(sig_list)):
		start,end,svlen = sig_list[i]
		dc_sig_start[start].append(i)
		dc_sig_end[end].append(i)
	sig_start_list_sorted = sorted(list(dc_sig_start.keys()))
	sig_end_list_sorted = sorted(list(dc_sig_end.keys()))


	#### intialize data for bed
	bed_call_list = [(sig[0]-flanking,sig[1]+flanking) for sig in call_list]
	bed_call_list_sorted = sort_region(bed_call_list,'start')
	dc_call_start = defaultdict(list)
	dc_call_end = defaultdict(list)
	for i in range(len(bed_call_list)):
		start,end = bed_call_list[i]
		dc_call_start[start].append(i)
		dc_call_end[end].append(i)
	call_start_list_sorted = sorted(list(dc_call_start.keys()))
	call_end_list_sorted = sorted(list(dc_call_end.keys()))


	# create dictionary for final match
	dc_from_bnd = defaultdict(list)

	### find the situation when at least one end of the sig is in between of the call region

	start_i = 0
	for j in range(len(bed_call_list_sorted)):
		lb,rb = bed_call_list_sorted[j]
		cnt = 0
		for i in range(start_i,len(sig_start_list_sorted)):
			pos = sig_start_list_sorted[i]
			if pos > rb:
				break
			if pos>= lb:
				cnt+=1
				if cnt == 1:
					real_i = i 
				dc_from_bnd[j].extend(dc_sig_start[pos])
		if cnt:
			start_i = real_i


	start_i = 0
	for j in range(len(bed_call_list_sorted)):
		lb,rb = bed_call_list_sorted[j]
		cnt = 0
		for i in range(start_i,len(sig_end_list_sorted)):
			pos = sig_end_list_sorted[i]
			if pos > rb:
				break
			if pos>= lb:
				cnt+=1
				if cnt == 1:
					real_i = i 
				dc_from_bnd[j].extend(dc_sig_end[pos])
		if cnt:
			start_i = real_i

	### find the situation when the whole call region is included by sig region

	start_i = 0
	for j in range(len(bed_sig_list)):
		lb,rb = bed_sig_list[j]
		cnt = 0
		for i in range(start_i,len(call_start_list_sorted)):
			pos = call_start_list_sorted[i]
			if pos > rb:
				break
			if pos>= lb:
				cnt+=1
				if cnt == 1:
					real_i = i 
				call_ids = dc_call_start[pos]
				for call_id in call_ids:
					dc_from_bnd[call_id].append(j)
		if cnt:
			start_i = real_i

	start_i = 0
	for j in range(len(bed_sig_list)):
		lb,rb = bed_sig_list[j]
		cnt = 0
		for i in range(start_i,len(call_end_list_sorted)):
			pos = call_end_list_sorted[i]
			if pos > rb:
				break
			if pos>= lb:
				cnt+=1
				if cnt == 1:
					real_i = i 
				call_ids = dc_call_end[pos]
				for call_id in call_ids:
					dc_from_bnd[call_id].append(j)
		if cnt:
			start_i = real_i


	dc_svlen = {}
	for key in dc_from_bnd:
		## remove redundancy
		dc_from_bnd[key] = list(set(dc_from_bnd[key]))

		### now calculate the depth
		total_svlen = sum([siglen_list[i] for i in dc_from_bnd[key]])

		### match back 
		bed_call = bed_call_list_sorted[key]
		del_bed = (bed_call[0]+flanking,bed_call[1]-flanking)
		dc_svlen[del_bed] = total_svlen

	# ### very basic, slow but stable algorithm, can be used as validation
	# dc_from_bnd1 = defaultdict(list)
	# for i in tqdm(range(len(bed_sig_list))):
	# 	bed = bed_sig_list[i]
	# 	start, end = bed
	# 	for j in range(len(bed_call_list_sorted)):
	# 		lb,rb = bed_call_list_sorted[j]
	# 		if lb>end:
	# 			break
	# 		overlap = min(rb,end)-max(start,lb)
	# 		if overlap>=0:
	# 			dc_from_bnd1[j].append(i)

	# for key in dc_from_bnd1:
	# 	## remove redundancy
	# 	dc_from_bnd1[key] = list(set(dc_from_bnd1[key]))

	# match = 0
	# for key in dc_from_bnd:
	# 	if len(dc_from_bnd[key]) == len(dc_from_bnd1[key]):
	# 		match+=1
	# 	else:
	# 		print(key,len(dc_from_bnd[key]),len(dc_from_bnd1[key]),set(dc_from_bnd1[key])-set(dc_from_bnd[key]))

	# print(match/len(dc_from_bnd))
	return dc_svlen






dc_sig_ins = load_sig(cutesv_dir+'/INS.sigs','INS')
dc_sig_del = load_sig(cutesv_dir+'/DEL.sigs','DEL')
dc_call_ins = load_vcf(vcffile,'INS')
dc_call_del = load_vcf(vcffile,'DEL')


if chr_num is None:
	ins_results = Parallel(n_jobs=n_thread)(delayed(calc_ins_call_cov)(dc_call_ins['chr'+str(i)],dc_sig_ins['chr'+str(i)]) for i in tqdm(range(1,23)))
	del_results = Parallel(n_jobs=n_thread)(delayed(calc_del_call_cov)(dc_call_del['chr'+str(i)],dc_sig_del['chr'+str(i)]) for i in tqdm(range(1,23)))
else:
	ins_results = calc_ins_call_cov(dc_call_ins['chr'+str(chr_num)],dc_sig_ins['chr'+str(chr_num)]) 
	del_results = calc_del_call_cov(dc_call_del['chr'+str(chr_num)],dc_sig_del['chr'+str(chr_num)]) 

##### add cov back
	
if chr_num is None:
	final_info = []
	for i in range(22):
		chrom = 'chr'+str(i+1)
		if chrom in dc_call_ins:
			call_list = dc_call_ins[chrom]
			dc_ins_cov = ins_results[i]
			for call in call_list:
				start,end,svlen,svid,gt,svtype = call
				if start in dc_ins_cov:
					cov = dc_ins_cov[start]
				else:
					cov = 0
				final_info.append([start,end,svlen,svid,gt,svtype,cov])

	for i in range(22):
		chrom = 'chr'+str(i+1)
		if chrom in dc_call_del:
			call_list = dc_call_del[chrom]
			dc_del_cov = del_results[i]
			for call in call_list:
				start,end,svlen,svid,gt,svtype = call
				if (start,end) in dc_del_cov:
					cov = dc_del_cov[(start,end)]
				else:
					cov = 0
				final_info.append([start,end,svlen,svid,gt,svtype,cov])

else:
	final_info = []
	chrom = 'chr'+str(chr_num)
	if chrom in dc_call_ins:
		call_list = dc_call_ins[chrom]
		dc_ins_cov = ins_results
		for call in call_list:
			start,end,svlen,svid,gt,svtype = call
			if start in dc_ins_cov:
				cov = dc_ins_cov[start]
			else:
				cov = 0
			final_info.append([start,end,svlen,svid,gt,svtype,cov])


	if chrom in dc_call_del:
		call_list = dc_call_del[chrom]
		dc_del_cov = del_results
		for call in call_list:
			start,end,svlen,svid,gt,svtype = call
			if (start,end) in dc_del_cov:
				cov = dc_del_cov[(start,end)]
			else:
				cov = 0
			final_info.append([start,end,svlen,svid,gt,svtype,cov])






df = pd.DataFrame(final_info,columns = ['start','end','svlen','svid','gt','svtype','cov'])
df['rel_cov'] = df['cov']/df['svlen']
input_dir = os.path.dirname(vcffile)
basename = os.path.basename(vcffile).split('.')[0]
csv = input_dir+'/'+basename+'_cutesv_sig_support_mins%d_fl%d.csv'%(min_size,flanking)
df.to_csv(csv, index = False)












