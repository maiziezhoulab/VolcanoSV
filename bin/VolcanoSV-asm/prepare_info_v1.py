import pysam 
import pickle
import pandas as pd
from argparse import ArgumentParser
from tqdm import tqdm
import logging
from joblib import Parallel, delayed
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i') ## phased_bam_file
parser.add_argument('--pb_info_dict','-o')
parser.add_argument('--unphased_reads_assign_dict')
parser.add_argument('--read_hp_og_dict')
parser.add_argument('--n_thread','-t', default = 10, type = int, help = "optional, default=10")
args = parser.parse_args()

bam_path = args.input_path
pb_info_dict = args.pb_info_dict
unphased_reads_assign_dict = args.unphased_reads_assign_dict
read_hp_og_dict = args.read_hp_og_dict
n_thread = args.n_thread

## set logger

logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger("Prepare phasing information")
logger.info("Start extracting information")
# ###################################################


#            loop over bam file


# ###################################################

samfile = pysam.AlignmentFile(bam_path)
samiter = samfile.fetch()
pb_inf = {}
uncover_dict = {}
read_hp_og={}
cnt =0
logger.info("Start looping over bam file")
for read in samiter:
	cnt+=1
	if cnt%1e4==0:
		logger.info("process %d reads"%cnt)
	tags = read.get_tags()
	phased_flag = 0
	hp = -1
	pb = -1
	for tag in tags:
		if tag[0]=='HP':
			phased_flag=1
			hp = tag[1]
		if tag[0]=='PS':
			pb = tag[1]
	if phased_flag:
		if pb in pb_inf:
			pb_inf[pb].append(read.pos)
		else:
			pb_inf[pb] = [read.pos]
		read_hp_og[read.qname] = (pb,hp)
	else:
		uncover_dict[read.qname]=[read.pos,read.reference_end]
logger.info("End looping")
pb_list = []
for pb in pb_inf:
	poss = pb_inf[pb]
	pb_inf[pb]=(min(poss),max(poss))
	pb_list.append((pb,min(poss),max(poss)))

				

df = pd.DataFrame(pb_list,columns = ['name','start','end'])
df.to_csv(pb_info_dict.replace(".p",".csv"),index=False)
pickle.dump(pb_inf,open(pb_info_dict,'wb'))

## map read_hp to real hp name

for read in tqdm(list(read_hp_og.keys()), desc = "Map read_hp_og",miniters =1):
	pb,hp = read_hp_og[read]
	pb_start = pb_inf[pb][0]
	pb_end = pb_inf[pb][1]
	hp = 'PS%d_%d_%d_hp%d'%(pb,pb_start,pb_end,hp)
	read_hp_og[read]=hp
pickle.dump(read_hp_og,open(read_hp_og_dict,'wb'))

# ###################################################


#            assign unphased reads


# ###################################################
# find 2 most close phase blocks for each uncovered read
def assign_unphased(read_start, read_end, pb_inf):
	dist_list = []
	pb_dist_dict = {}
	for pb in pb_inf:
		pb_start = pb_inf[pb][0]
		pb_end = pb_inf[pb][1]
		dist = max(pb_start,read_start)-min(pb_end,read_end)
		dist_list.append([pb,dist])
		pb_dist_dict[pb]=dist
	df_dist = pd.DataFrame(dist_list,columns = ['pb','dist'])
	df_dist_sorted = df_dist.sort_values('dist').reset_index(drop=True)
	if df_dist_sorted.shape[0]>1:
		asign_pbs = df_dist_sorted.pb[:2].tolist()
	else:
		asign_pbs = [df_dist_sorted.pb[0]]
	asign_list = []

	for pb in asign_pbs:
		asign_info = [int(pb), pb_dist_dict[pb],read_start,read_end,pb_inf[pb][0],pb_inf[pb][1]]
		asign_list.append(asign_info)

	if asign_list[0][0]<asign_list[1][0]:
		pb1 = asign_list[0]
		pb2 = asign_list[1]
	else:
		pb1 = asign_list[1]
		pb2 = asign_list[0]

	pb1s = "PS%d_%d_%d"%(pb1[0],pb1[-2],pb1[-1])
	pb2s = "PS%d_%d_%d"%(pb2[0],pb2[-2],pb2[-1])

	pbs = pb1s+'_and_'+pb2s

	return pbs



all_assign = Parallel(n_jobs=n_thread)(delayed(assign_unphased)(uncover_dict[read][0], uncover_dict[read][1], pb_inf)\
 for read in tqdm(list(uncover_dict.keys()), desc="Roughly assign unphased reads",miniters=1))
read_asign_dict = {}
i =0 
for read in uncover_dict:
	read_asign_dict[read] = all_assign[i]
	i+=1

pickle.dump(read_asign_dict,open(unphased_reads_assign_dict ,'wb'))



