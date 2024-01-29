import pandas as pd
from collections import Counter
import os
from tqdm import tqdm 
import pickle
from argparse import ArgumentParser
from joblib import Parallel, delayed
from sklearn import preprocessing
import numpy as np
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--unphased_dir','-ud')
parser.add_argument('--kmer_db','-kdb')
parser.add_argument('--output_dir','-o')
parser.add_argument('--prefix','-px')
parser.add_argument('--n_thread','-t', type = int, default = 10)
parser.add_argument('--significance_level','-sigl',type = float, default = 0.1)
args = parser.parse_args()
unphased_dir = args.unphased_dir
kmer_db =  args.kmer_db
output_dir = args.output_dir
prefix = args.prefix
n_thread = args.n_thread
significance_level = args.significance_level


import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")



def load_db(db_path):
    try:
        dc = pickle.load(open(db_path,'rb'))
        return dc
    except:
        return {}

def count_overlap(R,H):
    cnt = 0
    for k in R:
        if k in H:
            cnt+= (R[k] + H[k])
            # cnt+= R[k]
    return cnt

def filter_dc(dc,K):
    dc1 ={}
    for k in K:
        dc1[k]=dc[k]
    return dc1 

def get_unique_kmer(dc1,dc2,dc3,dc4):
    k1 = set(dc1.keys())
    k2 = set(dc2.keys())
    k3 = set(dc3.keys())
    k4 = set(dc4.keys())
    K1 = k1-k2-k3-k4
    K2 = k2-k1-k3-k4
    K3 = k3-k1-k2-k4
    K4 = k4-k1-k2-k3

    dc11 = filter_dc(dc1,K1)
    dc21 = filter_dc(dc2,K2)
    dc31 = filter_dc(dc3,K3)
    dc41 = filter_dc(dc4,K4)
    return dc11,dc21,dc31,dc41

def get_raw_overlap_count(hash_file,hp_dir):
    name_file = hash_file.replace(".hash",".name")
    pb1 = hash_file.split('/')[-1].split('_and_')[0]
    pb2 = hash_file.split('/')[-1].split('_and_')[1].split('.hash')[0]
    p1 = hp_dir+'/'+pb1+'_hp1_kc.p'
    p2 = hp_dir+'/'+pb1+'_hp2_kc.p'
    p3 = hp_dir+'/'+pb2+'_hp1_kc.p'
    p4 = hp_dir+'/'+pb2+'_hp2_kc.p'
    dc1 = load_db(p1)
    dc2 = load_db(p2)
    dc3 = load_db(p3)
    dc4 = load_db(p4)

    dc1,dc2,dc3,dc4 = get_unique_kmer(dc1,dc2,dc3,dc4)

    f = open(hash_file)
    ks = []
    for line in f:
        x = [int(cnt) for cnt in line.split()]
        dc5 = Counter(x)
        k1 = count_overlap(dc5,dc1)
        k2 = count_overlap(dc5,dc2)
        k3 = count_overlap(dc5,dc3)
        k4 = count_overlap(dc5,dc4)
        ks.append((k1,k2,k3,k4))
    f.close()
    with open(name_file,'r') as f:
        names = f.read().split('\n')[:-1]
    dc = dict(zip(names,ks))


    return dc,pb1,pb2



up_dir = unphased_dir
hp_dir = kmer_db

os.system("mkdir -p "+output_dir)

hash_files = [ up_dir + '/'+fname for fname in os.listdir(up_dir) if '.hash' in fname]
dcs= Parallel(n_jobs=n_thread)(delayed(get_raw_overlap_count)(hash_file,hp_dir) for hash_file in tqdm(hash_files))
dc_all = {}
for dc1,pb1,pb2 in dcs:
    dc_all[(pb1,pb2)] = dc1 

# dc_all = {}
# for hash_file in tqdm(hash_files):
#     name_file = hash_file.replace(".hash",".name")
#     pb1 = hash_file.split('/')[-1].split('_and_')[0]
#     pb2 = hash_file.split('/')[-1].split('_and_')[1].split('.hash')[0]
#     p1 = hp_dir+'/'+pb1+'_hp1_kmer_count.csv'
#     p2 = hp_dir+'/'+pb1+'_hp2_kmer_count.csv'
#     p3 = hp_dir+'/'+pb2+'_hp1_kmer_count.csv'
#     p4 = hp_dir+'/'+pb2+'_hp2_kmer_count.csv'
#     dc = get_raw_overlap_count(hash_file,name_file,p1,p2,p3,p4)
#     dc_all[(pb1,pb2)] = dc 

out_path = output_dir+'/'+prefix+'_raw_kc.p'
pickle.dump(dc_all,open(out_path,'wb'))


#### normalize data


## reverse data
dc_r_pbs = {}
for pbs,dc in dc_all.items():
    for r,_ in dc.items():
        dc_r_pbs[r] = pbs



dc_asg = {}
dc = dc_all
x = list(dc.keys())
names = []
X = []
for i in range(len(x)):
    dc1 = dc[x[i]]
    for r,ks in dc1.items():
        X.append(list(ks))
        names.append(r)
X = np.array(X)
X = preprocessing.normalize(X)
X_1d = X.flatten()
cut_off = np.quantile(X_1d, 1-significance_level )
logger.info("cut off %.2f"%cut_off)
dct = {0:(0,1),1:(0,1),2:(2,3),3:(2,3)}
# print(len(x))
# print(len(names))
# print(len(maxx))
suc = 0
for i in range(len(X)):
    r = names[i]
    mx = max(X[i])
    if mx >= cut_off:
        suc+=1
        asg = [np.argmax(X[i])]
    else:
        asg = list(dct[np.argmax(X[i])])
    pbs = dc_r_pbs[r]
    hps = [pbs[0]+'_hp1',
    pbs[0]+'_hp2',
    pbs[1]+'_hp1',
    pbs[1]+'_hp2',]

    hps_pick = [ hps[p] for p in asg]

    dc_asg[r] = hps_pick
logger.info("single assignment: %d"%suc)
logger.info("double assignment: %d"%(len(X)-suc))
logger.info("single assignment perct: %.2f"%(suc/len(X)))



out_path = output_dir+'/'+prefix+'_norm_asg.p'
pickle.dump(dc_asg,open(out_path,'wb'))











