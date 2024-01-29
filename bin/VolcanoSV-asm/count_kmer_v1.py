from joblib import Parallel, delayed
from collections import Counter
from tqdm import tqdm
import pickle
import os
from argparse import ArgumentParser
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--read_hp_og_path','-rho')
parser.add_argument('--hash_path','-hash')
parser.add_argument('--name_path','-name')
parser.add_argument('--output_dir','-o')
parser.add_argument('--batch_size','-b',type = int, default = 10000)
parser.add_argument('--n_threads','-t', type = int, default = 50)
args = parser.parse_args()
read_hp_og_path = args.read_hp_og_path
hash_path = args.hash_path
name_path = args.name_path
output_dir = args.output_dir
batch_size = args.batch_size
n_threads = args.n_threads


import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


def count_kmer(line):
    kmers = [int(k) for k in line.split()]
    return Counter(kmers)

def count_kmer_one_hap(lines):
    kmers = []
    for line in lines:
        kmers.extend([int(k) for k in line.split()])
    return Counter(kmers)

def reduce_kcs(kcs,hps):
    kc = kcs[0]
    dc = {hps[0]:kc}
    for i in range(1,len(kcs)):
        hp = hps[i]
        kcnew = kcs[i]
        if hp in dc:
            kc = dc[hp] 
            for k in kcnew:
                if k in kc:
                    kc[k]+=kcnew[k]
                else:
                    kc[k] = kcnew[k]
        else:
            dc[hp] = kcnew
    return dc

def update_dc_all(dc_all,dc):
    for hp,kcnew in dc.items():
        if hp in dc_all:
            kc = dc_all[hp]
            for k in kcnew:
                if k in kc:
                    kc[k]+=kcnew[k]
                else:
                    kc[k] = kcnew[k]
        else:
            dc_all[hp] = kcnew
    return dc_all








## load read hap dict
read_hp_og_dc = pickle.load(open(read_hp_og_path,'rb' ))

## load name file
with open(name_path,'r') as f:
    name_list = f.read().split('\n')[:-1]

## load hash file
fin = open(hash_path)
idx = 0
## collect lines
cnt = 0
lines = []
hps = []
dc_hp = {}
for line in fin:
    if name_list[idx] in read_hp_og_dc:
        cnt +=1
        hp = read_hp_og_dc[name_list[idx]]
        lines.append(line)
        hps.append(hp)
        if hp in dc_hp:
            dc_hp[hp].append(line)
        else:
            dc_hp[hp] = [line]
    idx+=1
    if cnt>=batch_size:
        break

dc_all = {}
batch_id = 0
while len(lines):
    logger.info("process batch %d"%batch_id)

    ## map
    kcs= Parallel(n_jobs=n_threads)(delayed(count_kmer_one_hap)(lines) for lines in dc_hp.values())

    ## reduce
    dc = dict(zip(list(dc_hp.keys()), kcs))

    # kcs= Parallel(n_jobs=n_threads)(delayed(count_kmer)(line) for line in tqdm(lines, desc = "map"))
    # dc = reduce_kcs(kcs,hps)
    dc_all = update_dc_all(dc_all,dc)
    logger.info("number of haps: %d"%len(dc_all))

    ## collect lines
    cnt = 0
    lines = []
    hps = []
    dc_hp = {}
    for line in fin:
        if name_list[idx] in read_hp_og_dc:
            cnt +=1
            hp = read_hp_og_dc[name_list[idx]]
            lines.append(line)
            hps.append(hp)
            if hp in dc_hp:
                dc_hp[hp].append(line)
            else:
                dc_hp[hp] = [line]
        idx+=1
        if cnt>=batch_size:
            break

    ## update batch id
    batch_id+=1

fin.close()

## dump kmer db
logger.info("write kmer database...")
os.system("mkdir -p "+output_dir )
for hap in tqdm(dc_all):
    output_path = output_dir+'/'+hap+'_kc.p'
    pickle.dump(dc_all[hap],open(output_path,'wb'))



















