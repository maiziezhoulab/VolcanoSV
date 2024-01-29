from argparse import ArgumentParser
import pickle
from tqdm import tqdm
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--name_file','-name')
parser.add_argument('--read_hp_dc_path','-dc')
parser.add_argument('--output_dir','-o')
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
name_file = args.name_file
read_hp_dc_path = args.read_hp_dc_path
output_dir = args.output_dir

import logging
## set logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

logger.info("load dict")
dc = pickle.load(open(read_hp_dc_path,'rb'))
logger.info("load name file")
with open(name_file,'r') as f:

    name_list = f.readlines()


logger.info("start splitting")
os.system("rm -r "+output_dir)
os.system("mkdir -p "+output_dir)
fh = open(input_path,'r')
for i in tqdm(range(len(name_list))):
    hash_line = fh.readline()
    read = name_list[i][:-1]
    if read in dc:
        hp = dc[read]
        with open(output_dir+'/%s.hash'%hp,'a+') as f:
            f.write(hash_line)
        with open(output_dir+'/%s.name'%hp,'a+') as f:
            f.write(name_list[i])














