import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcffile','-vcf')
parser.add_argument('--output_dir','-o')
parser.add_argument('--bamfile','-bam')
parser.add_argument('--n_thread','-t', type = int, default = 30 )
parser.add_argument('--flanking','-f',type = int, default = 1000, help = 'flanking region')
parser.add_argument('--min_support','-ms',type = int, default = 10, help = 'minimum support')
parser.add_argument('--min_mapq','-mq',type = int, default = 10, help = 'minimum mapping quality')
parser.add_argument('--max_dist','-d',type = int, default = 100, help = 'maximum distance to merge 2 breakends')

args = parser.parse_args()
vcffile= args.vcffile
output_dir = args.output_dir
bamfile = args.bamfile
n_thread = args.n_thread
flanking = args.flanking
min_support = args.min_support
min_mapq = args.min_mapq
max_dist = args.max_dist


import numpy as np
import os
from collections import defaultdict



def load_raw_vcf(vcffile):
    dc = defaultdict(list)
    header = []
    bnd_list1 = []
    bnd_list2 = []
    with open(vcffile,'r') as f:
        for line in f:
            if (line[0]!='#') and ('SVTYPE=BND' in line):
                data = line.split()
                chrom1 = data[0]
                pos1 = int(data[1])
                svid = data[2]
                chrom2, pos2 = data[4].replace(']','[') .split('[')[1].split(':')
                pos2 = int(pos2)
                if '[' in data[4]:
                    sep = '['
                    bnd_list1.append((chrom1,pos1,chrom2,pos2,sep))
                else:
                    sep = ']'
                    bnd_list2.append((chrom1,pos1,chrom2,pos2,sep))
                dc[(chrom1,pos1,chrom2,pos2,sep)].append(line)
                # dc[svid] = (chrom1,pos1,chrom2,pos2,line)
            elif line[0]=='#':
                header.append(line)

    return dc, header, bnd_list1, bnd_list2

def calculate_center(bnd_list):
    pos1_list = [bnd[1] for bnd in bnd_list]
    pos2_list = [bnd[3] for bnd in bnd_list]

    avg_pos1 = int(sum(pos1_list)/len(pos1_list))
    avg_pos2 = int(sum(pos1_list)/len(pos2_list))
    chrom1= bnd_list[0][0]
    chrom2= bnd_list[0][2]
    return (chrom1, avg_pos1, chrom2, avg_pos2,bnd_list[0][4])


def cluster_bnd(bnd_list, max_dist):
    clusters= [ [ bnd_list[0]]]
    for i in range(1,len(bnd_list)):
        new_bnd = bnd_list[i]
        old_bnd = clusters[-1][-1]
        # print(new_bnd, old_bnd)
        if (( new_bnd[0], new_bnd[2]) == ( old_bnd[0], old_bnd[2])) \
        &  ((new_bnd[1]-old_bnd[1])<=max_dist) &  ((new_bnd[3]-old_bnd[3])<=max_dist) & (new_bnd[4]==old_bnd[4]):
            clusters[-1].append(new_bnd)
        else:
            clusters.append([new_bnd])

    dc = {}

    for bnd_list in clusters:
        center = calculate_center(bnd_list)
        for bnd in bnd_list:
            dc[bnd] = center
    return dc






def merge_bnd(dc,output_path, header, max_dist,bnd_list1, bnd_list2):


    dc_merged = defaultdict(list)
    dc_map1 = cluster_bnd(bnd_list1, max_dist)
    dc_map2 = cluster_bnd(bnd_list2, max_dist)
    dc_map1.update(dc_map2)
    for key,val in dc.items():
        center = dc_map1[key]
        dc_merged[center].extend(val)


    with open(output_path, 'w') as f:
        f.writelines(header)
        for key,val in dc_merged.items():
            if len(val)==1:
                f.write(val[0])
            else:
                data= val[0].split()
                data[-1] = '1/1'
                line = '\t'.join(data)+'\n'
                f.write(line)
                
            




def merge_raw_vcf(vcffile, output_path, max_dist):
    

    if not os.path.exists(output_dir):
        os.system("mkdir -p " + output_dir)

    dc_var, header,bnd_list1, bnd_list2 = load_raw_vcf(vcffile)

    merge_bnd(dc_var,output_path , header, max_dist,bnd_list1, bnd_list2)


    return 




def load_merged_vcf(vcffile):
    dc = {}
    header = []
    with open(vcffile,'r') as f:
        for line in f:
            if (line[0]!='#') and ('SVTYPE=BND' in line):
                data = line.split()
                chrom1 = data[0]
                pos1 = int(data[1])
                svid = data[2]
                chrom2, pos2 = data[4].replace(']','[') .split('[')[1].split(':')
                pos2 = int(pos2)
                dc[svid] = (chrom1,pos1,chrom2,pos2,line)
            elif line[0]=='#':
                header.append(line)

    return dc, header


### create output folder
if not os.path.exists(output_dir):
    os.system("mkdir -p "+output_dir)

###### merge vcf
# merged_vcf = output_dir+'/TRA_merged.vcf'
merged_vcf = output_dir+'/TRA_final.vcf'
merge_raw_vcf(vcffile, merged_vcf , max_dist)






















