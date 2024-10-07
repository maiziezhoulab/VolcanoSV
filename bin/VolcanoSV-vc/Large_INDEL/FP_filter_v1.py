import numpy as np
from argparse import ArgumentParser
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--signature_dir','-sigd')
parser.add_argument('--output_path','-o')
# parser.add_argument('--chr_number','-chr',type = int)
parser.add_argument('--max_comp_svlen','-max_comp_svlen', type = int, default = 250)
parser.add_argument('--max_dist','-max_dist', type = int, default = 1000)
parser.add_argument('--max_shift','-max_shift', type = int, default = 500)
parser.add_argument('--min_size_sim','-min_size_sim', type = float, default = 0.5)
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
signature_dir = args.signature_dir
output_path = args.output_path
# chr_number = args.chr_number
max_comp_svlen = args.max_comp_svlen
max_dist = args.max_dist
max_shift = args.max_shift
min_size_sim = args.min_size_sim


# max_comp_svlen =250
# max_dist = 1000
# max_shift = 500
# min_size_sim = 0.5

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


def vcf_to_sig(vcf_path,target_chr_name):
    with open(vcf_path,'r') as f:
        s = f.readlines()
    sig_list =[]
        
    for line in s:
        if line[0]!='#':
            data = line.split()
            chr_name = data[0]
            pos = int(data[1])
            svlen = len(data[4])-len(data[3])
            if chr_name ==target_chr_name:
                if svlen<0:
                    svtype = 'DEL'
                else:
                    svtype = 'INS'
                sig = [chr_name, svtype, pos, abs(svlen)]
                sig_list.append(sig)
    return sig_list

def vcflines_to_sig(vcf_lines,target_chr_name):
    s = vcf_lines
    sig_list =[]
        
    for line in s:
        if line[0]!='#':
            data = line.split()
            chr_name = data[0]
            pos = int(data[1])
            svlen = len(data[4])-len(data[3])
            if chr_name ==target_chr_name:
                if svlen<0:
                    svtype = 'DEL'
                else:
                    svtype = 'INS'
                sig = [chr_name, svtype, pos, abs(svlen)]
                sig_list.append(sig)
    return sig_list

def load_sig(sig_path):
    with open(sig_path,'r') as f:
        s = f.readlines()
    sig_list =[]
    
    for line in s:
        data = line.split()
        data[2],data[3] = int(data[2]),int(data[3])
        sig_list.append(data)
    return sig_list
def compare_sigs(sig1,sig2,max_shift = 500,min_size_sim = 0.3):
#     max_shift = 100
    
    
    
    shift = abs(sig1[2]-sig2[2])
    # shift_ratio = shift/min(sig1[3],sig2[3])
    try:
        size_sim = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
    except:
        size_sim = 0
    
    if (shift<=max_shift) and (size_sim >= min_size_sim):
        return 1
    else:
        return 0



def eval_sig(sig_list, reads_sig_list,max_dist,max_comp_svlen =300,max_shift = 500,min_size_sim = 0.3) :
    support_list = []
    
    for sig1 in sig_list:
        if sig1[3]>max_comp_svlen:
            support_list.append(60)
        else:
            support=0
            for sig2 in reads_sig_list:
                shift = sig2[2]-sig1[2]
                if shift<-max_dist:
                    continue
                if shift>max_dist:
                    break
                elif compare_sigs(sig1,sig2,max_shift,min_size_sim )==1:
                    support+=1
            support_list.append(support)
    return support_list
    
                
def fp_filter(sig_list_tp,sig_list_fp,reads_sig_list,max_dist,max_comp_svlen,max_shift,min_size_sim):
    support_list_tp = eval_sig(sig_list_tp, reads_sig_list ,max_dist,max_comp_svlen,max_shift,min_size_sim)
    reduce_tp = (np.array(support_list_tp)<1).sum()
    print("Reduce TP %d"%reduce_tp)
    support_list_fp = eval_sig(sig_list_fp, reads_sig_list ,max_dist,max_comp_svlen,max_shift,min_size_sim)
    reduce_fp = (np.array(support_list_fp)<1).sum()
    print("Reduce FP %d"%reduce_fp)
    return reduce_tp,reduce_fp

def filter_vcf(vcf_lines,chr_name,sig_path,
               max_dist,max_comp_svlen,max_shift,min_size_sim):
    sig_chr = vcflines_to_sig(vcf_lines,chr_name)
    reads_sig_list = load_sig(sig_path)
            
    support_list = eval_sig(sig_chr, reads_sig_list ,max_dist,max_comp_svlen,max_shift,min_size_sim)
    
    assert len(support_list)==len(vcf_lines)
    support_list = np.array(support_list)
    idxs = np.where(support_list>0)[0]
            
    print("reduced %d lines"%((support_list==0).sum()))
    return np.array(vcf_lines)[idxs]
          
def load_wgs_vcf(vcf_path):
    with open(vcf_path,'r') as f:
        s = f.readlines()
    header = []
    dc ={}
    for line in s:
        if line[0]=='#':
            header.append(line)
        else:
            name = line.split()[0]
            if name not in dc:
                dc[name]=[line]
            else:
                dc[name].append(line)
    return header,dc
                
                
                       

vcf_path = input_path
header,dc = load_wgs_vcf(vcf_path)

final_lines = []
for i in range(1,23):
    chr_name = 'chr%d'%i
    if chr_name in dc:
        print(chr_name)
        sig_path = signature_dir + "/%s_reads_sig.txt"%chr_name
        vcf_lines_chr = filter_vcf(dc[chr_name],chr_name,sig_path,
                       max_dist,max_comp_svlen,max_shift,min_size_sim)
        final_lines.append(vcf_lines_chr)
    else:
        final_lines.append([])



out_path = output_path
fw = open(out_path,'w')
fw.writelines(header)
reduced_line  = 0
for i in range(22):
#     idxs = np.where(final_lines[i]>0)[0]
    chr_name = 'chr%d'%(i+1)
    if chr_name in dc:
        print(chr_name,len(dc[chr_name]),len(final_lines[i]))
        reduced_line += (len(dc[chr_name])-len(final_lines[i]))
    #     vcf_lines = list(np.array(dc[chr_name])[idxs])
        fw.writelines(final_lines[i])
fw.close()





print("final reduced lines %d"%reduced_line)