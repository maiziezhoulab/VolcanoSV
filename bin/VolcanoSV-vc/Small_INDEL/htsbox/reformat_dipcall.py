import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path


def write_one_var(data,ins_cnt, del_cnt, snp_cnt, fw):          
    chrom = data[0]
    pos = int(data[1])
    if 'SVTYPE=INS' in data[7]:
        svtype = 'INS'
        ins_cnt += 1
        svidx = '%s-%s-%d-%d'%(chrom, svtype, ins_cnt ,pos)
    elif 'SVTYPE=DEL' in data[7]:
        del_cnt +=1 
        svtype = "DEL"
        svidx = '%s-%s-%d-%d'%(chrom,svtype,  del_cnt,pos, )
    else:
        snp_cnt +=1
        svtype = 'SNP'
    
        svidx = '%s-%s-%d-%d'%( chrom,svtype, snp_cnt,pos, )

    data[2] = svidx
    fw.write('\t'.join(data)+'\n')
    return ins_cnt,del_cnt,snp_cnt

with open(output_path,'w') as fw:
    ins_cnt = 0
    del_cnt = 0
    snp_cnt = 0
    with open(input_path , 'r') as f:
        for line in f:

            if line[0]=='#':
                fw.write(line)
            else:
                data = line.split()
                if ',' not in data[4]:
                    ins_cnt,del_cnt,snp_cnt = write_one_var(data,ins_cnt, del_cnt, snp_cnt, fw)
                else:
                    alt_list = data[4].split(',')
                    info_list = data[7].split(',')
                    for i in range(len(alt_list)):
                        if info_list[i]!='*':
                            data1 = data.copy()
                            data1[4] = alt_list[i]
                            data1[7] = info_list[i]
                            ins_cnt,del_cnt,snp_cnt = write_one_var(data1,ins_cnt, del_cnt, snp_cnt, fw)





