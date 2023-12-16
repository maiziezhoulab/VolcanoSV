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
                    fw.write(line)
                else:
                    alt_list = data[4].split(',')
                    for i in range(len(alt_list)):
                        
                        data1 = data.copy()
                        data1[4] = alt_list[i]
                        if not ((len(data[3])==1) & (len(data[4]) ==1)):
                            line = '\t'.join(data1)+'\n'
                            fw.write(line)







