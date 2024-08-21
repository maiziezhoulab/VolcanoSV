import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--final_bed_file','-bed')
parser.add_argument('--pb_name_file','-pb')
parser.add_argument('--output_path','-o')


args = parser.parse_args()


pb_name_file = args.pb_name_file
final_bed_file = args.final_bed_file
outfile = args.output_path



import pandas as pd 
from subprocess import Popen

with open(pb_name_file,'r')as f:
	pbs = set(f.read().split('\n')[:-1])

new_bed_file = final_bed_file.replace(".bed",'.csv')
cmd = f"grep -v track {final_bed_file}|cut -f 1,2,3,4 > {new_bed_file}"
Popen(cmd, shell= True).wait()
# print(new_bed_file)
df = pd.read_csv(new_bed_file, header = None, names = ['contig','start','end','label'], sep = '\t')
# print(df.shape)
df['size'] = df['end'] - df['start']

df['pb'] = df['contig'].apply(lambda x: '_'.join(x.split('_')[:-2]))
# print(df.head())
df_subset = df[df['pb'].isin(pbs)].reset_index(drop=True)

df_agg = df_subset.groupby('label')['size'].sum().reset_index()


total_size = df_agg['size'].sum()

# Create a new row with the sum
sum_row = pd.DataFrame({'label': ['Total'], 'size': [total_size]})

# Append the sum row to the original DataFrame
df = pd.concat([df_agg, sum_row], ignore_index=True)
df['perct'] = round(df['size']/total_size * 100, 3)
df.to_csv(outfile, index = False)
