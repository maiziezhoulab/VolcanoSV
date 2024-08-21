import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--table_dir','-t')
parser.add_argument('--pb_name_path','-pb')
parser.add_argument('--output_path','-o')

args = parser.parse_args()
pb_name_file = args.pb_name_path
table_dir = args.table_dir
outfile = args.output_path



with open(pb_name_file,'r') as f:
	pbs = set(f.read().split('\n')[:-1])

import os, glob 
import pandas as pd

dfs = []
for in_file in glob.glob(table_dir+"/*table"):
	pb = '_'.join(os.path.basename(in_file).split('.')[-2].split('_')[:-4])
	if pb in pbs:
		df = pd.read_csv(in_file, sep = '\t')
		dfs.append(df)

df = pd.concat(dfs, axis = 0).reset_index(drop = True)

# df.columns = ['#coverage','raw_freq','raw_fit','raw_error',	'raw_duplicated',	'raw_haploid',	'raw_collapsed']

df1 = df.groupby('#coverage')[['freq', 'fit']].sum()

df2 = df.groupby('#coverage')[['error',	'duplicated',	'haploid',	'collapsed']].mean()
df1.reset_index(inplace=True)
df2.reset_index(drop=True, inplace=True)

# Calculate the sum of the specified columns for each row
row_sums = df2[['error', 'duplicated', 'haploid', 'collapsed']].sum(axis=1)

# Normalize each row by dividing each element by the row sum
df2_normalized = df2.div(row_sums, axis=0)

result = pd.concat([df1, df2_normalized], axis=1)



result.to_csv(outfile, sep = '\t', index = False)
