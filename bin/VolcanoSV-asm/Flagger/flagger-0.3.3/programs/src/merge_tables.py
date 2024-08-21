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

def bad_table(values):
	adjusted_values = []
	bad = 0
	# Iterate through the list of values
	for i in range(len(values)):
		# If the current value exceeds the threshold
		if i > 0 and i < len(values) - 1:
			if values[i] > 3 * values[i - 1] and values[i] > 3 * values[i + 1]:
				bad = 1
				break


	return bad

with open(pb_name_file,'r') as f:
	pbs = set(f.read().split('\n')[:-1])

import os, glob 
import pandas as pd

dfs = []
a = 0
b = 0
for in_file in glob.glob(table_dir+"/*table"):
	pb = '_'.join(os.path.basename(in_file).split('.')[-2].split('_')[:-4])
	if pb in pbs:
		df = pd.read_csv(in_file, sep = '\t')

		# # Calculate the sum of the specified columns for each row
		# row_sums = df[['error', 'duplicated', 'haploid', 'collapsed']].sum(axis=1)

		# # Normalize each row by dividing each element by the row sum
		# df = df.div(row_sums, axis=0)

		# df['fit'] = adjust_values(df['fit'])
		# if (df['fit'] > 10**12).sum() == 0:
		if not bad_table(df['fit']):
			# print('=================\n',df)
			dfs.append(df)
		else:
			a+=1
		b+=1
print(a,' bad tables out of ',b )

df = pd.concat(dfs, axis = 0).reset_index(drop = True)
print(df.head())
dfb = df.groupby('#coverage')[['freq', 'fit', 'error', 'duplicated',	'haploid',	'collapsed']].max()
print(dfb)
# df.columns = ['#coverage','raw_freq','raw_fit','raw_error',	'raw_duplicated',	'raw_haploid',	'raw_collapsed']
# List of columns to be multiplied by 'fit'
columns_to_multiply = ['error', 'duplicated', 'haploid', 'collapsed']

# Multiply each of the specified columns by the 'fit' column and replace the original values
df[columns_to_multiply] = df[columns_to_multiply].multiply(df['fit'], axis=0)
print(df.head())

dfa = df.groupby('#coverage')[['freq', 'fit', 'error', 'duplicated',	'haploid',	'collapsed']].max()
print(dfa)

df = df.groupby('#coverage')[['freq', 'fit', 'error', 'duplicated',	'haploid',	'collapsed']].sum()
df.reset_index(inplace=True)
print(df.head())
df[columns_to_multiply] = df[columns_to_multiply].div(df['fit'], axis=0)


# Calculate the sum of the specified columns for each row
row_sums = df[['error', 'duplicated', 'haploid', 'collapsed']].sum(axis=1)

# Normalize each row by dividing each element by the row sum
df = df.div(row_sums, axis=0)



df.to_csv(outfile, sep = '\t', index = False)
