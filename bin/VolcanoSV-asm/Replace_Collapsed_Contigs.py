from subprocess import Popen

def load_hap_names(filename):
    hap_names = []
    with open(filename, 'r') as file:
        for line in file:
            hap_names.append(line.strip())  # Remove leading/trailing whitespaces     
    return hap_names

def replace_col_contigs(og_file, new_file, out_file, hap_file):
    haps = set(load_hap_names(hap_file))
    with open(out_file,'w') as fw:
        with open(og_file,'r') as fin:
            for line in fin:
                if line[0] == '>':
                    w = 1

                    hap = '_'.join(line[1:-1].split('_')[:-1])

                    if hap in haps:
                        w = 0

                if w:
                    fw.write(line)

    cmd = f"cat {new_file} >> {out_file}"
    Popen(cmd, shell = True).wait()

import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--og_fasta','-og')
parser.add_argument('--new_fasta','-new')
parser.add_argument('--output_path','-o')
parser.add_argument('--hap_file','-hap')

args = parser.parse_args()
og_fasta = args.og_fasta
new_fasta = args.new_fasta
output_path = args.output_path
hap_file = args.hap_file

replace_col_contigs(og_fasta, new_fasta, output_path, hap_file)