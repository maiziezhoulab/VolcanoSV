import sys
from collections import defaultdict
import os 

# Define the forced order of header keys
HEADER_ORDER = [
    "##fileformat",
    "##contig",
    "##ALT",
    "##INFO",
    "##FILTER",
    "##FORMAT",
]

def parse_vcf(file):
    """Parse a VCF file into headers and records."""
    headers = []
    records = []

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line)
            else:
                records.append(line)
    
    return headers, records

def extract_contig_lines(headers):
    """Extract contig lines from VCF headers."""
    contigs = []
    for header in headers:
        if header.startswith("##contig"):
            contigs.append(header)
    return contigs

def deep_organize(headers):
    """Collect headers by their key (up to the first ',' sign)."""
    header_dict = defaultdict(list)
    for header in headers:
        key = header.split(',', 1)[0]
        header_dict[key].append(header)
        
    return header_dict
    
    
def collect_headers_deep(headers):
    """Collect headers by their key (up to the first '=' sign)."""
    header_dict = defaultdict(list)
    for header in headers:
        key = header.split('=', 1)[0]
        header_dict[key].append(header)
        
	# deep organize
    for key in header_dict:
        header_dict[key] = deep_organize(header_dict[key])
    return header_dict

def collect_headers(headers):
    """Collect headers by their key (up to the first '=' sign)."""
    header_dict = defaultdict(list)
    for header in headers:
        key = header.split('=', 1)[0]
        header_dict[key].append(header)
        
    return header_dict

def merge_contigs(contig_lists):
	"""Merge contig lines, taking the file with the most contigs as the basis."""
	# Find the file with the most contigs
	longest_contig_list = max(contig_lists, key=len)
	# flatten contigs
	all_contigs = []
	for contigs in contig_lists:
		all_contigs.extend(contigs)

	extra_contigs = set(all_contigs) - set(longest_contig_list)

	for contig in extra_contigs:
		longest_contig_list.append(contig)

	return longest_contig_list

def merge_headers(header_lists):
	"""Merge headers from multiple VCF files according to the forced order."""

	header_dicts = [collect_headers(headers) for headers in header_lists]

	# flatten headers
	all_header_lines = []
	for headers in header_lists:
		all_header_lines.extend(headers)
			
	all_header_lines = list(set(all_header_lines))

	header_dict_deep = collect_headers_deep( all_header_lines)

	merged_headers = []

	# Forced order processing
	for key in HEADER_ORDER:
		if key == "##contig":
			# Merge contig lines specially
			contig_lists = [header_dict[key] for header_dict in header_dicts if key in header_dict]
			merged_contigs = merge_contigs(contig_lists)
			merged_headers.extend(merged_contigs)
		else:
			# For other headers, just take the first occurrence and merge
			if key in header_dict_deep:
				dict_key = header_dict_deep[key]
				for key_ID in dict_key:
					vals = dict_key[key_ID]
					merged_headers.append(vals[0])
					if len(vals)>1:
						print(f"\nLogging: Duplicate {key_ID} found, keeping only the first one.")  # Logging duplicate headers
						for val in vals:
							print(val[:-1])


	return merged_headers

def merge_vcfs(vcf_files, output_file):
	"""Merge multiple VCF files into one."""
	header_lists = []
	all_records = []

	# Parse headers and records from each VCF
	for vcf_file in vcf_files:
		headers, records = parse_vcf(vcf_file)
		header_lists.append(headers)
		all_records.extend(records)

	chrom_line = headers[-1]

	# Merge headers using the forced order
	merged_headers = merge_headers(header_lists)

	# Write the merged VCF file
	with open(output_file+".temp", 'w') as f:
		# Write merged headers
		for header in merged_headers:
			f.write(header)
		
		f.write(chrom_line)
		# Write records
		for record in all_records:
			f.write(record)
			
	# sort file
	cmd = f"vcf-sort {output_file}.temp > {output_file}; rm {output_file}.temp"
	os.system(cmd)

	print(f"VCF files merged successfully into {output_file}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python merge_vcf.py output.vcf input1.vcf input2.vcf [inputN.vcf ...]")
        sys.exit(1)
    
    output_file = sys.argv[1]
    vcf_files = sys.argv[2:]
    
    merge_vcfs(vcf_files, output_file)

if __name__ == "__main__":
    main()
