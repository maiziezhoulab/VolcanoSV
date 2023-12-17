# :milky_way: VolcanoSV 🌋


## Dependencies for Github installation:
RegionIndel utilizes <a href="https://www.python.org/downloads/">Python3 (+ numpy, pysam, sortedcontainers, and scipy)</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to check their availability and install them if not, but make sure you have installed "python3", "conda" and "wget" first. 

# Install through Github:

```
git clone https://github.com/maiziex/RegionIndel.git
cd RegionIndel
chmod +x install.sh
./install.sh
```


## source folder:
After running "./install.sh", a folder "source" would be downloaded.

## Running The Code:
Put the "RegionIndel/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or just use the full path of "**RegionIndel_step1.py**", "**RegionIndel_step2.py**" and "**RegionIndel_step3.py**"

*We provide  <a href="https://github.com/maiziezhoulab/RegionIndel/blob/main/example_data/run_example_data.md">a test example dataset</a> to run the whole pipeline. 

### Step 0:
We added "orphan end reads" (OER) to our pipeline to generate more accurate assembly and boost the performance. Orphan end reads are defined as a read pair in which each mate is aligned to a different chromosome. The below code is how we prepare OER in a whole genome scale.
```
## specify the wgs bam path
wgs_bam=possorted_bam.bam

# create directory for by chromosme bam files
mkdir bam_by_chr

# extract bam file for chr1-chr22 and hs37d5, and add index
for i in {1..22}
do
samtools view -Sb $wgs_bam chr$i > bam_by_chr/possorted_bam_chr${i}.bam
samtools index bam_by_chr/possorted_bam_chr${i}.bam
done

samtools view -Sb $wgs_bam hs37d5 > bam_by_chr/possorted_bam_hs37d5.bam
samtools index bam_by_chr/possorted_bam_hs37d5.bam

# extract OER part1
for i in {1..22}
do
python3 RegionIndel/bin/OER_SCAN_part1.py \
-i bam_by_chr/possorted_bam_chr${i}.bam \
-o ./OER/part1
done

# extract OER part2
for i in {1..22}
do
python3 RegionIndel/bin/OER_SCAN_part2.py \
-i bam_by_chr/possorted_bam_chr${i}.bam \
-chrn chr${i} -o ./OER/part2 \
-rq ./OER/part1
done

python3 RegionIndel/bin/OER_SCAN_part2.py \
-i bam_by_chr/possorted_bam_hs37d5.bam \
-chrn hs37d5 -o ./OER/part2 \
-rq ./OER/part1
```
After running the above code, you will have output folder "./OER/part2" that contains important OER information. Now, you are good to run the RegionIndel pipeline.


### Step 1: 
```
python3 RegionIndel/bin/RegionIndel_step1.py  --bam_file selected.bam --vcf_file test_freebayes.vcf --chr_num 3 --out_dir test_sv --OER_dir ./OER/part2

```
#### *Required parameters

**--bam_file:** "selected.bam" is a bam file for a target region. How to get the bam file, you can also check <a href="https://github.com/maiziezhoulab/RegionIndel/blob/master/src/How_to_get_bam_and_vcf.md">here</a>.

**--vcf_file:** "test_freebayes.vcf" is a VCF file generated from variant caller like "FreeBayes". How to get the vcf file, you can also check <a href="https://github.com/maiziezhoulab/RegionIndel/blob/master/src/How_to_get_bam_and_vcf.md">here</a>. 

**--chr_num:** "3" is the chromosome number you need to define for the target region or the structural variant you are interested in.

**--OER_dir:**  ./OER/part2 is the OER folder you generated in step 0.


#### *Optional parameters
**--mole_boundary:** default = 50000 (50kb). We use 50kb to differentiate reads with the same barcode are drawn from different long molecules. 

**--out_dir:** default = ./RegionIndel_results. You can define your own folder name.

**--num_threads:** default = 8. 

**--num_threads_bwa_mem:** number of threads for bwa-mem, default = 20

**--clean:** default = 1. It will delete all assembly files from SPAdes and intermediate bam/fastq files from RegionIndel.



### Step 2: 
```
python3 RegionIndel/bin/RegionIndel_step2.py --out_dir test_sv --chr_num 3 

```
#### *Required parameters
**--chr_num:** "3" is the chromosome number you need to define for the target region or the structural variant you are interested in.

**--reference:** "genome_hg19.fa" is the human reference fasta file.

#### *Optional parameters
**--out_dir:** default = ./RegionIndel_results, make sure it's the same as "--out_dir" from ***Step1*** if you want to define your own output directory name.

**--num_threads:** default = 10, this determines the number of files assembled simultaneously by SPAdes.  

**--num_threads_spades:** default = 5, this is the "-t" for SPAdes. 



### Step 3: 
```
python3 RegionIndel/bin/RegionIndel_step3.py  --assembly_dir test_sv  -o_dir test_sv --ref_file genome_hg19.fa  --chr_num 3 

```
#### *Required parameters
**--assembly_dir:** folder to store assembly results from step1 and step2 (same as "--out_dir" for step1 and step2).

**--ref_file:** "genome_hg19.fa" is the human reference fasta file.

**--chr_num:** "3" is the chromosome number you need to define for the target region or the structural variant you are interested in.

#### *Optional parameters
**--out_dir:** default = ./RegionIndel_Step3_Results. Directory to store the final VCF file from RegionIndel, and you can define your own folder name

**--var_size:** default = 1, variant size, cut off size for indel and SV, 

**--num_of_threads:** number of threads, default = 1

**--clean:** default = 1. You can choose to delete intermediate files or no

### Step 4 (optional): 
```
#----------extract SVs
cat ./test_sv/RegionIndel_Step3_Result/RegionIndel_Contig_final_sorted.vcf \
	| awk '($1 ~ /^#/ || length($5) - length($(4)) > 30 || length($4) - length($(5)) > 30 )' \
	> ./test_sv/RegionIndel_Step3_Result/RegionIndel_Contig_final_sorted_sv.vcf

python3 RegionIndel/bin/remove_redundancy.py   \
-i ./test_sv/RegionIndel_Step3_Result/RegionIndel_Contig_final_sorted_sv.vcf  \
-o ./test_sv/Remove_redundancy/

```

Sometimes the result will have duplicate calls, and to remove these calls, you can use **remove_redundancy.py**. 

## Running for multiple regions:

If you have multiple interested regions, you can save the regions into a BED file and use our example script **bin/mt_region_RegionIndel.sh** to detect SVs in multiple regions. Make sure you change the parameters in the bash script before you run it. The parameters are explained below:

```
#-----------------Parameter------------

bed_file="test.bed"  # Specify your BED file
wgsbam=possorted_bam.bam #your whole genome bam file
wgsvcf=possorted_bam_freebayes.vcf #your whole genome VCF file
oerdir=./OER/part2   # the OER file generated in step 0
reference=genome_hg19.fa # reference file
RegionInde_dir=./RegionIndel/ # RegionIndel install directory (The full path the user installs RegionIndel)

```


#### Memory/Time Usage For RegionIndel
Coverage| Memory| Time for one SV on a single node 
--- | --- | --- | 
60X | 20GB | 00:10:32 |


## Final Output:
**test_sv/RegionIndel_step3_results:** RegionIndel_contig_final.vcf
```
test_sv
|
|-H5_for_molecules 
|   └-target_chr3_sorted.h5    --> (molecule files for target region including barcode, variants annotation (0: ref allele; 1: alt allele), coordinates for each molecule)
|
|-results_phased_probmodel
|   └-chr*.phased_final        --> (Phased molecule files)
|
|
|-Local_Assembly_by_chunks
|   └-chr*_files_cutPBHC
|       |-fastq_by_*_*_hp1.fastq                  --> (reads fastq file for haplotype 1)
|       |-fastq_by_*_*_hp2.fastq                  --> (reads fastq file for haplotype 2)
|       
|
|-Assembly_Contigs_files
|    |-RegionIndel_Contig_chr*.fasta                    --> (final contigs fasta file for the target region)
|    |-RegionIndel_Contig_chr*.bam                      --> (final contigs bam file for the target region)
|    |-RegionIndel_Contig_chr*_hp1.fasta                     --> (final contigs fasta file for haplotype 1)
|    └-RegionIndel_Contig_chr*_hp2.fasta                     --> (final contigs fasta file for haplotype 2)
|
└-RegionIndel_step3_results
     └-RegionIndel_Contig_final_sorted.vcf --> (final VCF including SNPs, small Indels and SVs)
     
     
```

## Final contig bam file (RegionIndel_Contig_chr*.bam) diplayed in IGV (SV + 25kb left and right flanking regions):
<p align="center">
	<img src="src/igv1.png"  width="650">
</p>


## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziezhoulab/RegionIndel/issues">RegionIndel</a>. 


