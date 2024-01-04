# :milky_way: VolcanoSV 🌋


## Dependencies for Github installation:
VolcanoSV utilizes Python3.8.6. 
The dependencies include

```
collections
edlib
joblib
math
multiprocessing
pysam
random
scipy
sklearn
statistics
subprocess
tqdm
argparse
edlib
gzip
jellyfish
logging
matplotlib
networkx
numpy
pandas
pickle
pysam
subprocess
pyvcf
```
You need to run `conda install <packagename>` to install all the above dependencies.

# Install through Github:

```
git clone https://github.com/maiziex/VolcanoSV.git
```




## Running The Code:
Put the "VolcanoSV/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or just use the full path of "**volcanosv-asm.py**", "**volcanosv-vc-large-indel.py**", "**volcanosv-vc-complex-sv.py**" and "**volcanosv-vc-small-indel.py**"



### VolcanoSV Assembly 


The VolcanoSV assembly is designed to be run by chromosomes. The main code is `bin/VolcanoSV-asm/volcanosv-asm.py`. The input arguments for this code are explained below:


```
  --inbam INBAM, -i INBAM
  --out_dir OUT_DIR, -o OUT_DIR
  --reference REFERENCE, -r REFERENCE
  --n_thread N_THREAD, -t N_THREAD
  --chrnum {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}, -chr {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
  --dtype {CCS,CLR,ONT}, -d {CCS,CLR,ONT}
  --prefix PREFIX, -px PREFIX

```



After running the above code, you will have output contigs in `<ouput_folder>/chr<chrnum>/assembly/final_contigs/final_contigs.fa` (for CLR and ONT data) or `<ouput_folder>/chr<chrnum>/assembly/final_contigs/final_contig.p_ctg.fa` (for Hifi data).


### VolcanoSV Variant Call: 

#### Large Indel detection

The main code is `bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py`. The input arguments for this code are explained below:


```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --data_type DATA_TYPE, -dtype DATA_TYPE
                        CCS;CLR;ONT
  --rbam_file RBAM_FILE, -rbam RBAM_FILE
                        reads bam file for reads signature extraction
  --reference REFERENCE, -ref REFERENCE
                        wgs reference file
  --n_thread N_THREAD, -t N_THREAD
  --n_thread_align N_THREAD_ALIGN, -ta N_THREAD_ALIGN
  --mem_per_thread MEM_PER_THREAD, -mempt MEM_PER_THREAD
                        Set maximum memory per thread for alignment; suffix K/M/G recognized; default = 768M

```
The input directory should be the output directory of VolcaoSV-asm.
After running the above code, you will have output VCF in `<ouput_folder>/variants_filter_DEL_corGT.vcf`.


#### Complex SV detection

The main code is `bin/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py`. The input arguments for this code are explained below:


```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --indelvcf INDELVCF, -vcf INDELVCF
  --bamfile BAMFILE, -bam BAMFILE
  --reference REFERENCE, -ref REFERENCE
  --datatype {CCS,CLR,ONT}, -d {CCS,CLR,ONT}
  --out_dir OUT_DIR, -o OUT_DIR
  --n_thread N_THREAD, -t N_THREAD

```
The input directory should be the output directory of VolcaoSV-asm.
After running the above code, you will have output VCF in `<ouput_folder>/complex_SV.vcf`.


#### Small Indel detection

The main code is `bin/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py`. The input arguments for this code are explained below:


```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --read_bam READ_BAM, -rbam READ_BAM
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --reference REFERENCE, -ref REFERENCE
  --bedfile BEDFILE, -bed BEDFILE
                        optional; a high confidence bed file (default: None)
  --region REGION, -r REGION
                        optional; exmaple: chr21:2000000-2100000 (default: None)
  --n_thread N_THREAD, -t N_THREAD
  --kmer_size KMER_SIZE, -k KMER_SIZE
  --ratio RATIO, -rt RATIO
                        maximum bad kmer ratio (default: 0.3)
  --min_support MIN_SUPPORT, -ms MIN_SUPPORT
                        maximum support for bad kmer (default: 5)

```

The input directory should be the output directory of VolcaoSV-asm.
After running the above code, you will have output VCF in `<ouput_folder>/indel_2_49.vcf`.








#### Memory/Time Usage For VolcanoSV
Coverage| Memory| Time for one SV on a single node 
--- | --- | --- | 
60X | 20GB | 00:10:32 |





## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziezhoulab/VolcanoSV/issues">VolcanoSV</a>. 


