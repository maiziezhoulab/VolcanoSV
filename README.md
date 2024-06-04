
# :milky_way: VolcanoSV 🌋
## Table of Content 
- [Installation](#install-through-github)
- [Single chromosome mode (test example included)](#single-chromosome-mode)
    - [Single chromosome mode VolcanoSV Assembly](#Single-chromosome-mode-VolcanoSV-Assembly)
    - [Single chromosome mode Large Indel detection](#Single-chromosome-mode-Large-Indel-detection)
    - [Single chromosome mode Small Indel detection](#Single-chromosome-mode-Small-Indel-detection)
- [WGS mode](#wgs-mode)
  -  [WGS mode VolcanoSV Assembly](#wgs-mode-VolcanoSV-Assembly)
  -  [WGS mode Large Indel detection](#wgs-mode-Large-Indel-detection)
  -  [WGS mode Complex SV detection](#wgs-mode-complex-sv-detection)
  -  [WGS mode Small Indel detection](#wgs-mode-small-Indel-detection)
- [Truvari evaluation](#Truvari-evaluation)
- [Computation resource usage](#Computation-resource-usage)
- [Troubleshooting](#Troubleshooting)
  
# Install through Github:

```
git clone https://github.com/maiziezhoulab/VolcanoSV.git
```

## Dependencies for Github installation:
VolcanoSV utilizes Python3.8.3. To set up the environment, you need to have conda installed. Then, simply run
```
conda env create -f VolcanoSV/requirement.yaml
```

Then you will have a virtual environment called "volcano” created. Before running any VolcanoSV commands, please activate this environment first.
```
conda activate volcano
```

<!-- Add VolcanoSV to `PATH` by
```
path_to_volcano=/path/to/VolcanoSV/bin
PATH+=":$path_to_volcano/VolcanoSV-asm"
for i in $path_to_volcano/VolcanoSV-vc/*/; do PATH+=":$i"; done
chmod +x $path_to_volcano/VolcanoSV-asm/volcanosv-asm.py $path_to_volcano/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py $path_to_volcano/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py $path_to_volcano/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py
```
You can also add the line above to `~/.bashrc` if you don't want to do it everytime you start a shell -->
<!-- Next, put the "VolcanoSV/bin" in the ".bashrc" file, and source the ".bashrc" file <br /> -->

You can set
```
path_to_volcano=/path/to/VolcanoSV/bin
```
for convenience or just use the full path of "**volcanosv-asm.py**", "**volcanosv-vc-large-indel.py**", "**volcanosv-vc-complex-sv.py**" and "**volcanosv-vc-small-indel.py**"


# Running The Code:

## Single chromosome mode

For the single chromosome mode, we provided the chr10 BAM file, contigs file and VCF file for Hifi, CLR and ONT data. You can download them from [zenodo](https://zenodo.org/records/10520476).
In the following sessions, we will provide the code to run the Hifi data. If you wish to reproduce the result for CLR data or ONT data, you can just simply change the input BAM file and the argument "dtype" to the corresponding data type (CLR/ONT).

The example data is aligned to hg19 reference. You can download the reference files(genome.fa and genome.fa.fai) from zenode(https://zenodo.org/records/10520476).

Alternatively, you can download the more complete reference data (including more indexing and meta data) using the command below
```
wget https://cf.10xgenomics.com/supp/genome/refdata-hg19-2.1.0.tar.gz
tar -xzvf refdata-hg19-2.1.0.tar.gz
```

Note: since translocation detection requires WGS BAM file as support, it does not make sense to run it on single chromsome level. Therefore, we only provide the complex SV pipeline in WGS mode.

### Single chromosome mode VolcanoSV Assembly (VolcanoSV-asm) 


The VolcanoSV assembly pipeline is designed to run by chromosomes. The main code is `$path_to_volcano/VolcanoSV-asm/volcanosv-asm.py`. The input arguments for this code are explained below:

```
  --inbam INBAM, -i INBAM, could be either wgs bam or single-chromosome bam file
  --out_dir OUT_DIR, -o OUT_DIR
  --reference REFERENCE, -r REFERENCE
  --n_thread N_THREAD, -t N_THREAD
  --chrnum {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}, -chr {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
  --dtype {CCS,CLR,ONT}, -d {CCS,CLR,ONT}
  --prefix PREFIX, -px PREFIX

```

After running the above code, you will have output contigs in `<ouput_folder>/chr<chrnum>/assembly/final_contigs/final_contigs.fa`.

For example, you can reproduce contigs file for Hifi data using the below scripts

```
python3 $path_to_volcano/VolcanoSV-asm/volcanosv-asm.py \
-i Hifi_L2_hg19_minimap2_chr10.bam \
-o volcanosv_asm_output \
-r refdata-hg19-2.1.0/fasta/genome.fa \
-t 10 \
-chr 10 \
-d CCS \
-px Hifi_L2
```
The final contig will be `volcanosv_asm_output/chr10/assembly/final_contigs/final_contigs.fa`. 
If the volcanosv-asm pipeline is executed successfully and completely, your final contig file should have roughly the same size as the Hifi_L2_contigs.fa from zenodo.



### Single chromosome mode Large Indel detection (VolcanoSV-vc) 

The main code is `$path_to_volcano/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py`. The input arguments for this code are explained below:

```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --data_type DATA_TYPE, -dtype DATA_TYPE
                        CCS;CLR;ONT
  --rbam_file RBAM_FILE, -rbam RBAM_FILE
                        reads bam file for reads signature extraction
  --reference REFERENCE, -ref REFERENCE
                        wgs reference file
  --chrnum {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}, -chr {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
  --n_thread N_THREAD, -t N_THREAD
  --n_thread_align N_THREAD_ALIGN, -ta N_THREAD_ALIGN
  --mem_per_thread MEM_PER_THREAD, -mempt MEM_PER_THREAD
                        Set maximum memory per thread for alignment; suffix K/M/G recognized; default = 768M

```
The input directory should be the output directory of VolcaoSV-asm. This code is compatible with either single chromosome mode or wgs mode: when the argument "chrnum" is provided, it will execute in single chromosome mode, otherwise it will assume the input_dir contains chr1-chr22 contigs and execute in wgs mode.
After running the above code, you will have output VCF in `<ouput_folder>/volcanosv_large_indel.vcf`.

For example, if you want to reproduce the large indel VCF file for Hifi_L2 data, you can use the following command:
```
python3 $path_to_volcano/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_large_indel_output/ \
-dtype CCS \
-rbam Hifi_L2_hg19_minimap2_chr10.bam \
-ref refdata-hg19-2.1.0/fasta/genome.fa \
-chr 10 -t 10
```
The VCF file will be `volcanosv_large_indel_output/volcanosv_large_indel.vcf`. 
If the volcanosv-vc-large-indel pipeline is executed successfully and completely, your VCF file should have roughly the same number of variants as the Hifi_L2_variants.vcf from zenodo.
Note that, due to the randomness in assembly and alignment procedure, your VCF file may have 1 or 2 variants more or less than the Hifi_L2_variants.vcf. If that happens, we may still consider the pipeline as executed successfully, as long as the difference is minor.


### Single chromosome mode Small Indel detection (VolcanoSV-vc) 

The main code is `$path_to_volcano/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py`. The input arguments for this code are explained below:


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

The example code is as below:
```
python3 $path_to_volcano/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_small_indel \
-rbam Hifi_L2_hg19_minimap2_chr10.bam \
-ref refdata-hg19-2.1.0/fasta/genome.fa \
-r chr10 \
-t 30
```
After running the above code, you will have output VCF in `volcanosv_small_indel/volcanosv_small_indel.vcf`.

## WGS mode

### WGS mode VolcanoSV Assembly (VolcanoSV-asm) 

The VolcanoSV assembly is designed to run by chromosomes. If you have a distributed computing system that allows you to submit multiple jobs, we recommend that you submit one job per chromosome and let them run simultaneously. You can follow the template below to construct your job script:
```
python3 $path_to_volcano/VolcanoSV-asm/volcanosv-asm.py \
-i <wgs_bam> \
-o volcanosv_asm_output \
-r <reference_file> \
-t 10 \
-chr <chromosome_number> \
-d <datatype> \
-px <prefix>
```
For different chromosomes, you just need to simply change the `<chromosome_number>` to make different job scripts. All scripts share the same output folder. After all jobs are finished, you will see 22 subfolders under the output folder you specified, which are the chr1-chr22 assembly results. This is the most efficient way to run volcanosv-asm pipeline on a distributed system.


However, if you do not have a distributed computing system, you may write a for loop to run different chromosomes:
```
for i in {1..22}
do
echo "***********************assembly for chr${i}*******************"
python3 $path_to_volcano/VolcanoSV-asm/volcanosv-asm.py \
-i <wgs_bam> \
-o volcanosv_asm_output \
-r <reference_file> \
-t 10 \
-chr ${i} \
-d <datatype> \
-px <prefix>
done
```
The chr1-chr22 will be saved under `volcanosv_asm_output`. This is slower but can still get the work done.






### WGS mode Large Indel detection (VolcanoSV-vc) 

The main code is `$path_to_volcano/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py`. This code is designed for both single chromosome mode and wgs mode. To run it in wgs mode, you should first finish running the volcanosv-asm pipeline and provide the assembly output folder as the input folder for this code. You should **not** provide the `chr' argument in WGS mode.

An example command is as below:
```
python3 $path_to_volcano/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_large_indel_output/ \
-dtype <dtype> \
-rbam <wgs_reads_bamfile> \
-ref <reference> \
-t 11
```
After running the above code, you will have output VCF in `volcanosv_large_indel_output/volcanosv_large_indel.vcf`.


### WGS mode Complex SV detection (VolcanoSV-vc) 

The main code is `$path_to_volcano/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py`. The input arguments for this code are explained below:


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

The example code is as below:
```
python3 $path_to_volcano/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py \
-i volcanosv_asm_output/ \
-vcf volcanosv_large_indel_output/volcanosv_large_indel.vcf \
-o volcanosv_complex_sv \
-dtype <dtype> \
-bam <wgs_reads_bamfile> \
-ref <reference> \
-t 11
```
After running the above code, you will have output VCF in `volcanosv_complex_sv/volcanosv_complex_SV.vcf`.


### WGS mode Small Indel detection (VolcanoSV-vc) 

The main code is `$path_to_volcano/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py`. The input arguments for this code are explained below:


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

The example code is as below:
```
python3 $path_to_volcano/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_small_indel \
-rbam <wgs_reads_bamfile> \
-ref <reference> \
-t 30
```

After running the above code, you will have output VCF in `volcanosv_small_indel/volcanosv_small_indel.vcf`.



## Truvari evaluation

We use truvari4.0.0 to perform benchmarking against genome in a bottle gold standard set in high confidence region. The parameter we use is 
```
p=0.5 P=0.5 r=500 S=30 O=0.01
```



## Computation resource usage

### Memory/Time Usage For VolcanoSV-asm
Data Type| Coverage| Memory| ncpu | Run time | CPU hours |
--- | --- | --- | ---|--- | --- | 
CCS | 56X | 168GB | 30| 2-01:07:59| 1474 |
CLR | 89x | 259GB | 20 | 5-09:05:59 | 2582|
ONT | 48x | 85GB| 30 | 6-02:33:59| 4397|


### Memory/Time Usage For VolcanoSV-VC-large-indel
Data Type| Coverage| Memory| ncpu | Run time | CPU hours |
--- | --- | --- | ---|--- | --- | 
CCS | 56X | 21GB | 50| 00:25:11| 21 |
CLR | 89x | 215GB | 50 | 02:27:35 | 123|
ONT | 48x | 34GB| 50 | 00:38:24| 32|


## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziezhoulab/VolcanoSV/issues">VolcanoSV</a>. 


