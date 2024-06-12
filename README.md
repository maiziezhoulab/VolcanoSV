
# :milky_way: VolcanoSV 🌋
## Table of Content 
- [Installation](#install-through-github)
- [Single chromosome mode (test example included)](#single-chromosome-mode)
    - [Single chromosome mode VolcanoSV Assembly](#Single-chromosome-mode-VolcanoSV-Assembly-volcanosv-asm)
    - [Single chromosome mode VolcanoSV Assembly-Single assembler](#Single-chromosome-mode-VolcanoSV-Assembly-Single-assembler)
    - [Single chromosome mode VolcanoSV Assembly-hybrid](#Single-chromosome-mode-VolcanoSV-Assembly-hybrid)
    - [Single chromosome mode Large Indel detection](#Single-chromosome-mode-Large-Indel-detection-volcanosv-vc)
    - [Single chromosome mode Small Indel detection](#Single-chromosome-mode-Small-Indel-detection-volcanosv-vc)
- [WGS mode](#wgs-mode)
  -  [WGS mode VolcanoSV Assembly](#wgs-mode-VolcanoSV-Assembly-volcanosv-asm)
  -  [WGS mode Large Indel detection](#wgs-mode-Large-Indel-detection-volcanosv-vc)
  -  [WGS mode Complex SV detection](#wgs-mode-complex-sv-detection-volcanosv-vc)
  -  [WGS mode Small Indel detection](#wgs-mode-small-Indel-detection-volcanosv-vc)
- [Improve assembly for regions enriched in segmental duplications](#sd-recovery)
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

Then you will have a virtual environment called `volcanosv` created. Before running any VolcanoSV commands, please activate this environment first.
```
conda activate volcanosv
```

<!-- Add VolcanoSV to `PATH` by
```
path_to_volcanosv=/path/to/VolcanoSV
PATH+=":${path_to_volcanosv}/bin/VolcanoSV-asm"
for i in ${path_to_volcanosv}/bin/VolcanoSV-vc/*/; do PATH+=":$i"; done
chmod +x ${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm.py ${path_to_volcanosv}/bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py ${path_to_volcanosv}/bin/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py ${path_to_volcanosv}/bin/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py
```
You can also add the line above to `~/.bashrc` if you don't want to do it everytime you start a shell -->
<!-- Next, put the "VolcanoSV/bin" in the ".bashrc" file, and source the ".bashrc" file <br /> -->

You can set
```
path_to_volcanosv=/path/to/VolcanoSV/bin
```
for convenience or just use the full path of `${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm.py`, `${path_to_volcanosv}/bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py`, `${path_to_volcanosv}/bin/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py` and `${path_to_volcanosv}/bin/VolcanoSV-vc/Small_Indel/volcanosv-vc-small-indel.py`.


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

#### Single chromosome mode VolcanoSV Assembly - Single assembler
The VolcanoSV assembly pipeline is designed to run by chromosomes. We integrated multiple state-of-the-art assemblers into the pipeline, including [hifiasm](https://github.com/chhylp123/hifiasm), [Flye](https://github.com/fenderglass/Flye), [wtdbg2](https://github.com/ruanjue/wtdbg2),[miniasm](https://github.com/lh3/miniasm),[Shasta](https://github.com/paoloshasta/shasta),[NextDenovo](https://github.com/Nextomics/NextDenovo),and [Hicanu/Canu](https://github.com/marbl/canu). Users can select the appropriate assembler based on the needs. The main script is `${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm.py`. The input arguments for this script are explained below:

```
  --bam_file INBAM, -bam INBAM, could be either wgs bam or single-chromosome bam file
  --output_dir output_dir, -o output_dir
  --reference REFERENCE, -ref REFERENCE
  --n_thread N_THREAD, -t N_THREAD
  --chrnum {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}, -chr {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
  --assembler {wtdbg2,canu,miniasm,shasta,nextdenovo,hifiasm,hicanu,flye}, -asm {wtdbg2,canu,miniasm,shasta,nextdenovo,hifiasm,hicanu,flye}
        optional; if not set, VolcanoSV use hifiasm for Hifi data and flye for CLR and ONT data by default.
  --data_type {CLR,ONT,Hifi}, -dtype {CLR,ONT,Hifi}
  --pacbio_subtype {CLR-rs,CLR-sq}, -pb {CLR-rs,CLR-sq}
                        must provide when using wtdbg2 on CLR data (default: None)
  --shasta_ont_config {Nanopore-OldGuppy-Sep2020}, -shacon {Nanopore-OldGuppy-Sep2020}
  --prefix PREFIX, -px PREFIX

```
Please select from hifiasm and hicanu for Hifi data, and the rest of the assemblers are for CLR and ONT data.
By default, VolcanoSV uses hifiasm for Hifi data and flye for CLR and ONT data.
After running the above code, you will have output contigs in `<ouput_folder>/chr<chrnum>/assembly/final_contigs/<prefix>_final_contigs.fa`.

For example, if you want to use hifiasm for hifi data, you can use the below scripts

```
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm.py \
-bam Hifi_L2_hg19_minimap2_chr10.bam \
-o volcanosv_asm_output \
-ref refdata-hg19-2.1.0/fasta/genome.fa \
-t 10 \
-chr 10 \
-dtype Hifi \
-px Hifi_L2 \
-asm hifiasm
```
The final contig will be `volcanosv_asm_output/chr10/assembly/final_contigs/Hifi_L2_final_contigs.fa`. 
If the volcanosv-asm pipeline is executed successfully, your final contig file should have roughly the same size as the Hifi_L2_contigs.fa from zenodo.
VolcanoSV-asm already includes the executable version of all assemblers, so you do not need to install them individually.
However, if you want more detailed information on these assemblers, you can [click here](Assemblers.md).
#### Single chromosome mode VolcanoSV Assembly - Hybrid mode
Different assemblers vary in their ability to assemble regions enriched in segmental duplications (SDs) and other complex regions. Therefore, it is often advantageous to utilize different assemblers for different genomic regions. We thus also provide a hybrid mode: users can input a BED file, and specify an "in-BED" assembler and an "out-BED" assembler. The phase blocks that overlap with the BED file will be assembled using the in-BED assembler, while the rest will be assembled by the out-BED assembler. The script for this mode is `${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm_hybrid.py`.

For example, if you provide a `segdups.bed`, and want to use hicanu for the segdup regions and hifiasm for the other rest regions, you can use the code below:

```
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm_hybrid.py \
-bam Hifi_L2_hg19_minimap2_chr10.bam \
-o volcanosv_asm_output \
-ref refdata-hg19-2.1.0/fasta/genome.fa \
-bed segdups.bed \
-inasm hicanu \
-outasm hifiasm \
-t 10 \
-chr 10 \
-dtype Hifi \
-px Hifi_L2 
```
The final contig will be `volcanosv_asm_output/chr10/assembly/final_contigs/Hifi_L2_final_contigs.fa`. 
If the volcanosv-asm pipeline is executed successfully, your final contig file should have roughly the same size as the Hifi_L2_contigs.fa from zenodo.


### Single chromosome mode Large Indel detection (VolcanoSV-vc) 

The main code is `${path_to_volcanosv}/bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py`. The input arguments for this code are explained below:

```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --data_type DATA_TYPE, -dtype DATA_TYPE
                        Hifi;CLR;ONT
  --bam_file RBAM_FILE, -bam RBAM_FILE
                        reads bam file for reads signature extraction
  --reference REFERENCE, -ref REFERENCE
                        wgs reference file
  --chrnum {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}, -chr {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
  --n_thread N_THREAD, -t N_THREAD
  --n_thread_align N_THREAD_ALIGN, -ta N_THREAD_ALIGN
  --mem_per_thread MEM_PER_THREAD, -mempt MEM_PER_THREAD
                        Set maximum memory per thread for alignment; suffix K/M/G recognized; default = 768M
  --prefix PREFIX, -px PREFIX

```
The input directory should be the output directory of volcanoSV-asm. This code is compatible with either single chromosome mode or wgs mode: when the argument "chrnum" is provided, it will execute in single chromosome mode, otherwise, it will assume the input_dir contains chr1-chr22 contigs and execute in wgs mode. Please note that `prefix` should remain consistent with what is set in volcanosv-asm.
After running the above code, you will have output VCF in `<ouput_folder>/volcanosv_large_indel.vcf`.

For example, if you want to reproduce the VCF file for large indels on Hifi_L2 data, you can use the following command:
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_large_indel_output/ \
-dtype Hifi \
-bam Hifi_L2_hg19_minimap2_chr10.bam \
-ref refdata-hg19-2.1.0/fasta/genome.fa \
-chr 10 -t 10 \
-px Hifi_L2
```
The VCF file will be `volcanosv_large_indel_output/Hifi_L2_volcanosv_large_indel.vcf`. 
If the volcanosv-vc-large-indel pipeline is executed successfully, your VCF file should have roughly the same number of variants as the Hifi_L2_variants.vcf from zenodo.
Note that, due to the randomness in assembly and alignment procedure, your VCF file may have 1 or 2 variants more or less than the Hifi_L2_variants.vcf. If that happens, we may still consider the pipeline as executed successfully, as long as the difference is minor.


### Single chromosome mode Small Indel detection (VolcanoSV-vc) 

The main script is `${path_to_volcanosv}/bin/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py`. The input arguments for this script are explained below:


```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --bam_file READ_BAM, -bam READ_BAM
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
  --prefix PREFIX, -px PREFIX


```

The input directory should be the output directory of volcanoSV-asm.

The example code is as below:
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_small_indel \
-bam Hifi_L2_hg19_minimap2_chr10.bam \
-ref refdata-hg19-2.1.0/fasta/genome.fa \
-r chr10 \
-t 30 \
-px Hifi_L2
```
After running the above code, you will have output VCF in `volcanosv_small_indel/Hifi_L2_volcanosv_small_indel.vcf`.

## WGS mode

### WGS mode VolcanoSV Assembly (VolcanoSV-asm) 

The VolcanoSV assembly is designed to operate on a per-chromosome basis. If you have access to a distributed computing system that supports multiple job submissions, we recommend submitting one job per chromosome and running them concurrently. You can follow the template below to construct your job script:
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm.py \
-bam <wgs_bam> \
-o volcanosv_asm_output \
-ref <reference_file> \
-t 10 \
-chr <chromosome_number> \
-dtype <datatype> \
-px <prefix>
```
To create job scripts for different chromosomes, simply change the <chromosome_number> in each script. All scripts share the same output folder. After all jobs are finished, you will see 22 subfolders under the output folder you specified, which are the chr1-chr22 assembly results. This is the most efficient way to run volcanosv-asm pipeline on a distributed system.


However, if you do not have a distributed computing system, you may write a for loop to run different chromosomes:
```
for i in {1..22}
do
echo "***********************assembly for chr${i}*******************"
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/volcanosv-asm.py \
-bam <wgs_bam> \
-o volcanosv_asm_output \
-ref <reference_file> \
-t 10 \
-chr ${i} \
-dtype <datatype> \
-px <prefix>
done
```
The chr1-chr22 will be saved under `volcanosv_asm_output`. This is slower but can still get the work done.






### WGS mode Large Indel detection (VolcanoSV-vc) 

The main script is `${path_to_volcanosv}/bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py`. This script is designed for both single chromosome mode and wgs mode. To run it in WGS mode, you should first finish running the volcanosv-asm pipeline and provide the assembly output folder as the input folder for this code. You should **not** provide the `chr' argument in WGS mode.

An example command is as below:
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-vc/Large_INDEL/volcanosv-vc-large-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_large_indel_output/ \
-dtype <datatype> \
-bam <wgs_reads_bamfile> \
-ref <reference> \
-t 11 \
-px <prefix>
```
After running the above code, you will have output VCF in `volcanosv_large_indel_output/<prefix>_volcanosv_large_indel.vcf`.


### WGS mode Complex SV detection (VolcanoSV-vc) 

The main code is `${path_to_volcanosv}/bin/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py`. The input arguments for this code are explained below:


```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --indelvcf INDELVCF, -vcf INDELVCF
  --bamfile BAMFILE, -bam BAMFILE
  --reference REFERENCE, -ref REFERENCE
  --datatype {Hifi,CLR,ONT}, -dtype {Hifi,CLR,ONT}
  --output_dir output_dir, -o output_dir
  --n_thread N_THREAD, -t N_THREAD
  --prefix PREFIX, -px PREFIX
```
The input directory should be the output directory of volcanoSV-asm.

The example code is as below:
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-vc/Complex_SV/volcanosv-vc-complex-sv.py \
-i volcanosv_asm_output/ \
-vcf volcanosv_large_indel_output/volcanosv_large_indel.vcf \
-o volcanosv_complex_sv \
-dtype <datatype> \
-bam <wgs_reads_bamfile> \
-ref <reference> \
-t 11 \
-px <prefix>
```
After running the above code, you will have output VCF in `volcanosv_complex_sv/<prefix>_volcanosv_complex_SV.vcf`.


### WGS mode Small Indel detection (VolcanoSV-vc) 

The main code is `${path_to_volcanosv}/bin/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py`. The input arguments for this code are explained below:


```
  --input_dir INPUT_DIR, -i INPUT_DIR
  --bam_file READ_BAM, -bam READ_BAM
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
  --prefix PREFIX, -px PREFIX


```

The input directory should be the output directory of volcanoSV-asm.

The example code is as below:
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-vc/Small_INDEL/volcanosv-vc-small-indel.py \
-i volcanosv_asm_output/ \
-o volcanosv_small_indel \
-bam <wgs_reads_bamfile> \
-ref <reference> \
-t 30 \
-px <prefix>
```

After running the above code, you will have output VCF in `volcanosv_small_indel/<prefix>_volcanosv_small_indel.vcf`.


## Improve assembly for regions enriched in segmental duplications (SDs)

After WGS assembly, if you would like to evaluate assembly for SDs and further achieve better assembly in SD-enriched regions, you can run the below pipeline, which includes 3 steps. 

### Step1
Align reads to the contig fasta file, and then utilize [Flagger](https://github.com/mobinasri/flagger) to annotate assembly for collapse components (collapsed SD regions).
To run this step, you need java and docker in your system.
```
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/Evaluate_Assembly.py \
  --input_dir <volcanosv_output> \
  --output_dir <SD_recovery_dir> \
  --fastq_file <raw_reads_fastq> \
  --data_type <DATA_TYPE> \
  --n_thread <t> \
  --mem_per_thread <mem> \
  --sample_name <sample> \
  --lib_name <lib>
```

After this step, you will have a `<sample>_<lib>_collapsed_hp_namex.txt` file generated in the output folder, which contains the haplotype names that contain collapsed SDs.

### Step2
Perform assembly only in those collapsed regions using a specified assembler.
The main script is `General_Assembly_Workflow.py`. The arguments are:
```
  --hap_file HAP_FILE, -haps HAP_FILE
  --fastq_dirs FASTQ_DIRS [FASTQ_DIRS ...], -fqds FASTQ_DIRS [FASTQ_DIRS ...]
                        the folder that includes FASTQ files. (default: None)
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --assemblers {wtdbg2,canu,miniasm,shasta,nextdenovo,hifiasm,hicanu,flye} [{wtdbg2,canu,miniasm,shasta,nextdenovo,hifiasm,hicanu,flye} ...], -asms {wtdbg2,canu,miniasm,shasta,nextdenovo,hifiasm,hicanu,flye} [{wtdbg2,canu,miniasm,shasta,nextdenovo,hifiasm,hicanu,flye} ...]
                        the assemblers used for fastq_dirs; order corresponds to the order of fastq_dirs (default: None)
  --data_type {CLR,ONT,Hifi}, -d {CLR,ONT,Hifi}
  --pacbio_subtype {rs,sq}, -pb {rs,sq}
                        must provide when using wtdbg2 on CLR data (default: None)
  --shasta_ont_config {Nanopore-OldGuppy-Sep2020}, -shacon {Nanopore-OldGuppy-Sep2020}
  --n_thread N_THREAD, -t N_THREAD
  --asm_thread ASM_THREAD, -ta ASM_THREAD
  --clean, -cl
  --prefix PREFIX, -px PREFIX
                        file prefix in the output folder (default: Sample)
```
Example Usage:

```
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/General_Assembly_Workflow.py \
-hap <sample>_<lib>_collapsed_hp_namex.txt \
-fqds <volcanosv_output>/chr*/Assembly/fastq_by_hap \
-o  <volcanosv_output>/SD_recovery
-asms <your_specified_assembler> \
-d <type> \
-t <t> 
```
You will have a `<volcanosv_output>/SD_recovery/final_contigs/final_contigs.fa` generated.

### Step3
Use the newly generated contigs to replace the previously collapsed contigs.

```
python3 ${path_to_volcanosv}/bin/VolcanoSV-asm/Replace_Collapsed_Contigs.py \
-og <volcanosv_output>/SD_recovery/assemblies.fa \
-new <volcanosv_output>/SD_recovery/final_contigs/final_contigs.fa \
-o <volcanosv_output>/SD_recovery/SD_recovered.fa \
-hap  <sample>_<lib>_collapsed_hp_namex.txt 
```
<volcanosv_output>/SD_recovery/SD_recovered.fa is the SD recovered contig file.

## Truvari evaluation

We use truvari4.0.0 to perform benchmarking against the Genome in a Bottle (GIAB) gold standard set in a high confidence region. The parameter we use is 
```
p=0.5 P=0.5 r=500 S=30 O=0.01
```



## Computation resource usage

### Memory/Time Usage For VolcanoSV-asm
Data Type| Coverage| Memory| ncpu | Run time | CPU hours |
--- | --- | --- | ---|--- | --- | 
Hifi | 56X | 168GB | 30| 2-01:07:59| 1474 |
CLR | 89x | 259GB | 20 | 5-09:05:59 | 2582|
ONT | 48x | 85GB| 30 | 6-02:33:59| 4397|


### Memory/Time Usage For VolcanoSV-VC-large-indel
Data Type| Coverage| Memory| ncpu | Run time | CPU hours |
--- | --- | --- | ---|--- | --- | 
Hifi | 56X | 21GB | 50| 00:25:11| 21 |
CLR | 89x | 215GB | 50 | 02:27:35 | 123|
ONT | 48x | 34GB| 50 | 00:38:24| 32|


## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziezhoulab/VolcanoSV/issues">VolcanoSV</a>. 


