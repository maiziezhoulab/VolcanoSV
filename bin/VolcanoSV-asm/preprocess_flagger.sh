reads=${1}
workdir=${2}
bam_prefix=${3}
dtype=${4} # hifi pb ont
t=${5} # num of threads
m=${6} # mem per threads



cd ${workdir}

temp_bam=${bam_prefix}_temp.bam
ref=assemblies.fa

mkdir -p temp

# ml  GCC/11.3.0  SAMtools/1.18

minimap2 -t ${t}  -Y -L -a -H --split-prefix temp -x  map-${dtype} ${ref} ${reads} | samtools index -@ ${t} -m ${m} -T temp/ -o  ${temp_bam}

samtools index -@ ${t} -m ${m}  ${temp_bam}

samtools calmd ${temp_bam} ./assemblies.fa -b -@ ${t} > ${bam_prefix}_MD.bam
rm ${temp_bam}*

# Compress the FASTA file
bgzip -c assemblies.fa > assemblies.fa.gz

# Index the compressed FASTA file
samtools faidx assemblies.fa.gz

# rm assemblies.fa