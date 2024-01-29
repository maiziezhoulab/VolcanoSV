bamfile=$1
vcffile=$2
dtype=$3 #ONT HiFi CLR
sigdir=$4 # cutesv signature dir
#ft_vtype=$4 # INS/DEL/INSDEL


ft_vtype=DEL

filtered_vcf=`dirname $vcffile`/`basename $vcffile .vcf`_filter_${ft_vtype}.vcf
dtype_lc=`echo "$dtype" | tr '[:upper:]' '[:lower:]'`
# sigdir=`dirname $vcffile`/cutesv_signature/
gtdir=`dirname $vcffile`/GT_Correction
final_vcf=`dirname $vcffile`/`basename $vcffile .vcf`_filter_${ft_vtype}_corGT.vcf
workdir=`dirname $vcffile`

script_dir=`dirname $0`
################################# filter out false positive #####################################


# # cuteSv extract reads signature
# if [ -d ${sigdir} ]; then
#     echo "${sigdir} exists. Deleting it and creating a new one..."
#     rm -r ${sigdir}
# fi

# mkdir -p ${sigdir}

# if [ "$dtype" = HiFi ]; then
# python3 /data/maiziezhou_lab/CanLuo/Software/cuteSV/src/cuteSV/cuteSV \
#     $bamfile \
#     /data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa \
#     ./test.vcf \
#     ${sigdir} \
#     --max_cluster_bias_INS      1000 \
#     --diff_ratio_merging_INS    0.9 \
#     --max_cluster_bias_DEL  1000 \
#     --diff_ratio_merging_DEL    0.5 \
#     --retain_work_dir -t 30
# elif [ "$dtype" = CLR ]; then
# python3 /data/maiziezhou_lab/CanLuo/Software/cuteSV/src/cuteSV/cuteSV \
#     $bamfile \
#     /data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa \
#     ./test.vcf \
#     ${sigdir} \
#     --max_cluster_bias_INS      100 \
#     --diff_ratio_merging_INS    0.3 \
#     --max_cluster_bias_DEL  200 \
#     --diff_ratio_merging_DEL    0.5 \
#     --retain_work_dir -t 30
# elif [ "$dtype" = ONT ]; then
# python3 /data/maiziezhou_lab/CanLuo/Software/cuteSV/src/cuteSV/cuteSV \
#     $bamfile \
#     /data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa \
#     ./test.vcf \
#     ${sigdir} \
#     --max_cluster_bias_INS      100 \
#     --diff_ratio_merging_INS    0.3 \
#     --max_cluster_bias_DEL  100 \
#     --diff_ratio_merging_DEL    0.3 \
#     --retain_work_dir -t 30
# else
#     echo "dtype can only be HiFi, CLR or ONT"
#     exit 1
# fi


# calculate signature support for each variant based on cuteSV signature
echo "------------------Calculate signature support for each variant"
python3 ${script_dir}/calculate_signature_support.py  \
    -v $vcffile \
    -ct ${sigdir} \

# Use empirical threshold to filter SV
echo "------------------Use empirical threshold to filter SV"
python3 ${script_dir}/filter_vcf_by_sig_cov_insdel.py \
    -i $vcffile \
    -d ${dtype_lc} \
    -a volcano \
    -v ${ft_vtype}




############################################## GT correction ###################################
echo "------------------Correct DEL GT"
mkdir -p ${gtdir}
python3 ${script_dir}/correct_gt_del_real_data.py \
    -i ${filtered_vcf}  \
    -o $gtdir/bnd_del_real.tsv \
    -bam $bamfile \
    -sig $sigdir/DEL.sigs \
    -t 7 \
    -d $dtype -v DEL
    
echo "------------------Correct INS GT"
python3 ${script_dir}/correct_gt_ins_real_data.py \
    -i ${filtered_vcf}  \
    -o ${gtdir}/bnd_ins_real.tsv \
    -bam $bamfile \
    -sig $sigdir/INS.sigs \
    -t 7 \
    -d $dtype -v INS

grep '#' ${filtered_vcf} > ${workdir}/header
cat ${workdir}/header ${filtered_vcf}.newgt.DEL ${filtered_vcf}.newgt.INS | vcf-sort > ${final_vcf}
