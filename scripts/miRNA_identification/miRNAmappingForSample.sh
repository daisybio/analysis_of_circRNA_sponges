#!/bin/bash
sample=$1
fastq=$2
adapter=$3
species=$4
ref_dir=$5
ref_prefix=$6
out_dir=$7

# prepare output folder
echo "$(date) Started miRNA identification for $sample"	
output="${out_dir}/${sample}/miRNA_identification/"
mkdir -p $output
cd $output
reads_collapsed="${output}/reads_collapsed.fa"
reads_vs_ref="${output}/reads_vs_ref.arf"

# preprocessing, adapter clipping and mapping
echo "$(date) Preprocessing and mapping for $sample"
mapper.pl $fastq -e -h -i -j -k $adapter -l 18 -m -p $ref_prefix -s $reads_collapsed -t $reads_vs_ref -v -o 4 > "${output}/preprocessing.log"
echo "$(date) Preprocessing done for sample $sample"

# miRNA identification
echo "$(date) Started miRNA identification for $sample"
miRDeep2.pl $reads_collapsed "${ref_prefix}.fa" $reads_vs_ref "${ref_dir}/mature_ref.fa" "${ref_dir}/mature_other.fa" "${ref_dir}hairpin_ref.fa" -t $species -d -v 2>"${output}/report.log"

echo "$(date) Finished miRNA identification for $sample"
echo "$(date) DONE $sample"