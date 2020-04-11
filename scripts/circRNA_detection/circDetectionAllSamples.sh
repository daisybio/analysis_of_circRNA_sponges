#!/bin/bash
dataset=$1
reference_prefix=$2
output_dir=$3
scripts_dir=$4
paired_end=$5

mkdir -p $output_dir
while read line
do
sample=$(cut -f1 <<< $line)
totalRNAfastq=$(cut -f2 <<< $line)
if [ "$paired_end" = true ] ; then
    totalRNAfastq2=$(cut -f3 <<< $line)
    bash ${scripts_dir}/circRNA_detection/circDetectionForPairedEndSample.sh $sample $totalRNAfastq $totalRNAfastq2 $reference_prefix $output_dir
else
    bash ${scripts_dir}/circRNA_detection/circDetectionForSample.sh $sample $totalRNAfastq $reference_prefix $output_dir
fi
echo $sample
done < $dataset

