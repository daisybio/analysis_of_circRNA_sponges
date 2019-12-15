#!/bin/bash
dataset=$1
adapter=$2
species=$3
ref_dir=$4
ref_prefix=$5
out_dir=$6
scripts_dir=$7

mkdir -p $out_dir

while read line
do
	sample=$(cut -f1 <<< $line)
	miRNAfastq=$(cut -f3 <<< $line)
	bash ${scripts_dir}/miRNAidentification/miRNAmappingForSample.sh $sample $miRNAfastq $adapter $species $ref_dir $ref_prefix $out_dir
done < $dataset
