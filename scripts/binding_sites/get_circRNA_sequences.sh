#!/bin/bash
ref_prefix=$1
output_dir=$2
input="${output_dir}/results/circRNA/circRNA_counts_all_samples_filtered.tsv"
fasta="${ref_prefix}.fa"
output="${output_dir}/results/binding_sites/input/circRNAs.fa"
while read line
do
	chr=$(echo $line | awk '{ print $1 }')
	from=$(echo $line | awk '{ print $2 }')
	to=$(echo $line | awk '{ print $3 }')
	strand=$(echo $line | awk '{ print $4 }')
	sam_out=$(samtools faidx $fasta "$chr:$from-$to")
	header=">${chr}:${from}-${to}_${strand}"
	echo $header >> $output
	if [ $strand == "+" ]; then
		sequence=$(echo "$sam_out" 2>&1 | tail -n +2 )
	elif [ $strand == "-" ]; then
		sequence=$(echo "$sam_out" 2>&1 | tail -n +2 | rev | tr {AGTCagtc} {TCAGtcag})

	fi		
		echo "$sequence" >> $output
done < <(tail -n +2 $input)


