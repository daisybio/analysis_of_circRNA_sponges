#!/bin/bash
sample=$1
fastq=$2
reference_prefix=$3
output_dir=$4

# prepare output directory
echo "Started processing $sample"
sample_out="$output_dir/$sample/"
mkdir "$sample_out"
sample_out="$sample_out/circRNA_detection"
mkdir "$sample_out"

# mapping with TopHat2
echo "$(date) $sample : Starting mapping for $sample"
tophat2 -a 6 --microexon-search -m 2 -p 16 -o "$sample_out/tophat2/" "$reference_prefix" $fastq

# process unmapped reads
echo "$(date) $sample : Processing unmapped reads for $sample"
bamToFastq -i "$sample_out/tophat2/unmapped.bam" -fq "$sample_out/tophat2/unmapped.fastq"

# map unmapped reads with TopHat Fusion
echo "$(date) $sample : Mapping unmapped reads for $sample"
tophat2 -o "$sample_out/tophat_fusion" -p 15 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search "$reference_prefix" "$sample_out/tophat2/unmapped.fastq"

mkdir "$sample_out/CIRCexplorer2"

# CIRCExplorer2 parsing module
echo "$(date) $sample Preparing circRNA detection for $sample"
CIRCexplorer2 parse -b "$sample_out/CIRCexplorer2/back_spliced_junction.bed" -t TopHat-Fusion "$sample_out/tophat_fusion/accepted_hits.bam" > "$sample_out/CIRCexplorer2/CIRCexplorer2_parse.log"

# CIRCExplorer2 annotation module
echo "$(date) $sample Detection of circRNAs for $sample"
CIRCexplorer2 annotate -r "${reference_prefix}_ref.txt" -g "${reference_prefix}.fa" -b "$sample_out/CIRCexplorer2/back_spliced_junction.bed" -o "$sample_out/CIRCexplorer2/circularRNA_known.txt" > "$sample_out/CIRCexplorer2/CIRCexplorer2_annotate.log"

echo "$(date) $sample Annotation done"
