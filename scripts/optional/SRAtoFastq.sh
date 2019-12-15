#/bin/bash
sraFiles=$1    # directory containing all sra files
output_dir=$2  # output directory
sra_bin=$3     # bin directory of sratoolkit

for f in $sraFiles
do
	  echo "Processing $f file..."
	${sra_bin}/fastq-dump -I --split-files -O $output_dir $f
done
