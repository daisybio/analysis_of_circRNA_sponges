# Pipeline for the analysis of circRNA sponges
This pipeline covers the systematical detection and analysis of circRNAs that sponge miRNAs from total RNA and small RNA sequencing data.

## Overview
![pipeline_overview](https://user-images.githubusercontent.com/51077615/70857669-3b5c1c80-1ef3-11ea-9a12-3d8427781cf6.PNG)



## Dependencies

Thie pipeline is written in bash and R and requires several external tools:
* Bowtie version 1.0.0
* Bowtie2 version 2.1.0
* TopHat2 version 2.0.9
* CIRCExplorer2 version 2.3.6
* miRanda version 3.3a
* miRDeep2 version 2.0.1.2
* samtools
* R packages: ggplot2, plyr, dplyr, data.table, gridExtra, grid, ggpubr

This pipeline is also available as a Docker container.

## Input

##### Requirements
For this analysis, a minimum of 5 samples coming from the same organism are required. For each sample, two fastq are needed: one obtained through total RNA-seq (or rRNA-depleted) and one obtained through small RNA-seq (miRNA). 
Various reference files are also needed: a GTF file, a genome fasta file, miR-Base reference files for hairpin and mature miRNA sequences. The exact input format and a way for obtaining the needed reference files are described below.

### References preparation
##### Get genome fasta file, gtf file. 
This example is for mouse annotation mm10. If you use another organism or annotation, please replace it in the following commands.
For this step you will need genePredToGtf from UCSC Utilities (http://hgdownload.soe.ucsc.edu/admin/exe/).
For further information, please refer to: https://circexplorer2.readthedocs.io/en/latest/tutorial/setup/#installation

```bash
fetch_ucsc.py mm10 fa mm10.fa
fetch_ucsc.py mm10 ref mm10_ref.txt
fetch_ucsc.py mm10 kg mm10_kg.txt
cut -f2-11 mm10_kg.txt|./data/reference/genePredToGtf file stdin mm10_kg.gtf
```

##### Get miRNA reference files from mirbase. 
For further information, please refer to: https://drmirdeep.github.io/mirdeep2_tutorial.html 
Get mirbase file:
``` bash
git clone https://github.com/rajewsky-lab/mirdeep2.git
git clone https://github.com/Drmirdeep/mirdeep2_patch.git
cd mirdeep2
perl install.pl
miRDeep2.pl 
cd mirdeep2_patch
bash patchme.sh

perl mirbase.pl 22
perl mirbase.pl 22 1
```

Extract the mature sequences from the mirbase file downloaded before and get the hairpin sequences. This example applies to mouse (mmu). For other organism, pleas use the repsective three letter code  (for example hsa for human).  
Also extract mature miRNA information from related species(here rat(rno) and human (hsa)).
``` bash
extract_miRNAs.pl ~/mirbase/22/mature.fa.gz mmu > mature_ref.fa
extract_miRNAs.pl ~/mirbase/22/hairpin.fa.gz mmu > hairpin_ref.fa
extract_miRNAs.pl ~/mirbase/22/mature.fa.gz rno,hsa > mature_other.fa
```

##### Build Bowtie index from genome fasta file
For this step, Bowtie1 and Bowtie2 are needed (see version in section Dependencies). 
```bash
bowtie-build mm10.fa mm10
bowtie2-build mm10.fa mm10
```

### Input data preparation

#### Folder structure
```bash
├─── input_folder
│   ├─── parameters.txt
│   ├─── dataset.tsv
│   └─── data
│       ├─── miRNA_fastq
│       |   ├─── <sample1>.fastq
│       |   ├─── <sample2>.fastq
│       |   └─── ...
│       └─── circRNA_fastq
│           ├─── <sample1>.fastq
│           ├─── <sample2>.fastq
│           └─── ...
|
├─── references
│   ├─── <mm10>.fa
│   ├─── <mm10>_ref.txt
│   ├─── <mm10>_kg.gtf
│   ├─── hairpin_ref.fa
│   ├─── mature_other.fa
│   ├─── mature_ref.fa
│   ├─── <mm10>.*.ebwt
│   ├─── <mm10>.rev.*.ebwt
│   ├─── <mm10>.*.bt2
│   ├─── <mm10>.rev.*.bt2
|
└───  output_folder (empty)
```

#### Parameters
Replace every parameter in the file ```input/parameters.txt``` with information suitable for your dataset.
```bash
dataset=input_folder/dataset.tsv
adapter="TGGAATTCTCGGGTGCCAAGG"
species=mmu
ref_dir=references/
ref_prefix=references/mm10
out_dir=output_folder
scripts_dir=/bin/scripts/
```

##### dataset
The dataset file contains the name of the sample (first column), the path to the circRNA fastq file corresponding to the mentioned sample (second column), the path to the miRNA fastq file corresponding to the mentioned sample (third column). The file should be tab-separated and contain no header. The dataset file should look like this:
```
cerebellum_rep1 input_folder/data/circRNA_fastq/<sample1>.fastq	input_folder/data/miRNA_fastq/<sample1>.fastq
cerebellum_rep2	input_folder/data/circRNA_fastq/<sample2>.fastq	input_folder/data/miRNA_fastq/<sample1>.fastq
cerebellum_rep3	input_folder/data/circRNA_fastq/<sample3>.fastq	input_folder/data/miRNA_fastq/<sample1>.fastq
hippocampus_rep1  input_folder/data/circRNA_fastq/<sample4>.fastq	input_folder/data/miRNA_fastq/<sample1>.fastq
...
```

##### adapter
This depends on the sequencing library you used. In order to find out which library adapter was used, check: https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
For example TruSeq Small RNA uses the adapter "TGGAATTCTCGGGTGCCAAGG".

##### species
Three letter code for the species used in this experiment. For example mmu (mouse), hsa(human).

##### ref_dir and ref_prefix
ref_dir is the directory containing all references. ref_prefix is the reference directory followed by the prefix used for building the bowtie index.

##### scripts_dir
```scripts_dir=/bin/scripts``` for the use of Docker Container. If you prefer to run the pipeline manually, please specify the location of the scripts.


### Build Docker
Download the Dockerfile and the folders scripts and soft. In the same directory run:
```docker build -t pipeline_docker .```
