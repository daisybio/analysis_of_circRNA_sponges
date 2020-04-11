# Pipeline for the analysis of circRNA sponges
This pipeline covers the systematical analysis of circRNAs that sponge miRNAs. It requires total RNA and small RNA sequencing data and is based on the hypothesis that the negatively correlating expression of miRNAs and circRNAs (having miRNA binding sites) is an indicator of sponging. This pipeline is available as a Docker container.

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
* samtools version 0.1.19
* R packages: ggplot2, plyr, dplyr, data.table, gridExtra, grid, ggpubr

## Input

For this analysis, a minimum of 5 samples coming from the same organism are required. For each sample, two fastq files are needed: one obtained through total RNA-seq (or rRNA-depleted) and one obtained through small RNA-seq. 
Various reference files are also needed: a GTF file, a genome fasta file, miR-Base reference files for hairpin and mature miRNA sequences. The exact input format and how to get the needed reference files is described below.

### References preparation
##### Get genome fasta file, gtf file. 
In this example, the mouse annotation mm10 is used. If you are using another organism or annotation, please adapt the following commands. For further information, please refer to: https://circexplorer2.readthedocs.io/en/latest/tutorial/setup/#installation
1. Download gene annotation file:
```bash
fetch_ucsc.py mm10 fa mm10.fa
fetch_ucsc.py mm10 ref mm10_ref.txt # optional
fetch_ucsc.py mm10 kg mm10_kg.txt # optional
```
2. Convert gene annotation file to GTF file using genePredToGtf. For this step, you will need genePredToGtf from UCSC Utilities (http://hgdownload.soe.ucsc.edu/admin/exe/).
```
cut -f2-11 mm10_kg.txt|./data/reference/genePredToGtf file stdin mm10_kg.gtf
```
3. If you downloaded multiple gene annotation files in step 1,  you can concatenate all of them into one file.
```
cat hg19_ref.txt hg19_kg.txt hg19_ens.txt > hg19_ref_all.txt
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

Extract the mature sequences from the mirbase file downloaded before and get the hairpin sequences. This example applies to mouse (mmu). For other organisms, please use the respective three letter code  (for example hsa for human).  
Also extract mature miRNA sequences from related species (here rat (rno) and human (hsa)).
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
│   └─── dataset.tsv
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
Replace every parameter in the file ```input/parameters.txt``` with information suitable for your dataset. If your totalRNA sequencing data is paired-end, set the paramerer ```paired_end=true```, if it is single-end set ```paired_end=false```.
```bash
dataset=input_folder/dataset.tsv
adapter="TGGAATTCTCGGGTGCCAAGG"
species=mmu
ref_dir=references/
ref_prefix=references/mm10
out_dir=output_folder
scripts_dir=/bin/scripts/
paired_end=false
```

##### dataset
The ```input/dataset.tsv``` file contains the paths to the totalRNA and smallRNA fastq files. The file should be tab-separated without header. The file structure depends on whether your totalRNA sequencing data is single-end or paired-end. Make sure to specify this in the ```parameters.tsv``` file.

####### single-end
If your totalRNA sequencing data is single-end, the ```dataset.tsv``` file should look like this:
Column 1: sample name
Column 2: path to totalRNA fastq file corresponding to the sample mentioned in column 1
Column 3: path to the miRNA fastq file corresponding to the sample mentioned in column 1
```
cerebellum_rep1 path/to/<totalRNA_sample1>.fastq	path/to/<smallRNA_sample1>.fastq
cerebellum_rep2 path/to/<totalRNA_sample2>.fastq	path/to/<smallRNA_sample2>.fastq
cerebellum_rep3	path/to/<totalRNA_sample3>.fastq	path/to/<smallRNA_sample3>.fastq
hippocampus_rep1  path/to/<totalRNA_sample4>.fastq path/to/<smallRNA_sample4>.fastq
...
```

####### paired-end
If your totalRNA sequencing data is single-end, the ```dataset.tsv``` file should look like this:
Column 1: sample name
Column 2: path to totalRNA read1 fastq file (corresponding to the sample mentioned in column 1)
Column 3: path to totalRNA read2 fastq file (corresponding to the sample mentioned in column 1)
Column 4: path to the miRNA fastq file (corresponding to the sample mentioned in column 1)
```
cerebellum_rep1 path/to/<totalRNA_sample1_R1>.fastq	path/to/<totalRNA_sample1_R2>.fastq path/to/<smallRNA_sample1>.fastq
cerebellum_rep2 path/to/<totalRNA_sample2_R1>.fastq	path/to/<totalRNA_sample2_R2>.fastq  path/to/<smallRNA_sample2>.fastq
cerebellum_rep3	path/to/<totalRNA_sample3_R1>.fastq	path/to/<totalRNA_sample3_R2>.fastq path/to/<smallRNA_sample3>.fastq
hippocampus_rep1  path/to/<totalRNA_sample4_R1>.fastq path/to/<totalRNA_sample4_R2>.fastq path/to/<smallRNA_sample4>.fastq
...
```

##### adapter
The adapter depends on the sequencing library used for your data. In order to find out which library adapter was used, check: https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
For example TruSeq Small RNA uses the adapter "TGGAATTCTCGGGTGCCAAGG".

##### species
Three letter code for the species used in this experiment. For example mmu (mouse), hsa(human).

##### ref_dir and ref_prefix
ref_dir is the directory containing all references. ref_prefix is the reference directory followed by the prefix used for building the bowtie index.

##### scripts_dir
Path to the ```scripts``` folder. ```scripts_dir=/scripts/``` for the use of Docker Container. If you prefer to run the pipeline manually, please specify the location of the scripts.

## Usage 
```bash
bash pipeline.sh /path/to/parameters.txt
```

## Build Docker
Download the ```Dockerfile``` and the folders ```scripts``` and ```soft```. In the same directory run:

```docker build -t pipeline_docker .```
