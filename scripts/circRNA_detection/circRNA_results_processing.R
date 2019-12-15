#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
output_dir = args[2]

dataset <- read.table(dataset_path, sep = "\t", header=F, stringsAsFactors = F)
samples <- dataset$V1

alldata <- NULL
finaldata <- NULL
finaldata_filtered <- NULL
N_of_circRNAs_raw <- c()
N_of_circRNAs_filtered <- c()
cutoff <- 5
for (i in 1:length(samples)){
  sample <- samples[i]
  path <- paste0(output_dir,"/samples/", sample , "/circRNA_detection/CIRCexplorer2/circularRNA_known.txt")
  CIRCexplorer2_output <- read.table(path, sep = "\t")
  
  N_of_circRNAs_raw[i] <- nrow(CIRCexplorer2_output)
  names(N_of_circRNAs_raw)[i] <- sample
  
  # consider just circRNAs with more than 5 reads
  CIRCexplorer2_output_filtered <- CIRCexplorer2_output[CIRCexplorer2_output[,13]>cutoff,]
  N_of_circRNAs_filtered[i] <- nrow(CIRCexplorer2_output_filtered)
  names(N_of_circRNAs_filtered)[i] <- samples[i]
  
  expression_raw <- CIRCexplorer2_output[,c(1,2,3,6,13,15,16)]
  colnames(expression_raw) <- c("chr", "start", "stop", "strand", sample, paste("gene_", sample, sep = ""), paste("isoform_", sample, sep = ""))
  expression_raw <- expression_raw[,c(1,2,3,4,5)]
  
  expression_filtered <- CIRCexplorer2_output_filtered[,c(1,2,3,6,13,15,16)]
  colnames(expression_filtered) <- c("chr", "start", "stop", "strand", sample, paste("gene_", sample, sep = ""), paste("isoform_", sample, sep = ""))
  expression_filtered <- expression_filtered[,c(1,2,3,4,5)]
  if(is.null(finaldata)){
    finaldata <- expression_raw
    finaldata_filtered <- expression_filtered
    alldata <- expression_raw[,sample]
  } else {
    finaldata <- merge(finaldata, expression_raw, by = c("chr", "start", "stop", "strand"), all = T)
    finaldata_filtered<- merge(finaldata_filtered, expression_filtered, by = c("chr", "start", "stop", "strand"), all = T)
    alldata <- append(alldata, expression_raw[,sample])
    }
  
}
finaldata[is.na(finaldata)] <- 0
finaldata_filtered[is.na(finaldata_filtered)] <- 0
write.table(finaldata, paste0(output_dir, "/results/circRNA/circRNA_counts_all_samples_raw.tsv"), quote = F, sep = "\t", row.names = F)


# filter data (counts > 5) and circRNA present in minimum 15% of samples
if(length(samples)>=7){
sample_nr_cutoff <- floor(0.15*length(samples))
} else {
  sample_nr_cutoff <- 1
}
rows_to_keep <- c()
for (i in 1:nrow(finaldata_filtered)){
  number_of_samples_containing_this_circRNA <- 0
  for (j in 5:ncol(finaldata_filtered)){
    if(finaldata_filtered[i,j]>0){
      number_of_samples_containing_this_circRNA <- number_of_samples_containing_this_circRNA + 1
    }
  }
  if(number_of_samples_containing_this_circRNA >= sample_nr_cutoff){
    rows_to_keep <- append(rows_to_keep, i)
  }
  
}
filtered_data <- finaldata_filtered[rows_to_keep,]
write.table(filtered_data, paste0(output_dir, "/results/circRNA/circRNA_counts_all_samples_filtered.tsv"), quote = F, sep = "\t", row.names = F)

