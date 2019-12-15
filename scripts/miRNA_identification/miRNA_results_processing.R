#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
output_dir = args[2]

library(ggplot2)

dataset <- read.table(dataset_path, sep = "\t", header=F, stringsAsFactors = F)
samples <- dataset$V1

dataset_counts_raw <- NULL
dataset_counts_norm <- NULL
for(i in 1:length(samples)){
  # get results
  sample <- samples[i]
  sample_folder <- paste0(output_dir, "/samples/",  sample, "/miRNA_identification/")
  miRDeep2_output <- read.table(dir(sample_folder, full.names=T, pattern="miRNAs_expressed"), sep = "\t", header=F, stringsAsFactors = F)
  names(miRDeep2_output) <- c("miRNA", "counts", "precursor", "total", "seq", "norm")
  raw_counts <- miRDeep2_output[,c(1,2)]
  norm_counts <- miRDeep2_output[,c(1,6)]
  
  # sum up counts which come from the same miRNA because of different precursors
  raw_counts <- aggregate(raw_counts$counts, by=list(miRNA=raw_counts$miRNA), FUN=sum)
  names(raw_counts) <- c("miRNA", sample)
  norm_counts <- aggregate(norm_counts$norm, by=list(miRNA=norm_counts$miRNA), FUN=sum)
  names(norm_counts) <- c("miRNA", sample)
  if(is.null(dataset_counts_raw)){
    dataset_counts_raw <- raw_counts
    dataset_counts_norm <- norm_counts
  } else {
    dataset_counts_raw <- merge(dataset_counts_raw, raw_counts, all = T, by = "miRNA")
    dataset_counts_norm <- merge(dataset_counts_norm, norm_counts, all = T, by = "miRNA")
    
  }
  
   nonzero_counts <- raw_counts[raw_counts[,sample] > 0,]
   colnames(nonzero_counts)[2] <- "reads"
   qplot(nonzero_counts$reads,
         geom="histogram",
         fill=I("red"), 
         col=I("red"),
         alpha=I(.2),
         binwidth=100,
         main=paste0("Read distribution for miRNA mapping (", sample, ")"),
         xlab="Read unnormalized counts (> 0)",
         ylab="Number of miRNAs")
   
  ggsave(paste0(output_dir, "results/miRNA/read_distributions/miRNA_read_distribution_", sample, ".png"), width = 8, height = 4)

}
dataset_counts_raw[is.na(dataset_counts_raw)]  <- 0
write.table(dataset_counts_raw, paste0(output_dir, "results/miRNA/miRNA_counts_all_samples_raw.tsv"), quote = F, sep = "\t", row.names = F)

dataset_counts_norm[is.na(dataset_counts_norm)]  <- 0
write.table(dataset_counts_norm, paste0(output_dir, "results/miRNA/miRNA_counts_all_samples_norm.tsv"), quote = F, sep = "\t", row.names = F)

 
  



