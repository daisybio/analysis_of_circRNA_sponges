library(plyr)
library(dplyr)
library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
output_dir = args[2]

dataset <- read.table(dataset_path, sep = "\t", header=F, stringsAsFactors = F)
samples <- dataset$V1

raw_bindSites <- read.table(paste0(output_dir, "/results/binding_sites/output/bindsites_25%_filtered.tsv"), header = T, sep = "\t", stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]
allBindSites <- count(bindSites, Target, name="freq")
bindsitDT  <- data.table(bindSites)

zero_counts_sample_cutoff <- ceiling(0.2*nrow(dataset))

miRNA_expression_raw <- read.table(paste0(output_dir, "/results/miRNA/miRNA_counts_all_samples_raw.tsv"), header = T, stringsAsFactors = F)
miRNA_expression <- miRNA_expression_raw[rowSums(miRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]

circRNA_expression_raw <- read.table(paste0(output_dir, "/results/circRNA/circRNA_counts_all_samples_filtered.tsv"),header = T, stringsAsFactors = F)
circRNA_expression <- circRNA_expression_raw[rowSums(circRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]

correlations <- data.frame(circRNA = character(),
                           miRNA = character(), 
                           circRNA_miRNA_ratio = numeric(),
                           miRNA_binding_sites = numeric(),
                           pearson_R = numeric(), 
                           corr_pval = numeric(),
                           RSS_norm = numeric(),
                           intercept = numeric(),
                           intercept_pval = numeric(),
                           slope = numeric(),
                           slope_pval = numeric(),
                           adj_r_squared = numeric(),
                           stringsAsFactors=FALSE) 
for (i in 1:nrow(circRNA_expression)){
  # get coordinations of current circRNA
  chr <- as.character(circRNA_expression[i,1])
  start <- as.numeric(as.character(circRNA_expression[i,2]))
  end <- as.numeric(as.character(circRNA_expression[i,3]))
  strand <- as.character(circRNA_expression[i,4])
  circRNA <- paste(chr,":", start, "-", end, "_", strand, sep="")
  
  # get sample counts for current circRNA
  circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:(5+length(samples)-1))]))
  colnames(circRNA_counts) <- "circRNA_counts"
  circRNA_counts$sample <- row.names(circRNA_counts)
  circRNA_counts$circRNA_counts <- as.numeric(as.character(circRNA_counts$circRNA_counts))
  
  for (j in 1:nrow(miRNA_expression)){
    miRNA <-as.character(miRNA_expression[j,1])
    
    # get sample counts for current miRNA
    miRNA_counts <- t(miRNA_expression[miRNA_expression$miRNA == miRNA,])
    miRNA_counts <- miRNA_counts[-1, ] 
    miRNA_counts <- as.data.frame(miRNA_counts)
    colnames(miRNA_counts) <- "miRNA_counts"
    miRNA_counts$sample <- row.names(miRNA_counts)
    miRNA_counts$miRNA_counts <- as.numeric(as.character(miRNA_counts$miRNA_counts))
    
    # compute circRNA expression vs. miRNA expression
    joined_counts <- merge(miRNA_counts, circRNA_counts, by="sample")
    
    # analyse circRNA/miRNA ratio
    mean_circRNA_counts <- mean(joined_counts$circRNA_counts)
    mean_miRNA_counts <- mean(joined_counts$miRNA_counts)
    circRNA_miRNA_ratio <- mean_circRNA_counts/mean_miRNA_counts
    
    # compute number of miRNA binding sites on circRNA
    mirna <-miRNA
    binding_sites <- nrow(bindsitDT[miRNA == mirna & Target == circRNA])
    
    # compute circRNA-miRNA correlation for all samples
    cor_res <- cor.test(joined_counts$miRNA_counts, joined_counts$circRNA_counts,  method = "pearson", use = "complete.obs")
    corr_R <- as.numeric(as.character(cor_res$estimate))
    corr_pval <- as.numeric(as.character(cor_res$p.value))
    
    # compute linear regression
    regression_model <- lm(miRNA_counts~circRNA_counts, data = joined_counts)
    intercept <- summary(regression_model)$coefficients[1,1]
    intercept_pval <- summary(regression_model)$coefficients[1,4]
    slope <- summary(regression_model)$coefficients[2,1]
    slope_pval <- summary(regression_model)$coefficients[2,4]
    adj_r_squared <- summary(regression_model)$adj.r.squared
    
    # compute residuals sum of squares
    # normalize counts for residuals sum of squares
    normalized_counts <- joined_counts[,c("circRNA_counts", "miRNA_counts")]
    min_circRNA_counts <- min(normalized_counts$circRNA_counts)
    max_circRNA_counts <- max(normalized_counts$circRNA_counts)
    normalized_counts[,"circRNA_counts"] <- (normalized_counts[,"circRNA_counts"] - min_circRNA_counts)/(max_circRNA_counts - min_circRNA_counts)
    min_miRNA_counts <- min(normalized_counts$miRNA_counts)
    max_miRNA_counts <- max(normalized_counts$miRNA_counts)
    normalized_counts[,"miRNA_counts"] <- (normalized_counts[,"miRNA_counts"] - min_miRNA_counts)/(max_miRNA_counts - min_miRNA_counts)
    norm_reg_model <- lm(miRNA_counts~circRNA_counts, data = normalized_counts)
    RSS_norm <- sum(norm_reg_model$residuals^2)
    
    # write correlation info in correlations data frame
    correlations[nrow(correlations) + 1,] <- c(circRNA, miRNA, circRNA_miRNA_ratio, binding_sites, corr_R, corr_pval, RSS_norm, intercept, intercept_pval, slope, slope_pval, adj_r_squared)
  }
}
write.table(correlations, file=paste0(output_dir, "/results/sponging/filtered_circRNA_miRNA_correlation.tsv"), sep = "\t", quote = F, row.names = F)
