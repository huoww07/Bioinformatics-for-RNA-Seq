## R scripts for differential expression
## These scripts are used to calculate differential expression using featurecounts data

```markdown
# set to work directory: make sure to set directory to where you store featurecounts result
setwd("./featurecounts/")

# load shared Tufts bio library path; if you don't have access to tufts HPC, skip this step
.libPaths('/cluster/tufts/bio/tools/R_libs/3.5')

# load required libraries
library(DESeq2)
library(vsn)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(DEGreport)
library(pheatmap)


# load data.
meta <- read.table("raw_data/sample_info.txt", header=TRUE)
feature_count <- read.table("./featurecounts/featurecounts_results.txt",
                            header=TRUE, row.names = 1)
# check to make sure that all rows labels in meta are columns in data
all(colnames(data) == rownames(meta))

# Differntial expression using DESeq2
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
dds <- DESeq(dds)
# Creating contrasts and running a Wald test
contrast <- c("condition", "SNF2", "WT")
res_unshrunken <- results(dds, contrast=contrast)
summary(res_unshrunken)
# Shrinkage of the log2 fold changes
res <- lfcShrink(dds, contrast=contrast, res=res_unshrunken)
summary(res)
# Filtering to find significant genes using FDR cutoff of 0.05
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results <- res[which(res$padj < padj.cutoff),]
# save results using customized file_name
file_name = 'significant_padj_0.05.txt'
write.table(significant_results, file_name, quote=FALSE)

# Visualization
# simple plot for a single gene YOR290C (SNF2)
plotCounts(dds, gene="YOR290C", intgroup="condition")
# plot multiple genes in a heatmap: top 50 with most significant padj value
significant_results_sorted <- significant_results[order(significant_results$padj), ]
significant_genes_50 <- rownames(significant_results_sorted[1:50, ])
# extract the counts from the rlog transformed object
rld_counts <- assay(rld)
# select by row name using the list of genes:
rld_counts_sig <- rld_counts[significant_genes_50, ]
# Plot multiple genes in a heatmap:
pheatmap(rld_counts_sig,
         cluster_rows = T,
         show_rownames = T,
         annotation = meta,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 8,
         height = 20)
```
