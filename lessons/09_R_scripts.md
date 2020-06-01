## R scripts for differential expression
## These scripts are used to calculate differential expression using featurecounts data

```markdown
# set to work directory: make sure to set directory to project folder. Change whuo01 to your own username below
setwd("/cluster/tufts/bio/tools/training/intro-to-rnaseq/users/whuo01/intro-to-RNA-seq/")

# load shared Tufts bio library path; if you don't have access to tufts HPC, skip this step
.libPaths('/cluster/tufts/bio/tools/R_libs/3.5')

# load required libraries
library(DESeq2)
library(vsn)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(DEGreport)
library(pheatmap)


# load data.
meta <- read.table("raw_data/sample_info.txt", header=TRUE)
feature_count <- read.table("./featurecounts/featurecounts_results.mod.txt",
                            header=TRUE, row.names = 1)
# remove first 6 columns by select the column 6 to 19
data <- feature_count[,6:19]

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

# Quality control: Principle Components Analysis
# regularized log transformation (rlog)
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition") + geom_text(aes(label=name))

# Filtering to find significant genes using FDR cutoff of 0.05
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results <- res[which(res$padj < padj.cutoff),]
# save results using customized file_name
file_name = 'significant_padj_0.05.txt'
write.table(significant_results, file_name, quote=FALSE)



# Visualization

# simple plot for a single gene YOR290C (SNF2)
plotCounts(dds, gene="YOR290C", intgroup="condition")

# heatmap  
# plot top 50 genes in a heatmap: top 50 with most significant padj value
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


# volcano plot
# load necessary library ggplot2
library(ggplot2)
# add another column in the results table to label the significant genes using threshold of padj<0.05 and absolute value of log2foldchange >=1
res_table <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_table <- res_table %>%
  mutate(threshold_OE =  padj < 0.05 & abs(log2FoldChange) >= 1)
# you can view the modified table
view(res_table)
# make volcano plot, the significant genes will be labeled in red
ggplot(res_table) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  scale_color_manual(values=c("black", "red")) +  # black v.s. red dots
  ggtitle("SNF2 against WT") +                       # this line defines the title of the plot
  xlab("log2 fold change") +                      # this line defines the name of the x-axis
  ylab("-log10 adjusted p-value") +               # name of y-axis
  scale_x_continuous(limits = c(-7.5,7.5)) +      # the axis range is set to be from -7.5 to 7.5
  theme(legend.position = "none", #c(0.9, 0.9),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))




# functional analysis using clusterprofiler
# load library
library(org.Sc.sgd.db)
library(clusterProfiler)

## Run GO enrichment analysis for the top 500 genes
significant_results_sorted <- res[order(res$padj), ]
significant_genes_500 <- rownames(significant_results_sorted[1:500, ])
ego <- enrichGO(gene = significant_genes_500,
                         keyType = "ENSEMBL",
                         OrgDb = org.Sc.sgd.db)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
## Dotplot
dotplot(ego, showCategory=50)
## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50)
```

- [Review bash scripts](08_bash_scripts.md)

## Go back to workshop schedule
- [Introduction](../README.md)
- [Setup using Tufts HPC](01_Setup.md)
- [Process Raw Reads](02_Quality_Control.md)
- [Read Alignment](03_Read_Alignment.md)
- [Gene Quantification](04_Gene_Quantification.md)
- [Differential Expression](05_Differential_Expression.md)
- [Pathway Enrichment](06_Pathway_Enrichment.md)
