### A Tufts University Research Technology Workshop

## Description:
This tutorial will cover the basics of bioinformatics for RNA sequencing analysis using command line tools and R on the Tufts  High Performance Compute Cluster (HPC). The analysis pipeline includes quality control, read alignment, feature quantification, differential expression and pathway analysis. There was a [1-hour introductory session](slides/pipeline_walthrough.pdf), after which students will follow this self-guided online material. The material is designed to run on [Tufts High Performance Compute (HPC) Cluster](https://access.tufts.edu/research-cluster-account). If you don't have access to Tufts HPC, you will need to install [all required modules on your own computer](lessons/07_dependencies.md).

## Goals
- Intro to RNA sequencing logistics
- Write and run bash scripts
- Intro to command line bioinformatics tools: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [STAR](https://github.com/alexdobin/STAR), [SAMtools](http://samtools.sourceforge.net/)
- Perform differential expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in [Rstudio](https://rstudio.com/)  

<img src="img/workflow.png" width="500">

## Prerequisites

### Basic understanding of computational skills
- Introduction to [Linux](slides/Intro_To_Basic_Linux_SHARED.pdf)
- Introduction to [HPC](slides/Tufts_HPC_Cluster_New_User_Guide.pdf)
- Introduction to R [workshop materials](https://tufts.app.box.com/v/IntroR) and [recording](https://sites.tufts.edu/datalab/recorded-data-lab-workshops/)

### Materials Needed
If you are a Tufts member and have access to Tufts HPC:
- Chrome web browser
- Account on [Tufts HPC](https://access.tufts.edu/research-cluster-account)
- [VPN](https://access.tufts.edu/vpn) if accessing the HPC from off campus

If you don't have access to Tufts HPC:
- [Here is the list of full dependencies you will need to install](lessons/07_dependencies.md)

All scripts you will be using:
- [Command line scripts](lessons/08_bash_scripts.md)
- [R scripts](lessons/09_R_scripts.md)

1-hr introduction session slides:
- [Introduction](slides/RNAseq_intro_RB_28May20.pdf)
- [Pipeline walkthrough](slides/pipeline_walthrough.pdf)

## Schedule
- Currently at: Course Home
- [Introduction](slides/RNAseq_intro_RB_28May20.pdf)
- [Setup using Tufts HPC](lessons/01_Setup.md)
- [Quality Control](lessons/02_Quality_Control.md)
- [Read Alignment](lessons/03_Read_Alignment.md)
- [Gene Quantification](lessons/04_Gene_Quantification.md)
- [Differential Expression](lessons/05_Differential_Expression.md)
- [Pathway Enrichment](lessons/06_Pathway_Enrichment.md)


## Acknowledgement
Much of this workshop was adapted from [Bioinformatics @ Tufts ](https://sites.tufts.edu/biotools/tutorials/) and the [HBC Training DGE workshop](https://github.com/hbctraining/DGE_workshop) with the help of Dr. Rebecca Batorsky and Dr. Albert Tai at Tufts University.
