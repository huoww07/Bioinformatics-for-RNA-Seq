## A Tufts University Research Technology Workshop

## Description:
This course will cover the basics for RNA sequencing and analysis using command line tools and R on Tufts High Performance Compute Cluster (HPC). Topics covered include quality control, read alignment, feature counts, differential expression and pathway analysis. There will be a 1-hour introductory session after which students will follow online self-guided online material. The material is designed to run on [Tufts High Performance Compute (HPC) Cluster](https://access.tufts.edu/research-cluster-account). If you don't have access to Tufts HPC, you will need to install [all required modules on your own computer](lessons/07_dependencies.md).

## Goals
- Writing and running bash scripts
- Intro to RNA sequencing logistics
- Intro to bioinformatics toos: BWA, samtools
- Differential expression analysis using R and DESeq2

<img src="img/workflow.png" width="400">

## Prerequisites

### Computational skills needed
- Introduction to [Linux](https://tufts.box.com/s/x9aflewr2qw59pcbgcghbo9muykbi4ju)
- Introduction to [HPC](https://tufts.box.com/s/yubnzxnpih14hd80mbfxqrkdri8s2nws)
- Introduction to [R](https://learn.datacamp.com/courses/free-introduction-to-r)

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


## Schedule
- Introduction
- [Setup using Tufts HPC](lessons/01_Setup.md)
- [Process Raw Reads](lessons/02_Process_Raw_Reads.md)
- [Read Alignment](lessons/03_Read_Alignment.md)
- [Gene Quantification](lessons/04_Gene_Quantification.md)
- [Differential Expression](lessons/05_Differential_Expression.md)
- [Pathway Enrichment](lessons/06_Pathway_Enrichment.md)


## Acknowledgement
Much of this workshop was adapted from [Bioinformatics @ Tufts ](https://sites.tufts.edu/biotools/tutorials/) and the [HBC Training DGE workshop](https://github.com/hbctraining/DGE_workshop) with the help of Dr. Rebecca Batorsky and Dr. Albert Tai at Tufts University.
