## Dependencies for RNA sequencing workshop

## Description:
To follow along with the practice in this course, you will use command line tools and R installed on your computer. Additionally, you will have to install the softwares and R libraries we will be using.

## Full list of Dependencies:
1. Download [Raw data for this course](https://tufts.box.com/v/intro-to-RNA-seq-material)

2. Command line tools for preprocessing reads and read alignment:
[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[STAR](https://github.com/alexdobin/STAR)
[samtools](http://www.htslib.org/download/)
[subread](http://subread.sourceforge.net)

 3. Download and install [Interactive Genome Viewer (IGV)](https://software.broadinstitute.org/software/igv/download)

 4. R and its libraries for differential expression analysis:
[R studio](https://www.google.com/search?client=safari&rls=en&q=r+studio&ie=UTF-8&oe=UTF-8)

    Bioconductr libraries (check out each library for installation guide):
    - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    - [vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)
    - [DEGreport](https://bioconductor.org/packages/release/bioc/html/DEGreport.html)
    - [org.Sc.sgd.db](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)
    - [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

    Regular libraries (install by entering `install.packages("name")`):

    - [tidyverse](https://ggplot2.tidyverse.org)
    - [ggplot2](https://ggplot2.tidyverse.org)
    - [dplyr](https://dplyr.tidyverse.org)
    - [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html#installation)
    - pheatmap

## Scripts you will be using:
- [Command line scripts](08_bash_scripts.md)
- [R scripts](09_R_scripts.md)


## Workshop Schedule
- [Introduction](../README.md)
- [Setup using Tufts HPC](01_Setup.md)
- [Process Raw Reads](02_Quality_Control.md)
- [Read Alignment](03_Read_Alignment.md)
- [Gene Quantification](04_Gene_Quantification.md)
- [Differential Expression](05_Differential_Expression.md)
- [Pathway Enrichment](06_Pathway_Enrichment.md)
