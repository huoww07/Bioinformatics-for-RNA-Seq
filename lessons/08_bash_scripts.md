## Bash scripts for RNA sequencing
## These scripts are used to process raw reads, align processed reads and quantify gene expression using feature counts

```markdown
# For Tufts HPC users who attended our 1-hr zoom workshop, please run following command to get a compute node
srun -t 3:00:00 --mem 16G -N 1 -n 4 -p preempt --reservation bioworkshop --pty bash

# For other HPC users, please use below command
srun -t 3:00:00 --mem 16G -N 1 -n 4 -p preempt --pty bash

# Navigate to project folder
cd /cluster/tufts/bio/tools/training/intro-to-rnaseq/users/

# make a new directory using your username for your practice, and enter that directory
mkdir YOUR_USERNAME
cd YOUR_USERNAME
# copy course material to your directory, unzip it and enter the course material directory
cp /cluster/tufts/bio/tools/training/intro-to-rnaseq/intro-to-RNA-seq-May-2020.tar.gz ./
tar -xvzf intro-to-RNA-seq-May-2020.tar.gz
cd intro-to-RNA-seq

# load fastqc module
module load fastqc/0.11.8
mkdir fastqc

# Run fastqc on raw sequencing reads. Each fastq file will be analyzed individually. * is a wild card. 
fastqc raw_data/WT/*.fastq.gz -o fastqc --extract
fastqc raw_data/SNF2/*.fastq.gz -o fastqc --extract

# run multiqc to compile individual fastqc files, this helps visualization of fastqc reports
module load multiqc/1.7.0
mkdir multiqc
multiqc fastqc/ -o multiqc

# Read alignment Step 1: prepare reference genomes. Your input will be genome.fa.
module load STAR/2.6.1d
mkdir genome
STAR --runMode genomeGenerate --genomeDir ./genome --genomeFastaFiles /cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa --runThreadN 4

# Read alignment Step 2: align. You will need your own annotation file in gtf format. You will run this step for individual samples. Here we are using ERR458493.fastq.gz as an example.
mkdir STAR
STAR --genomeDir ./genome \
--readFilesIn raw_data/WT/ERR458493.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix STAR/ \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 4 \
--alignIntronMin 1 \
--alignIntronMax 2500 \
--sjdbGTFfile annotation.gtf \
--sjdbOverhang 49

# Read alignment Step 3: generate bam index. The output.bam is output file from step 2. You will run this for individual samples following Step 2.
module load samtools/1.2
samtools index STAR/output.bam


# Before moving on to next step, make sure your STAR folder contains 14 bam files, one for each replicate.
# You can do so by running
sh ./scripts/star_align_individual.sh


# Gene quantification using featureCounts - This step compiles all alignment results together. This is done after alignment is finished for all samples.
# Gene quantification step 1: load subread module
module load subread/1.6.3
# Gene quantification step 2: create output directory
mkdir featurecounts
# Gene quantification step 3: Run featurecounts
featureCounts \
-a /cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/sacCer3.gtf \
-o featurecounts/featurecounts_results.txt \
STAR/*.bam


# Check feature count results
cat featurecounts/featurecounts_results.txt.summary
head featurecounts/featurecounts_results.txt
# Clean the column names in featurecounts_results.txt
cat featurecounts/featurecounts_results.txt |sed "2s/STAR\///g" | sed "2s/\_Aligned.sortedByCoord.out.bam//g" > featurecounts/featurecounts_results.mod.txt

# now you are ready to move on to R scripts.



# Optional step 1. Visualize number of mapped reads v.s. unmapped reads in all samples using barplot. This code will generate a pdf file named Mapping_stat.pdf.
module load R/3.5.0
Rscript ./scripts/mapping_percentage.R

# Optional step 2. Visualize number of assigned reads in all samples using barplot. This code will generate a pdf file named Featurecount_stat.pdf.
module load R/3.5.0
Rscript ./scripts/featurecount_stat.R

```

- [Next: R scripts](09_R_scripts.md)

## Go back to workshop schedule
- [Introduction](../README.md)
- [Setup using Tufts HPC](01_Setup.md)
- [Process Raw Reads](02_Quality_Control.md)
- [Read Alignment](03_Read_Alignment.md)
- [Gene Quantification](04_Gene_Quantification.md)
- [Differential Expression](05_Differential_Expression.md)
- [Pathway Enrichment](06_Pathway_Enrichment.md)
