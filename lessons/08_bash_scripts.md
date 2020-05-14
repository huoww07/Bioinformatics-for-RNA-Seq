## Bash scripts for RNA sequencing
## These scripts are used to process raw reads, align processed reads and quantify gene expression using feature counts

```markdown
# load fastqc module
module load fastqc/0.11.8
mkdir fastqc

# Run fastqc on raw sequencing reads on input1.fastq.gz and input2.fastq
fastqc input1.fastq.gz -o fastqc --extract
fastqc input2.fastq -o fastqc --extract

# run multiqc to compile individual fastqc files
module load multiqc/1.7.0
mkdir multiqc
multiqc fastqc/ -o multiqc

# Read alignment Step 1: prepare reference genomes. Your input will be genome.fa.
module load STAR/2.6.1d
mkdir genome
STAR --runMode genomeGenerate --genomeDir ./genome --genomeFastaFiles genome.fa --runThreadN 12

# Read alignment Step 2: align. You will need your own annotation file in gtf format. You will run this step for individual samples.
mkdir STAR
STAR --genomeDir ./genome \
--readFilesIn input1.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix STAR/ \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 4 \
--alignIntronMin 1 \
--alignIntronMax 2500 \
--sjdbGTFfile annotation.gtf \
--sjdbOverhang 49

# Read alignment Step 3: generate bam index. The output.bam is output file from step 2.
module load samtools/1.2
samtools index STAR/output.bam


# Gene quantification using featureCounts
# Gene quantification step 1: load subread module
module load subread/1.6.3
# Gene quantification step 2: create output directory
mkdir featurecounts
# Gene quantification step 3: Run featurecounts
featureCounts \
-a /directory/to/where/you/store/annotation.gtf \
-o featurecounts/featurecounts_results.txt \
STAR/*bam

# now you are ready to move on to R scripts.
```
