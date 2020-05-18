Approximate time: 20 minutes

## Learning Objectives

- Use the `featureCounts` function from the `Subread` package to perform feature quantification

<img src="../img/workflow_gene_quant.png" width="400">

## Counting reads: featureCounts

The mapped coordinates of each read are compared with the features in the GTF file. 
Reads that overlap with a gene by >=1 bp are counted as belonging to that feature. 
Ambiguous reads will be discarded and the output will be a matrix of genes and samples.

<img src="../img/featurecount_ambiguous.png" width="400">

By default featurecounts will 
1) count reads in features labeled as 'exon' in the GTF and 
2) group all exons with a given 'gene_id'. 

An example of a transcript with multiple exons:

<img src="../img/featurecount_multi_exons.png" width="600">

## Counting reads: running the script

Create a new script called `featureccounts.sh` using the nano text editor `nano featurecounts.sh` and enter the following 
content:

```
## load subread module
module load subread/1.6.3

## create output directory
mkdir featurecounts

## reference directory
REF_DIR=/cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3

## Run featurecounts
featureCounts \
-a ${REF_DIR}/Annotation/Genes/sacCer3.gtf \
-o featurecounts/featurecounts_results.txt \
STAR/*bam
```

To run the script, type in:
```
sh scripts/featurecounts.sh
```

Result:
```

        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
	  v1.6.3

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 1 BAM file                                       ||
||                           S WT_1_Aligned.sortedByCoord.out.bam             ||
||                                                                            ||
||             Output file : featurecounts_results.txt                        ||
||                 Summary : featurecounts_results.txt.summary                ||
||              Annotation : sacCer3.gtf (GTF)                                ||
||      Dir for temp files : featurecounts                                    ||
||                                                                            ||
||                 Threads : 1                                                ||
||                   Level : meta-feature level                               ||
||              Paired-end : no                                               ||
||      Multimapping reads : not counted                                      ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\===================== http://subread.sourceforge.net/ ======================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file sacCer3.gtf ...                                       ||
||    Features : 7050                                                         ||
||    Meta-features : 6692                                                    ||
||    Chromosomes/contigs : 17                                                ||
||                                                                            ||
|| Process BAM file WT_1_Aligned.sortedByCoord.out.bam...                     ||
||    Single-end reads are included.                                          ||
||    Assign alignments to features...                                        ||
||    Total alignments : 6014703                                              ||
||    Successfully assigned alignments : 5395145 (89.7%)                      ||
||    Running time : 0.12 minutes                                             ||
||                                                                            ||
||                                                                            ||
|| Summary of counting results can be found in file "featurecounts/featureco  ||
|| unts_results.txt.summary"                                                  ||
||                                                                            ||
\\===================== http://subread.sourceforge.net/ ======================//
```

The output files will contain results and results summary.
```
featurecounts/
├── featurecounts_results.txt
└── featurecounts_results.txt.summary
```

## Counting reads: Result summary
To look at the result, type:
```
cat featurecounts/featurecounts_results.txt.summary | column -t
```

The code after the bash pipe `|` serves to provide nice tab-separated formatting.
To see the effect of `| column -t`, try removing it and compare the output.

Result:
```
Status                         STAR/WT_1_Aligned.sortedByCoord.out.bam
Assigned                       5395145
Unassigned_Unmapped            0
Unassigned_MappingQuality      0
Unassigned_Chimera             0
Unassigned_FragmentLength      0
Unassigned_Duplicate           0
Unassigned_MultiMapping        0
Unassigned_Secondary           0
Unassigned_Nonjunction         0
Unassigned_NoFeatures          232980
Unassigned_Overlapping_Length  0
Unassigned_Ambiguity           386578
```

We see that 5395145 reads have been assigned to the features present in our gtf file.

## Tracking read numbers
As the analysis progresses you should keep track of the following:

<img src="../img/featurecount_read_summary.png" width="400">

## Workshop Schedule
- [Introduction](../README.md)
- [Setup using Tufts HPC](01_Setup.md)
- [Process Raw Reads](02_Quality_Control.md)
- [Read Alignment](03_Read_Alignment.md)
- Currently at: Gene Quantification
- Next: [Differential Expression](05_Differential_Expression.md)
- [Pathway Enrichment](06_Pathway_Enrichment.md)
