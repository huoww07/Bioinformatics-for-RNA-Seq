Approximate time: 20 minutes

## Goals
- Connect to the HPC cluster via On Demand Interface
- Download data

# Log into the HPC cluster's On Demand interface
1. Open a Chrome browser visit [ondemand.cluster.tufts.edu](ondemand.cluster.tufts.edu)
2. Log in with your Tufts Credentials
3. On the top menu bar choose Clusters->HPC Shell Access

<img src="../img/od_terminal.png" width="400">

4. Type your password at the prompt (the password will be hidden for security purposes):

`whuo01@login.cluster.tufts.edu's password:`

5. You'll see a welcome message and a bash prompt, for example for user `whuo01`:

`[whuo01@login001 ~]$`

This indicates you are logged in to the login node.

6. Type `clear` to clear the screen

# Set up for the analysis

## Find 500M storage space

1. Check how much available storage you have in your home directory by typing `showquota`.

Result:
```
Home Directory Quota
Disk quotas for user whuo01 (uid 31394):
     Filesystem  blocks   quota   limit   grace   files   quota   limit   grace
hpcstore03:/hpc_home/home
                  1222M   5120M   5120M            2161   4295m   4295m        


Listing quotas for all groups you are a member of
Group: facstaff	Usage: 16819478240KB	Quota: 214748364800KB	Percent Used: 7.00%
```

Under `blocks` you will see the amount of storage you are using, and under quota you see your quota.
Here, the user has used 1222M/5120M and has enough space for our 500M analysis.

2. If you do not have 500M available, you may have space in a project directory for your lab.
These are located in `/cluster/tufts` with names like `/cluster/tufts/labname/username/`.
If you don't know whether you have project space, please email [tts-research@tufts.edu](mailto:tts-research@tufts.edu).

## Download the data
1. Get an interaction session on a compute node by typing:

`srun --pty -t 3:00:00  --mem 16G  -N 1 -n 4 bash`

Once you hit enter, you will see something like below showing it is requesting for the session:
```
(base) [whuo01@login001 ~]$ srun --pty -t 3:00:00  --mem 16G  -N 1 -n 4 bash
srun: job 55918493 queued and waiting for resources
```
If wait times are very long, you can try a different partitions by adding, e.g. -p preempt or -p interactive before bash.

```
(base) [whuo01@login001 ~]$ srun --pty -t 3:00:00  --mem 16G  -N 1 -n 4 -p preempt bash
(base) [whuo01@pcomp45 ~]$
```

The success is indicated by the change of environment name after your username. Here it was changed from `login001` to `pcomp45`. This is an indication that you may proceed to the next step.

Note: If you go through this workshop in multiple steps, you will have to rerun this step each time you log in.


2. Change to your home directory
`cd ~`
Or, if you are using a project directory:
`cd /cluster/tufts/labname/username/`

3. Copy the course directory:
`cp /cluster/tufts/bio/tools/training/bioinformatics-for-rna-seq/intro-to-RNA-seq.tar.gz ./`

4. Unzip the course directory:
`tar -xvzf intro-to-RNA-seq.tar.gz`

5. Take a look at the contents by typing:
`tree intro-to-RNA-seq`

You'll see a list of all files
```
intro-to-RNA-seq/
├── ERP004763_info.txt                 <-- sample description
├── raw_data                           <-- Folder with fastq files
│   ├── sample_info.txt
│   ├── SNF2
│   │   ├── ERR458500.fastq.gz         <-- gzip compressed fastq files
│   │   ├── ERR458501.fastq.gz
│   │   ├── ERR458502.fastq.gz
│   │   ├── ERR458503.fastq.gz
│   │   ├── ERR458504.fastq.gz
│   │   ├── ERR458505.fastq.gz
│   │   └── ERR458506.fastq.gz
│   └── WT_1
│       ├── ERR458493.fastq.gz
│       ├── ERR458494.fastq.gz
│       ├── ERR458495.fastq.gz
│       ├── ERR458496.fastq.gz
│       ├── ERR458497.fastq.gz
│       ├── ERR458498.fastq.gz
│       └── ERR458499.fastq.gz
└── scripts                           <-- Folder with all commands
    ├── fastqc.sh
    ├── featurecounts.sh
    ├── intro.R
    ├── sbatch_star_align_individual.sh
    ├── sbatch_star_align.sh
    └── sbatch_star_align_SNF2.sh

4 directories, 22 files
```

## Data for the class

Project: European Nucleotide Archive Project number PRJEB5348

Samples: WT v.s. SNF2 (a snf2 knock-out mutant cell line)

Organism: Saccharomyces cerevisiae

Sequencing: Illumina HiSeq, Single End, 50bp read length


## Workshop Schedule
- [Introduction](../README.md)
- Currently at: Setup using Tufts HPC
- Next: [Quality Control](02_Quality_Control.md)
- [Read Alignment](03_Read_Alignment.md)
- [Gene Quantification](04_Gene_Quantification.md)
- [Differential Expression](05_Differential_Expression.md)
- [Pathway Enrichment](06_Pathway_Enrichment.md)
