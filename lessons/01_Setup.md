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

2. Change to your home directory
`cd ~`
Or, if you are using a project directory:
`cd /cluster/tufts/labname/username/`

3. Copy the course directory:
`cp /cluster/tufts/bio/tools/intro-to-ngs.tar.gz .`

4. Unzip the course directory:
`tar -xvzf intro-to-ngs.tar.gz`

5. Take a look at the contents by typing:
`tree intro-to-ngs`

You'll see a list of all files
```
intro-to-ngs
├── all_commands.sh          <-- Bash script with all commands
├──intro_to_ngs_Dec2019.pdf  <-- Course slides
├── raw_data                 <-- Folder with paired end fastq files
│   ├── na12878_1.fq         
│   └── na12878_2.fq
├── README.md                <-- Instructions
└── ref_data                 <-- Folder with reference sequence
    └── chr10.fa
2 directories, 5 files
```

## Data for the class

Project: European Nucleotide Archive Project number PRJEB5348

Samples: WT v.s. SNF2 (a snf2 knock-out mutant cell line)

Organism: Saccharomyces cerevisiae

Sequencing: Illumina HiSeq, Single End, 50bp read length

[Previous: Introduction](00_Introduction.md)

[Next: Quality Control](02_Quality_Control.md)
