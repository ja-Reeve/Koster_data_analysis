#!/bin/bash

### Mergeing fastq data from 99 Littorina saxatilis samples. Each sample is split
### into multiple files, depending on the sequencing lane. This has been a hassel
### when read mapping, so these seperate reads are concatonated together to create
### a forward read and reverse read file for each individual, that is labeled with
### with the sample ID. This job was run on the server 'Rackham', part of the 
### National Academic Infrastructure for Supercomputing in Sweden.
### James Reeve - Göteborgs Universitet
### 16/09/2021

### Job parameters
#SBATCH -A snic2021-5-142
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --array=1-99
#SBATCH -J KT_merge_fq
#SBATCH --output=KT_logfiles/KT_merge_fq.out

### The script
# 1. Filepaths:
KT_DIR=/proj/snic2020-6-155

# 2. Create list of samples:
basename -a $KT_DIR/fastq/KT* > $KT_DIR/snailIDs.txt
samplesheet=$KT_DIR/snailIDs.txt
SnailID=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`

# 3. Merge fastq:
echo "\nProgress: Merging fastq files for $SnailID"
## Error check: read numbers prior to merger
zcat $KT_DIR/fatq/$SnailID/*_1.fq.gz | echo "Total forward read count for $SnailID = $((`wc -l`/4))"
zcat $KT_DIR/fatq/$SnailID/*_2.fq.gz | echo "Total reverse read count for $SnailID = $((`wc -l`/4))"

## Merge forward reads
cat $KT_DIR/fastq/$SnailID/*_1.fq.gz > $KT_DIR/fastq/${SnailID}_1.fq.gz
zcat $KT_DIR/fatq/${SnailID}_1.fq.gz | echo "Forward read count after merger = $((`wc -l`/4))"

## Merge reverse reads
cat $KT_DIR/fastq/$SnailID/*_2.fq.gz > $KT_DIR/fastq/${SnailID}_2.fq.gz
zcat $KT_DIR/fatq/${SnailID}_2.fq.gz | echo "Reverse read count after merger = $((`wc -l`/4))"
