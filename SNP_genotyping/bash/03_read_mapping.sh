#!/bin/bash

### Read mapping of WGS data from 99 Littorina saxatilis samples. Each sample was
### seqeunced on seperate lanes of BGI Genomic's DNBSEQ platform. This script 
### uses a nested for loop to navigate to direcotories containing fastq files that
### are split between forward and reverse stands and different lanes. The next loop
### finds all files that were run on the same lane and runs BWA-MEM to call reads.
### This job was run on the server 'Rackham', part of the National Academic 
### Infrastructure for Supercomputing in Sweden.

### James Reeve - Göteborgs Universitet
### 09/08/2021

### Job parameters
#SBATCH -A snic2021-5-142
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 30:00:00
#SBATCH --array=1-99
#SBATCH -J KT_read_map_v2
#SBATCH --output=KT_logfiles/mapping_run4/KT_read_mapping.out

### Required software
module load bioinfo-tools
module load samtools/1.12
module load bwa/0.7.17

### The script:  
# 1. Filepaths:
KT_DIR=/proj/snic2020-6-155

# 2. Index reference genome:
# Note: this line only needs to be run once!
#bwa index $KT_DIR/ref_genome/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta

# 3. Set array flag to sample names:
samplesheet=$KT_DIR/snailIDs.txt
SnailID=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`

# 4. Run BWA-MEM on each file
### Read group information
lib=$(ls $KT_DIR/fastq/$SnailID | head -n 1 | cut -d "_" -f1)
pltunit=$(ls $KT_DIR/fastq/$SnailID | head -n 1 | cut -d "_" -f3)
###     
### Read mapping
bwa mem -M -t 2 -R "@RG\tSM:${SnailID}\tLB:${lib}\tPL:DNDSEQ\tPU:$pltunit" \
        $KT_DIR/ref_genome/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta \
        $KT_DIR/fastq/$SnailID/${SnailID}_1.fq.gz \
        $KT_DIR/fastq/$SnailID/${SnailID}_2.fq.gz | \
        samtools view -b > $TMPDIR/${SnailID}_PE.bam


# 4.5 Error check
for i in $TMPDIR/${SnailID}_PE.bam
do
        samtools view -H $i | head
done

# 5. Sort BAM and generate summary file
samtools sort -o $KT_DIR/bam/${SnailID}_sorted.bam $TMPDIR/${SnailID}_PE.bam
samtools flagstat $KT_DIR/bam/${SnailID}_sorted.bam > $KT_DIR/bam/${SnailID}.flagstat
