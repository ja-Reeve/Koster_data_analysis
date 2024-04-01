#!/bin/bash

### Filtering of sequence alignemnt files (.bam) of 99 Littorina saxatilis
### snails collected from the Kosterharvets national park in Sweden. Two
### filters are used, 1) removing duplicates and 2) removing unpaired reads.
### This job was run on the server 'Rackham', part of the National Academic 
### Infrastructure for Supercomputing in Sweden.

### James Reeve - Göteborgs Universitet
### 09/08/2021

### Job parameters
#SBATCH -A snic2021-22-384
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 30:00:00
#SBATCH --array=1-99
#SBATCH -J KT_mark_duplicates
#SBATCH --output=KT_logfiles/duplicates/KT_mark_duplicates_%a.out

### Required software
module load bioinfo-tools
module load samtools/1.12
module load picard/2.23.4

### The script:  
# 1. Filepaths:
KT_DIR=/proj/snic2020-6-155

# 2. Set array flag to sample names:
samplesheet=$KT_DIR/snailIDs.txt
SnailID=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`

# 2. Index bam
samtools index $KT_DIR/bam/${SnailID}_sorted.bam
# Print flagstats
sed -n '1p;4p;9p' $KT_DIR/bam/${SnailID}.flagstat

# 3. Mark duplicates
java -Xmx16G -jar $PICARD_ROOT/picard.jar MarkDuplicates \
        I=$KT_DIR/bam/${SnailID}_sorted.bam \
        O=$KT_DIR/bam/${SnailID}_markDups.bam \
        M=$KT_DIR/bam/${SnailID}_mark_duplicate_stats.txt \
        TMP_DIR=$TMPDIR

# Print number of duplicates
samtools flagstats $KT_DIR/bam/${SnailID}_markDups.bam | sed -n '4p'

# 4. Filter BAM
### Keeping properly paired alignments (-f 2)
### Remove duplicated alignments (-F 1024)
samtools view -b -f 2 -F 1024 $KT_DIR/bam/${SnailID}_markDups.bam > $KT_DIR/bam/${SnailID}_nodups_paired.bam

# 5. Index BAMs
samtools index $KT_DIR/bam/${SnailID}_nodups_paired.bam

# 6. Create flagstats
samtools flagstats $KT_DIR/bam/${SnailID}_nodups_paired.bam > $KT_DIR/bam/${SnailID}_nodups_paired.flagstat
head -n 1 $KT_DIR/bam/${SnailID}_nodups_paired.flagstat
