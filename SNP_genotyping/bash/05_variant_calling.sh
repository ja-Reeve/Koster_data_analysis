#!/bin/bash

### SNP calling for 99 BAM files created from the WGS of Littorina saxatilis sampled
### from Western Sweden. BAMs were generated with BWA-MEM, then sorted and filtered
### to remove duplicates and imporperlly paired reads. This script converts these
### BAMs together into an mpileup file which is then variant called with bcftools.
### This job was run on the server 'Rackham', part of the National Academic 
### Infrastructure for Supercomputing in Sweden.

### James Reeve - Göteborgs Universitet
### 13/10/2021

### Job parameters
#SBATCH -A snic2021-5-142
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 3-00:00:00
#SBATCH --array=1-17
#SBATCH -J KT_var_call_LG%a
#SBATCH --output=KT_logfiles/KT_var_call_LG%a.out

### Required software
module load bioinfo-tools
module load samtools/1.12
module load bcftools/1.12

### The script:  
# 1. Filepaths:
KT_DIR=/proj/snic2020-6-155/James
KT_REF=$KT_DIR/ref_genome

# 2. Create list of bams:
ls $KT_DIR/bam/*_nodups_paired.bam > $KT_DIR/bam_list.txt

# 3. Variant call:
bcftools mpileup -b $KT_DIR/bam_list.txt \
        -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR \
        -f $KT_REF/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta \
        -R $KT_REF/LG${SLURM_ARRAY_TASK_ID}.bed | \
        bcftools call -mOz -v --threads 4 -f GQ,GP -o $KT_DIR/vcf/KT_${SLURM_ARRAY_TASK_ID}.vcf.gz
