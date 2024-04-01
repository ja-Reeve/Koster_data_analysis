#!/bin/bash

### Filter SNPs from VCF files generated from WGS of 99 Littorina saxatilis sampled 
### Western Sweden. Variants were called with a bcftools-mpileup pipeline. Filtering 
### was conducted by Eva Koch. This is a transcription of her description. This job 
### was run on the server 'Rackham', part of the National Academic Infrastructure for 
### Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 27/09/2022

### Job parameters
#SBATCH -A snic2022-22-823
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:30:00
#SBATCH --array=1-17
#SBATCH -J KT_filter_vcfs
#SBATCH --output=KT_logfiles/Filter_vcf_LG%a.out

### Required software
module load bioinfo-tools
module load vcftools/0.1.16
module load bcftools/1.14

### Script:
# Filepath
KTDIR=/proj/snic2020-6-155

# Index VCF
bcftools sort -Ob $KTDIR/vcf/KT_LG${SLURM_ARRAY_TASK_ID}.vcf.gz | \
bcftools index -o $KTDIR/vcf/KT_LG${SLURM_ARRAY_TASK_ID}.csi

# Filter vcf files
vcftools --gzvcf $KTDIR/vcf/KT_LG${SLURM_ARRAY_TASK_ID}.vcf.gz \
        --min-alleles 2 --max-alleles 2 --minGQ 20 \
        --min-meanDP 5.0 --max-meanDP 20.0 --minQ 20 \
        --max-missing 0.2 --remove-indels \
        --mac 2 --maf 0.03 --remove-filtered-all \
        --recode -c > $KTDIR/vcf/KT_LG${SLURM_ARRAY_TASK_ID}_filtered_MAF003.vcf

