#!/bin/bash

### Find the identity (i.e. A/C/T/G) of the reference and alternative alleles
### in 1000 SNPs which will be used in a SNP panel. This job was run on the 
### server 'Rackham', part of the National Academic Infrastructure for 
### Supercomputing in Sweden.

### James Reeve - Göteborgs Universitet
### 31/01/2022

### Job parameters
#SBATCH -A snic2021-22-384
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0:30:00
#SBATCH --array=1-17
#SBATCH -J KT_Ref_and_Alt_LG%a
#SBATCH --output=KT_logfiles/Ref_and_Alt_LG%a.out

### Required software
module load bioinfo-tools
module load vcftools/0.1.16
module load bcftools/1.14

### The script:  
# Filepath to data
KTDIR=/proj/snic2020-6-155

# Extract contig and position from BED file
# This BED file was made with the R script 'SNP_panel_selection.R'
awk 'BEGIN {OFS = "\t"} {print $1,$3}' $KTDIR/Koster_SNP_panel_YYYY-MM-DD.bed > $TMPDIR/SNP_panel.bed

# Find VCF entry for each SNP marker
# and extract CHROM, POS, REF & ALT
vcftools --gzvcf $KTDIR/vcf/KT_LG${SLURM_ARRAY_TASK_ID}_filtered_MAF003.vcf.gz \
  --positions $TMPDIR/SNP_panel.bed --recode -c > $TMPDIR/KT_SNPpanel_LG${SLURM_ARRAY_TASK_ID}.vcf
bcftools query -f '%CHROM %POS %REF %ALT\n' $TMPDIR/KT_SNPpanel_LG${SLURM_ARRAY_TASK_ID}.vcf > \
  $KTDIR/KT_SNPpanel_RefAlt_LG${SLURM_ARRAY_TASK_ID}.txt
