#!/bin/bash

### Merge vcf.gz files that were called seperately for different linkage groups
### on the Littorina saxatilis genome. Note, this file will still be missing 
### any SNP that don't fall on the linkage map. This job was run on the server 
### 'Rackham', part of the National Academic Infrastructure for Supercomputing 
### in Sweden.

### James Reeve - Göteborgs Universitet
### 21/06/2022

### Job parameters
#SBATCH -A snic2021-22-384
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J concat_vcfs
#SBATCH --output=KT_logfiles/concat_vcf.out

### Required software
module load bioinfo-tools
module load bcftools/1.12

### The script:
# 1. Filepaths:
VCF=/proj/snic2020-6-155/vcf

# 2. Make list of vcfs
ls $VCF/*filtered.MAF003.vcf.gz > $TMPDIR/vcf_list.txt

# 3. Concatenate VCFs
bcftools concat -a -d all -Oz -o $VCF/KT_filtered.vcf.gz \
       $VCF/KT_LG1_filtered_MAF003.vcf.gz \
       $VCF/KT_LG2_filtered_MAF003.vcf.gz \
       $VCF/KT_LG3_filtered_MAF003.vcf.gz \
       $VCF/KT_LG4_filtered_MAF003.vcf.gz \
       $VCF/KT_LG5_filtered_MAF003.vcf.gz \
       $VCF/KT_LG6_filtered_MAF003.vcf.gz \
       $VCF/KT_LG7_filtered_MAF003.vcf.gz \
       $VCF/KT_LG8_filtered_MAF003.vcf.gz \
       $VCF/KT_LG9_filtered_MAF003.vcf.gz \
       $VCF/KT_LG10_filtered_MAF003.vcf.gz \
       $VCF/KT_LG11_filtered_MAF003.vcf.gz \
       $VCF/KT_LG12_filtered_MAF003.vcf.gz \
       $VCF/KT_LG13_filtered_MAF003.vcf.gz \
       $VCF/KT_LG14_filtered_MAF003.vcf.gz \
       $VCF/KT_LG15_filtered_MAF003.vcf.gz \
       $VCF/KT_LG16_filtered_MAF003.vcf.gz \
       $VCF/KT_LG17_filtered_MAF003.vcf.gz
