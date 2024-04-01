!/bin/bash

### Short script to extract genotype scores from the filtered VCF from the WGS Koster Data.
### This job was run on the server 'Rackham', part of the National Academic Infrastructure 
### for Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 29/08/2023

#SBATCH -A snic2022-22-823
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J KT_SNPpanel_extract_GT_from_VCF
#SBATCH --array=1-17
#SBATCH --output=KT_logfiles/KT_extract_GT_from_LG%a_VCF.out
#SBATCH --error=KT_logfiles/KT_extract_GT_from_LG%a_VCF.err

### Load modules
module load bioinfo-tools
module load bcftools/1.14

### VCF directory path
VCF=/proj/snic2020-6-155/vcf

### List of linkage groups
LG=`echo "LG"$SLURM_ARRAY_TASK_ID`

### Extract SNPs from BGI data
bcftools query -f '%CHROM %POS [ %GT]\n' $VCF/KT_${LG}_filtered_MAF003.vcf.gz \
        > /proj/snic2020-6-155/genotypes/KT_${LG}_genotypes_BGI.txt
