##################################### Pehnotype - inversion models #########################################
### Running Bayesian logistic models to determin how the inversion frequencies from 1488 Littorina saxatilis 
### samples influence phenotypic traits. One model is tested for each pehnotypic trait. This job was run on 
### the server 'Rackham', part of the National Academic Infrastructure for Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 2023-06-23

### Preparation
options(stringsAsFactors = FALSE)
options(repos=structure(c(CRAN="https://ftp.acc.umu.se/mirror/CRAN/")))

### Filepaths
PATH <- "/path/to/data/"

### Packages
library(tidyverse) # lib.loc may need to be set if you create your own R package library
library(brms)

### Parameters
Iter <- 20000
Chain <- 4
Params <- commandArgs(trailingOnly = TRUE) # This draws parameters from a wraper bash script
Trait <- Params[1]
Inversion <- Params[2]



#### A: Read in the data ####

Pheno <- read.csv(paste0(PATH, "phenoPCA/Pheno_PC1-9.v2.csv"),
                      col.names = c("snail_ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9"))

### Dissection record
Diss <- read.csv(paste0(PATH, "KT_dissection_log.csv"))

### Genotypes
GT <- read.csv(paste0(PATH, "Koster_SNP_imputed_20221031.csv"), header = TRUE)
GT <- GT[grep("KT", GT$snail), ] # Filter to just Koster data
# Drop duplicates
GT <- GT[!(duplicated(GT$snail)),]

### Get SNP info
SNP_info <- read.csv(paste0(PATH, "Koster_SNP_info_FINAL.csv"))



#### B: Data wrangeling ####

### 1: Separate SNP and inversion genotypes
GT_inv <- GT[,c(1, grep("LGC", colnames(GT)))]
GT_snp <- GT[,c(1, grep("Contig", colnames(GT)))]

# Add LGC14.1 to GT_inv using the values in Contig5977_84235
GT_inv$LGC14.1 <- GT$Contig5977_84235
GT_inv$LGC14.1[GT_inv$LGC14.1 == 0] <- "RR"
GT_inv$LGC14.1[GT_inv$LGC14.1 == 1] <- "RA"
GT_inv$LGC14.1[GT_inv$LGC14.1 == 2] <- "AA"

# Further sample 'GT_snp' to just genetic background
SNPs <- SNP_info %>%
  filter(inv == "background") %>%
  summarise("SNP" = paste(CHROM, POS, sep = "_")) %>%
  unlist()
GT_snp <- GT_snp[,SNPs]

# Add rownames
rownames(GT_snp) <- GT$snail

# Clean-up files
rm(GT, SNPs, SNP_info)




### 2: Add sex to 'Pheno'
Pheno <- Diss %>% select(snail_ID, site, sex, adult) %>%
  right_join(Pheno, by = c("snail_ID" = "snail_ID"))

# Rescore juveniles to females
Pheno[Pheno$sex == "J", "sex"] <- "F"
# Rescore as c("adult", "juvenile")
Pheno$adult <- ifelse(Pheno$adult, "adult", "juvenile")


### 3: Add inversion frequencies
# Filter to just snails that are in the SNP data
Pheno <- Pheno[Pheno$snail_ID %in% GT_inv$snail,]
# Filter the genotypes
GT_inv <- GT_inv[Pheno$snail_ID %in% GT_inv$snail,]
GT_snp <- GT_snp[Pheno$snail_ID %in% rownames(GT_snp),]
# Merge the data
PhenoDat <- inner_join(Pheno, GT_inv, by = c("snail_ID" = "snail"))


### 4: Get genetic covaraince from SNP genotypes

# Filter to snails in Phenotypic data
GT_snp <- GT_snp[PhenoDat$snail_ID,]

# Transpose and divide by 2
GT_snp <- t(GT_snp / 2)

# Covariances
GT_cov <- cov(GT_snp)

# Add a tiny value to the within-group covariance (i.e. variance)
# to avoid a 'within-group covariance' error in brms
diag(GT_cov) <- diag(GT_cov) + min(GT_cov)/100



#### C: select key parameters ####

### Rename the phenotype of interest - this is needed to automate the scripts
colnames(PhenoDat)[which(colnames(PhenoDat) == Trait)] <- "Pheno"
colnames(PhenoDat)[which(colnames(PhenoDat) == Inversion)] <- "Inv"



#### D: brms model ####

### 1: Null model

mod_null <- brm(Pheno ~ sex + adult + (1|gr(snail_ID, cov = GD)),
                   data = PhenoDat, data2 = list(GD = GT_cov),
                   iter = Iter, chains = Chain, cores = Chain)
mod_null <- add_criterion(mod_null, criterion = "loo")
# Save model fit
save(mod_null, file = paste0(PATH, "brms_results/brms_fit.", Trait, "_", Inversion, "_null"))


### 2: Tested model

mod <- brm(Pheno ~ sex + adult + LGC1.1 + LGC1.2 + LGC2.1 + LGC4.1 + LGC6.1.2 + LGC7.1 +
                        LGC7.2 + LGC9.1 + LGC10.1 + LGC10.2 + LGC11.1 + LGC12.1 + LGC12.2 +
                        LGC12.3 + LGC12.4 + LGC14.1 + LGC14.3 + LGC17.1 +
                        (1|gr(snail_ID, cov = GD)),
                   data = PhenoDat, data2 = list(GD = GT_cov),
                   iter = Iter, chains = Chain, cores = Chain)

mod <- add_criterion(mod, criterion = "loo")
save(mod, file = paste0(PATH, "brms_results/brms_fit.", Trait, "_", Inversion))
