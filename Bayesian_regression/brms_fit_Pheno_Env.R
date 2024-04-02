###################################### Pehnotype - environment models ###################################
### Running Bayesian logistic models to determin how envrionmental variables from 196 Littorina saxatilis 
### sampling sites influence phenotypic traits. Two models are tested for each pehnotypic trait. The among 
### habitat model, and three within habitat models. A seprate null model with just the site effect is run 
### to determine the crediability of each tested model. This job was run on the server 'Rackham', part of 
### the National Academic Infrastructure for Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 2023-07-03

### Preparation
options(stringsAsFactors = FALSE)
options(repos=structure(c(CRAN="https://ftp.acc.umu.se/mirror/CRAN/")))

### Filepaths
PATH <- "/path/to/data/"

### Packages
library(tidyverse) # lib.loc may need to be set if you create your own R package library
library(brms)

### Parameters
Iter <- 40000
Chain <- 4
Params <- commandArgs(trailingOnly = TRUE) # This draws parameters from a wraper bash script
Trait <- Params[1]
# Possible values: "PC1.pheno", "PC2.pheno", and "PC3.pheno"
Habitat <- Params[2]
# Possible values: "rock", "boulder", "mud-sand" and "hab.only"



#### A: Read in the data ####

Pheno <- read.csv(paste0(PATH, "phenoPCA/Pheno_PC1-9.v2.csv"),
                      col.names = c("snail_ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9"))

### Dissection record
Diss <- read.csv(paste0(PATH, "KT_dissection_log.csv"))

# Drone photo measurements
DPh <- read.csv(paste0(PATH,"KT_drone_photos_v2.csv"))


### Envrionmental PCA
# Rocky shores
RS <- read.csv(paste0(PATH, "envPCA/Env_RS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
# Boulder shores
BS <- read.csv(paste0(PATH, "envPCA/Env_BS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
# Muddy or sandy shores
MSS <- read.csv(paste0(PATH, "envPCA/Env_MSS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
# Merge data into single dataframe
Env <- rbind.data.frame(RS, BS, MSS)
Env <- Env[order(Env$ID),]
rm(RS,BS,MSS)

### Genotypes
GT <- read.csv(paste0(PATH, "Koster_SNP_imputed_20221031.csv"), header = TRUE)
GT <- GT[grep("KT", GT$snail), ] # Filter to just Koster data
# Drop duplicates
GT <- GT[!(duplicated(GT$snail)),]

### Get SNP info
SNP_info <- read.csv(paste0(PATH, "Koster_SNP_info_FINAL.csv"))



#### B: Data wrangeling ####

### 1: Filter to SNP genotypes
GT_snp <- GT[,c(1, grep("Contig", colnames(GT)))]

# Further sample 'GT_snp' to just genetic background markers
SNPs <- SNP_info %>%
  filter(inv == "background") %>%
  summarise("SNP" = paste(CHROM, POS, sep = "_")) %>%
  unlist()
GT_snp <- GT_snp[,SNPs]

rownames(GT_snp) <- GT$snail

# Clean-up files
rm(GT, SNPs, SNP_info)


### 2: Add sex to 'Pheno'
Pheno <- Diss %>% select(snail_ID, site, sex, adult) %>%
  right_join(Pheno, by = "snail_ID")

# Rescore juveniles to females
Pheno[Pheno$sex == "J", "sex"] <- "F"

# Rescore as c("adult", "juvenile")
Pheno$adult <- ifelse(Pheno$adult, "adult", "juvenile")


### 3: Add habitat to 'Pheno'
Pheno <- DPh %>% select(Site, Habitat) %>%
  right_join(Pheno, by = c("Site" = "site"))


### 4: Add envrionmental varaibles
PhenoDat <- inner_join(Pheno, Env, by = c("Site" = "ID"), suffix = c(".pheno", ".env"))
# Filter to snails in genotypes
PhenoDat <- PhenoDat[PhenoDat$snail_ID %in% rownames(GT_snp),]


### 5: Subset based on habitat score
if(Habitat %in% unique(PhenoDat$Habitat)){
  PhenoDat <- PhenoDat[PhenoDat$Habitat == Habitat, ]
} else { # Drop the environmental PCA scores
  PhenoDat <- PhenoDat[, grep(".env", colnames(PhenoDat), invert = TRUE)]
}


### 6: Get genetic covaraince from SNP genotypes

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

### Rename phenotype of interest
colnames(PhenoDat)[which(colnames(PhenoDat) == Trait)] <- "Pheno"



#### D: brms model ####

### 1: Null model

mod_null <- brm(Pheno ~ sex + adult + (1|gr(snail_ID, cov = GD)),
                data = PhenoDat, data2 = list(GD = GT_cov),
                iter = Iter, chains = Chain, cores = 4L)
mod_null <- add_criterion(mod_null, criterion = "loo")
# Save model fit
save(mod_null, file = paste0(PATH, "brms_results/brms_fit.", Trait, "_", Habitat, "_null"))


### 2: Tested model

if(Habitat == "hab.only"){
  mod <- brm(Pheno ~ Habitat + sex + adult + (1|gr(snail_ID, cov = GD)),
             data = PhenoDat, data2 = list(GD = GT_cov),
             iter = Iter, chains = Chain, cores = 4L)
} else {
  mod <- brm(Pheno ~ sex + adult + PC1.env + PC2.env +
               PC3.env + (1|gr(snail_ID, cov = GD)),
             data = PhenoDat, data2 = list(GD = GT_cov),
             iter = Iter, chains = Chain, cores = 4L)
}

mod <- add_criterion(mod, criterion = "loo")
save(mod, file = paste0(PATH, "brms_results/brms_fit.", Trait, "_", Habitat))
