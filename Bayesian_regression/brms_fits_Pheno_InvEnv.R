##################################### Pehnotype - inversion + environment models ##################################
### Running Baeysian regression models to determine how inversions and habitat effect phenotypic traits from 1488 
### Littorina saxatilis. Two models are run for each trait: the among habitat model and three within habitat models. 
### These are each compared to a Null model which excludes the inversions and envrionmental variables. This job was 
### run on the server 'Rackham', part of the National Academic Infrastructure for Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 2024-01-12

### Preparation
options(stringsAsFactors = FALSE)
options(repos=structure(c(CRAN="https://ftp.acc.umu.se/mirror/CRAN/")))

### Packages
library(tidyverse) # lib.loc may need to be set if you create your own R package library
library(brms)

### Filepaths
PATH <- "/path/to/data/"

### Parameters
Params <- commandArgs(trailingOnly = TRUE) # This draws parameters from a wraper bash script
Trait <- Params[1]
# Possible values = "PC1.pheno", "PC2.pheno", and "PC3.pheno"
Habitat <- Params[2]
# Possible values = "rock", "boulder", "mud-sand", and "hab.only"
Iter <- Params[3]
Chain <- Params[4]



#### A: Read in the data ####

### Phenotyp PCs
Pheno <- read.csv(paste0(PATH, "phenoPCA/Pheno_PC1-9.v2.csv"), col.names = c("snail_ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9"))

### Dissection records
Diss <- read.csv(paste0(PATH, "KT_dissection_log.csv"))

### Drone photo measurements
DPh <- read.csv(paste0(PATH,"KT_drone_photos_v2.csv"))

### Envrionmental PCs
if(Habitat == "hab.only"){
        RS <- read.csv(paste0(PATH, "envPCA/Env_RS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
        BS <- read.csv(paste0(PATH, "envPCA/Env_BS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
        MSS <- read.csv(paste0(PATH, "envPCA/Env_MSS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
# Merge data into single dataframe
        Env <- rbind.data.frame(RS, BS, MSS)
        Env <- Env[order(Env$ID),]
        rm(RS,BS,MSS)}

if(Habitat == "rock") {
        Env <- read.csv(paste0(PATH, "envPCA/Env_RS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))}
if(Habitat == "boulder") {
        Env <- read.csv(paste0(PATH, "envPCA/Env_BS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))}
if(Habitat == "mud-sand") {
        Env <- read.csv(paste0(PATH, "envPCA/Env_MSS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))}

### Genotypes
GT <- read.csv(paste0(PATH, "Koster_SNP_imputed_20221031.csv"), header = TRUE)
GT <- GT[grep("KT", GT$snail), ] # Filter to just Koster data
# Drop duplicates
GT <- GT[!(duplicated(GT$snail)),]

### SNP info
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


### 3: Add habitat to 'Pheno'
Pheno <- DPh %>% select(Site, Habitat) %>%
  right_join(Pheno, by = c("Site" = "site"))


### 4: Add inversion frequencies to 'Pheno'
# Re-score inversions as the count of alternative alleles
tmp <- sapply(GT_inv[,-1],
        function(X){
                X[X == "RR"] <- 0
                X[X == "RA"] <- 1
                X[X == "RB"] <- 1
                X[X == "AA"] <- 2
                X[X == "AB"] <- 2
                X[X == "BB"] <- 2
                return(as.numeric(X))
        })
tmp <- data.frame("snail" = GT_inv$snail, tmp)

# Add frequencies for the A vs. B comparisons in LGC6.1/2
tmp$LGC6.1.2b <- 0
tmp[GT_inv$LGC6.1.2 == "RB", "LGC6.1.2b"] <- 1
tmp[GT_inv$LGC6.1.2 == "AB", "LGC6.1.2b"] <- 1
tmp[GT_inv$LGC6.1.2 == "BB", "LGC6.1.2b"] <- 2

# Add to 'Pheno'
Pheno <- inner_join(Pheno, tmp, by = c("snail_ID" = "snail"))


### 5: Add envrionmental variables
PhenoDat <- inner_join(Pheno, Env, by = c("Site" = "ID"), suffix = c(".pheno", ".env"))
# Filter to snails in genotypes
PhenoDat <- PhenoDat[PhenoDat$snail_ID %in% rownames(GT_snp),]


### 6: Get genetic covaraince from SNP genotypes

# Filter to snails in Phenotypic data
tmp2 <- GT_snp[PhenoDat$snail_ID,]
# Transpose and divide by 2
tmp2 <- t(tmp2 / 2)
# Covariances
GT_cov <- cov(tmp2)

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
save(mod_null, file = paste0(PATH, "brms_results/brms_fit.", Trait, "_InvEnv_", Habitat, "_null"))


### 2: Tested model
if(Habitat == "hab.only"){
  mod <- brm(Pheno ~ Habitat + sex + adult +
                        LGC1.1 + LGC1.2 + LGC2.1 + LGC4.1 + LGC6.1.2 + LGC6.1.2b + LGC7.1 + LGC7.2 +
                        LGC9.1 + LGC10.1 + LGC10.2 + LGC11.1 + LGC12.1 + LGC12.2 + LGC12.3 + LGC12.4 +
                        LGC14.1 + LGC14.3 + LGC17.1 + (1|gr(snail_ID, cov = GD)),
             data = PhenoDat, data2 = list(GD = GT_cov),
             iter = Iter, chains = Chain, cores = 4L)
} else {
  mod <- brm(Pheno ~ sex + adult + PC1.env + PC2.env + PC3.env +
                        LGC1.1 + LGC1.2 + LGC2.1 + LGC4.1 + LGC6.1.2 + LGC6.1.2b + LGC7.1 + LGC7.2 +
                        LGC9.1 + LGC10.1 + LGC10.2 + LGC11.1 + LGC12.1 + LGC12.2 + LGC12.3 + LGC12.4 +
                        LGC14.1 + LGC14.3 + LGC17.1 + (1|gr(snail_ID, cov = GD)),
             data = PhenoDat, data2 = list(GD = GT_cov),
             iter = Iter, chains = Chain, cores = 4L)
}

mod <- add_criterion(mod, criterion = "loo")
save(mod, file = paste0(PATH, "brms_results/brms_fit.", Trait, "_InvEnv_", Habitat))
