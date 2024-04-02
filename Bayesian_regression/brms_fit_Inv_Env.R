###################################### Inversion - environment models ##################################
### Running Bayesian logistic models to determin how envrionmental variables from 196 Littorina saxatilis 
### sampling sites influence inversion frequencies. Four models are tested for each inversion. The among 
### habitat model, and three within habitat models. A seprate null model with just the site effect is run 
### to determine the crediability of each tested model. This job was run on the server 'Rackham', part of 
### the National Academic Infrastructure for Supercomputing in Sweden.

### James Reeve - University of Gothenburg
### 2023-06-05

### Preparation
options(stringsAsFactors = FALSE)
options(repos=structure(c(CRAN="https://ftp.acc.umu.se/mirror/CRAN/")))

### Filepaths
PATH <- "/path/to/data/"

### Packages
library(tidyverse) # lib.loc may need to be set if you create your own R package library
library(brms)

### Parameters
Params <- commandArgs(trailingOnly = TRUE) # This draws parameters from a wraper bash script
Inv <- Params[1]
# Possible values: "LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", "LGC6.1.2b", "LGC7.1", "LGC7.2", "LGC9.1", 
# "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3", and "LGC17.1"
Habitat <- Params[2]
# Possible values: "rock", "boulder", "mud-sand" and "hab.only"
Iter <- Params[3]
Chain <- Params[4]


#### A: Read in the data ####

### Envrionmental PCA
# Rocky shores
RS <- read.csv(paste0(PATH, "envPCA/Env_RS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
RS$Habitat <- "rock"
# Boulder shores
BS <- read.csv(paste0(PATH, "envPCA/Env_BS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
BS$Habitat <- "boulder"
# Muddy or sandy shores
MSS <- read.csv(paste0(PATH, "envPCA/Env_MSS_PC1-5.csv"), col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5"))
MSS$Habitat <- "mud-sand"
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

### Dissection records to get sex
Diss <- read.csv(paste0(PATH, "KT_dissection_log.csv"))



#### B: Data wrangeling ####

### 1: Separate SNP and inversion genotypes
GT_inv <- GT[,c(1, grep("LGC", colnames(GT)))]
GT_snp <- GT[,c(1, grep("Contig", colnames(GT)))]

# Add LGC14.1 to GT_inv using the values in Contig5977_84235
GT_inv$LGC14.1 <- GT$Contig5977_84235
GT_inv$LGC14.1[GT_inv$LGC14.1 == 0] <- "RR"
GT_inv$LGC14.1[GT_inv$LGC14.1 == 1] <- "RA"
GT_inv$LGC14.1[GT_inv$LGC14.1 == 2] <- "AA"

# Add trails parameter to LGC6.1-2. This adjust brms to ignore the 'R' arragement
# when running LGC6.1-2b.
GT_inv$LGC6.1.2b_trials <- NA
GT_inv[GT_inv$LGC6.1.2 == "RR", "LGC6.1.2b_trials"] <- 0
GT_inv[GT_inv$LGC6.1.2 == "RA", "LGC6.1.2b_trials"] <- 1
GT_inv[GT_inv$LGC6.1.2 == "RB", "LGC6.1.2b_trials"] <- 1
GT_inv[GT_inv$LGC6.1.2 == "AB", "LGC6.1.2b_trials"] <- 2
GT_inv[GT_inv$LGC6.1.2 == "AA", "LGC6.1.2b_trials"] <- 2
GT_inv[GT_inv$LGC6.1.2 == "BB", "LGC6.1.2b_trials"] <- 2

# Further sample 'GT_snp' to just genetic background
SNPs <- SNP_info %>%
  filter(inv == "background") %>%
  summarise("SNP" = paste(CHROM, POS, sep = "_")) %>%
  unlist()
GT_snp <- GT_snp[,SNPs]

rownames(GT_snp) <- GT$snail

# Clean-up files
rm(GT, SNPs, SNP_info)


### 2: Combine enviromental data and inversion frequencies
# Add site ID to inversions
GT_inv$site <- substr(GT_inv$snail, 1, 5)
# Merge data by site
InvData <- inner_join(GT_inv, Env, by = c("site" = "ID"))


### 3: Add sex and maturity to 'InvData'
InvData <- Diss %>% select(snail_ID, sex, adult) %>%
  inner_join(InvData, by = c("snail_ID" = "snail"))

# Rescore juveniles to females
InvData[InvData$sex == "J", "sex"] <- "F"

# Rescore as c("adult", "juvenile")
InvData$adult <- ifelse(InvData$adult, "adult", "juvenile")


### 4: Get genetic covaraince from SNP genotypes
# Trnspose genotype table and divide by 2
GT_snp <- t(GT_snp / 2)

# Covariances
GT_cov <- cov(GT_snp)

# Add a tiny value to the within-group covariance (i.e. variance)
# to avoid a 'within-group covariance' error in brms
diag(GT_cov) <- diag(GT_cov) + min(GT_cov)/100


### 4: Change inversion frequnecies to be numeric response
# Select focal inversion and name it 'Inv'
if(Inv == "LGC6.1.2b"){
        colnames(InvData)[which(colnames(InvData) == "LGC6.1.2")] <- "Inv"
} else {
        colnames(InvData)[which(colnames(InvData) == Inv)] <- "Inv"
}

# Rescore inversion frequencies
if(Inv == "LGC6.1.2b"){
        # Change so R|A = 0 and B = 1
        InvData$Inv[InvData$Inv == "RR"] <- 0
        InvData$Inv[InvData$Inv == "RA"] <- 0
        InvData$Inv[InvData$Inv == "RB"] <- 1
        InvData$Inv[InvData$Inv == "AA"] <- 0
        InvData$Inv[InvData$Inv == "AB"] <- 1
        InvData$Inv[InvData$Inv == "BB"] <- 2
} else {
        # Change so R = 0 and A|B = 1
        InvData$Inv[InvData$Inv == "RR"] <- 0
        InvData$Inv[InvData$Inv == "RA"] <- 1
        InvData$Inv[InvData$Inv == "RB"] <- 1
        InvData$Inv[InvData$Inv == "AA"] <- 2
        InvData$Inv[InvData$Inv == "AB"] <- 2
        InvData$Inv[InvData$Inv == "BB"] <- 2
}
InvData$Inv <- as.numeric(InvData$Inv)



#### C: select key parameters ####

### An additonal filtering step to select parameters that were specified at the top of the script

### 1: Remove individuals not in covariance matrix
InvData <- InvData[InvData$snail_ID %in% rownames(GT_cov),]

### 2: Subset data based on habitat score
if(Habitat %in% unique(InvData$Habitat)){
  InvData <- InvData[InvData$Habitat == Habitat, ]
} else { # Drop the envrionmental PCA scores
  InvData <- InvData[, grep(".env", colnames(InvData), invert = TRUE)]
}



#### D: brms model ####

### 1: Null model

if(Inv == "LGC6.1.2b"){
        # Change so N trials RR = 0, RA = 1 & RB = 1
        mod_null <- brm(Inv | trials(LGC6.1.2b_trials) ~ sex + (1 | gr(snail_ID, cov = GD)),
                        data = InvData, data2 = list(GD = GT_cov),
                        family = binomial(link = "logit"),
                        iter = Iter, chains = Chain, cores = Chain)
} else {
        # Otherwise N trials = 2
        mod_null <- brm(Inv | trials(2) ~ sex + adult + (1 | gr(snail_ID, cov = GD)),
                        data = InvData, data2 = list(GD = GT_cov),
                        family = binomial(link = "logit"),
                        iter = Iter, chains = Chain, cores = Chain)
}

mod_null <- add_criterion(mod_null, criterion = "loo")

# Save model fit
save(mod_null, file = paste0(PATH, "brms_results/brms_fit.", Inv, "_", Habitat, "_null"))


### 2: Tested model

if(Habitat == "hab.only"){
        if(Inv == "LGC6.1.2b"){
                # hab.only LGC6.1/2b
                mod <- brm(Inv | trials(LGC6.1.2b_trials) ~ Habitat + sex + adult + (1|gr(snail_ID, cov = GD)),
                        data = InvData, data2 = list(GD = GT_cov),
                        family = binomial(link = "logit"),
                        iter = Iter, chains = Chain, cores = Chain)
        } else {
                # hab.only all other inversions
                mod <- brm(Inv | trials(2) ~ Habitat + sex + adult + (1|gr(snail_ID, cov = GD)),
                        data = InvData, data2 = list(GD = GT_cov),
                        family = binomial(link = "logit"),
                        iter = Iter, chains = Chain, cores = Chain)
        }
} else {
        if(Inv == "LGC6.1.2b"){
                # RS|BS|MSS LGC6.1/2b
                mod <- brm(Inv | trials(LGC6.1.2b_trials) ~ sex + adult +
                                PC1 + PC2 + PC3 + (1|gr(snail_ID, cov = GD)),
                        data = InvData, data2 = list(GD = GT_cov),
                        family = binomial(link = "logit"),
                        iter = Iter, chains = Chain, cores = Chain)
        } else {
                # RS|BS|MSS all other inversions
                mod <- brm(Inv | trials(2) ~ sex + adult +
                                PC1 + PC2 + PC3 + (1|gr(snail_ID, cov = GD)),
                        data = InvData, data2 = list(GD = GT_cov),
                        family = binomial(link = "logit"),
                        iter = Iter, chains = Chain, cores = Chain)
        }
}

mod <- add_criterion(mod, criterion = "loo")
save(mod, file = paste0(PATH, "brms_results/brms_fit.", Inv, "_", Habitat))
