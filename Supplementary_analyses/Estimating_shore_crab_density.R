#### Crab sampling trial analysis #####
### This script analyses crab observations taken in August 2020, to determine if crab danisity
### changes in different Littorina saxatilis habitat types. These habitats are grouped into 
### 6 categories based on the avaliability of carb habitat;
###   i.  H1: highly exposed; no seaweed
###   ii. H2: moderately exposed; Fucus vesiculosus without bladders
###   iii.H3: semi-sheltered; Fucus vesiculosus with bladders but no Ascophyllum
###   iv. H4: sheltered; Fucus vesiculosus with bladders and Ascophyllum
###   v.  H5: highly sheltered; mud flat with small patches of seaweed, often grass on shore
###   vi. H6: sandy shore; small patches of seaweed, often no Littorina saxatilis

### James Reeve - University of Gothenburg
### 24/08/2020

#### 1. Upload the data ####
## Preparation
setwd("C:/Users/James/")
rm(list = ls())

## Packages
library("ggplot2")
library("dplyr")
library("lme4")

## Upload data
crab <- read.csv("file/path/crabbing_trial_data_08-2020.csv", header = T)



#### 2. Inspect the data ####
## Means and totals
sum(crab$N_large) # Total number of large crabs
sum(crab$N_small) # Total number of small crabs
mean(crab$N_large) # Mean large crabs per site
mean(crab$N_small) # Mean small crabs per site

## Plot number per site
ps1 <- ggplot(crab)+
  geom_col(aes(Site, N_small))+
  theme_bw()

ps2 <- ggplot(crab)+
  geom_col(aes(Site, N_large))+
  theme_bw()

ps3 <- ggplot(crab)+
  geom_col(aes(Site, N_small+N_large))+
  theme_bw()

## Plot number per habitat type
ps4 <- ggplot(crab)+
  geom_col(aes(Habitat_category, N_small))+
  theme_bw()

ps5 <- ggplot(crab)+
  geom_col(aes(Habitat_category, N_large))+
  theme_bw()

ps6 <- ggplot(crab)+
  geom_col(aes(Habitat_category, N_small+N_large))+
  theme_bw()

## Plot depth vs number
ps7 <- ggplot(crab)+
  geom_point(aes(Depth, N_small, colour = Habitat_category))+
  theme_bw()

ps8 <- ggplot(crab)+
  geom_point(aes(Depth, N_large, colour = Habitat_category))+
  theme_bw()

ps9 <- ggplot(crab)+
  geom_point(aes(Depth, N_small+N_large, colour = Habitat_category))+
  theme_bw()


## Plot distance from transect vs number
# This is an error check, ideally distance shouldn't matter
ggplot(crab)+
  geom_boxplot(aes(x = as.character(Position_m), y = N_small))+
  facet_grid(cols = vars(crab$Habitat_category))+
  theme_bw()

ggplot(crab)+
  geom_boxplot(aes(x = as.character(Position_m), y = N_large))+
  facet_grid(cols = vars(crab$Habitat_category))+
  theme_bw()

ggplot(crab)+
  geom_boxplot(aes(x = as.character(Position_m), y = N_small+N_large))+
  facet_grid(cols = vars(crab$Habitat_category))+
  theme_bw()



#### 3. Analyse the data ####
## Chi-squared test
# Start out simple by looking for overall significant difference between habitats
chisq.test(x = crab$N_small, y = crab$Habitat_category)
chisq.test(x = crab$N_large, y = crab$Habitat_category)
chisq.test(x = crab$Habitat_category, y = crab$N_small+ crab$N_large)
### Conclusion: looks like there's a significant difference in the number of small crabs, but not large crabs.

## Poisson regression (generalized linear model)
# Simplest model number of crabs is a product of habitat
m1 <- glm((N_small+N_large) ~ Habitat_category + Depth, family = "poisson", data = crab)
summary(m1)

# Model 2: add trial and depth as covariates
m2 <- glm((N_small+N_large) ~ Habitat_category + Trial + Depth, family = "poisson", data = crab)
summary(m2)

# Model 3: add Trial and Time as random effects - using lme4 package
m3 <- glmer((N_small+N_large) ~ Habitat_category + (1|Trial) + Depth + (1|Site), 
            family = "poisson", data = crab)
summary(m3)

### Conclusion: no effect of habitat and a significant effect of depth. But, it's not over. The model might be effected
### by the large number of 0s. I will try again after reordering the variables.



#### 4. Correlations between PCs in main data and crab densities ####

### Access PC scores

# Phenotypes
Pheno <- read.csv("file/path/Phenotype/Pheno_PC1-8.v3.csv")
Pheno$Site <- substr(Pheno$X, 1, 5)

# Environment
Env.RS <- read.csv("file/path/Envrionment/Env_perSnail_RS_PC1-5.csv")
Env.BS <- read.csv("file/path/Environment/Env_perSnail_BS_PC1-5.csv")
Env.MSS <- read.csv("file/path/Environment/Env_perSnail_MSS_PC1-5.csv")

Env.RS$Site <- substr(Env.RS$X, 1, 5)
Env.BS$Site <- substr(Env.BS$X, 1, 5)
Env.MSS$Site <- substr(Env.MSS$X, 1, 5)

### Access shellLength
Diss <- read.csv("file/path/Phenotype/KT_ShellShaper_parameters.v2.csv")
# Subset to only snailID and shellLength
Diss <- Diss[, c("snailID", "shellLength")]
# Drop the "_ac" suffix from snail IDs
Diss$snailID <- substr(Diss$snailID, 1, 6)
Diss$Site <- substr(Diss$snailID, 1, 5)

### Comparing phenotypes and crab densities
tmp <- crab %>% group_by(Site) %>% 
  summarise(Habitat_category, "N" = sum(N_small + N_large)) %>% 
  unique()

# Site average PCs
tmp2 <- Pheno %>% group_by(Site) %>%
  summarise("PC1" = mean(PC1),
            "PC2" = mean(PC2)) %>%
  right_join(., tmp, by = join_by(Site)) %>%
  filter(Site != "KT027") # Remove site missing snails

# Add shellLength
tmp3 <- Diss %>% group_by(Site) %>%
  summarise("shellLength" = mean(shellLength)) %>%
  right_join(., tmp2, by = join_by(Site))


### Spearman's correlation
cor.test(tmp3$PC1, tmp3$N, method = "spearman")
cor.test(tmp3$PC2, tmp3$N, method = "spearman")
cor.test(tmp3$shellLength, tmp3$N, method = "spearman")



### Plot
p1 <- ggplot(tmp3, aes(x = N, y = PC1))+
  geom_point(aes(colour = Habitat_category))+
  geom_smooth(method = lm, se = FALSE)+
  theme_classic()
  
p2 <- ggplot(tmp3, aes(x = N, y = PC2))+
  geom_point(aes(colour = Habitat_category))+
  geom_smooth(method = lm, se = FALSE)+
  theme_classic()

p3 <- ggplot(tmp3, aes(x = N, y = shellLength))+
  geom_point(aes(colour = Habitat_category))+
  geom_smooth(method = lm, se = FALSE)+
  theme_classic()

ggpubr::ggarrange(p1, p2, p3, ncol = 1)


### Merge env and crab
# Rocky shores
tmp.RS <- Env.RS %>% group_by(Site) %>%
  summarise("Dim1" = mean(dim.1),
            "Dim2" = mean(dim.2),
            "Dim3" = mean(dim.3)) %>%
  inner_join(., tmp, by = join_by(Site)) 

cor.test(tmp.RS$Dim1, tmp.RS$N, method = "spearman")
cor.test(tmp.RS$Dim2, tmp.RS$N, method = "spearman")
cor.test(tmp.RS$Dim3, tmp.RS$N, method = "spearman")

# Boulder fields
tmp.BS <- Env.BS %>% group_by(Site) %>%
  summarise("Dim1" = mean(dim.1),
            "Dim2" = mean(dim.2),
            "Dim3" = mean(dim.3)) %>%
  inner_join(., tmp, by = join_by(Site)) 

cor.test(tmp.BS$Dim1, tmp.BS$N, method = "spearman")
cor.test(tmp.BS$Dim2, tmp.BS$N, method = "spearman")
cor.test(tmp.BS$Dim3, tmp.BS$N, method = "spearman")

# Muddy-sandy shores
tmp.MSS <- Env.MSS %>% group_by(Site) %>%
  summarise("Dim1" = mean(dim.1),
            "Dim2" = mean(dim.2),
            "Dim3" = mean(dim.3)) %>%
  inner_join(., tmp, by = join_by(Site)) 

cor.test(tmp.MSS$Dim1, tmp.MSS$N, method = "spearman")
cor.test(tmp.MSS$Dim2, tmp.MSS$N, method = "spearman")
cor.test(tmp.MSS$Dim3, tmp.MSS$N, method = "spearman")

