##################### Analysing the effect of snail height on shore to brms models #######################
### The position of a Littorina saxatilis snail up shore is one known factor driving local adaptation. This
### was not included as an envrionmental variable, as the measurements were not consistently taken among 
### sites. This extra analysis was run to see how much of an effect the snail height could have on the brms
### fits. To save computational power, only inversions previously associated with shore height adaptation
### (LGC9.1 and LGC12.2; Morales et al. 2019) were tested. Note: after Morales et al., 2019 was published, 
### LGC12.2 was renamed. The original publication uses the name 'LGC12.1'.

### James Reeve - DTU Aqua
### 2026-07-02

### Preparation
rm(list = ls())
try(dev.off())
setwd("/Users/james")

### Packages
library(brms)
library(ggplot2)
library(ggpubr)

### Filepath
PATH <- "file/path/"

### Phenotype to investigate
Inv <- ...
# Values = {LGC9.1; LGC12.2}

### Model fits (Pheno~Inv+Env models)
load(paste0(PATH, "brms_fits/brms_fit.", Inv, "_rock_v2")) # Full model

### Run new models for sex, maturity and genetic covariance independently
dat <- mod$data
GCov <- mod$data2

### Access per snail-level environmental data
# There are mainly three relevant columns
#     i) snail_ID
#     ii) adj_height = snail height on shore adjusted by annual mean sea level (cm)
#     iii) substrate = category of material snail was found on - note this may be quite imprecise as the snails can move
SnailEnv <- read.csv("file/path/Environmental_data/Env_data_snail_details.v2.csv")

### Merge the data
dat2 <- merge(dat, SnailEnv[,c("snail_ID", "adj_height", "Substrate")], by = "snail_ID")


### Plot relationship
ggplot(dat2)+
  geom_point(aes(x = adj_height, y = Inv))+
  labs(x = "Snail height from sea level (cm)", y = Inv)


### Fit the model with height rank ####
dat2$height_cat <- substr(dat2$snail_ID, 6, 6)

dat2$height_rank <- 0
dat2[dat2$height_cat == "a", "height_rank"] <- 1
dat2[dat2$height_cat == "b", "height_rank"] <- 2
dat2[dat2$height_cat == "c", "height_rank"] <- 3
dat2[dat2$height_cat == "d", "height_rank"] <- 4
dat2[dat2$height_cat == "e", "height_rank"] <- 5
dat2[dat2$height_cat == "f", "height_rank"] <- 6
dat2[dat2$height_cat == "g", "height_rank"] <- 7
dat2[dat2$height_cat == "h", "height_rank"] <- 8


mod_height.rank <- brm(
  formula = update.formula(formula(mod), ~. + height_rank),
  data = dat2, data2 = GCov,
  family = binomial(link = "logit"),
  iter = 40000, chains = 4, cores = 2L)
mod_height.rank <- add_criterion(mod_height.rank, criterion = "loo")  


### Plot posteriors
mcmc_plot(mod_height.rank, type = "areas",
          pars = c("b_sexM", "b_adultjuvenile", 
                   "b_PC1", "b_PC2", "b_PC3", "b_height_rank"))+
  geom_vline(xintercept = 0, col = "red")+
  labs(x = Inv)
# 400 x 525

## Residual variance
res_mod <- VarCorr(mod)$residual$sd[1]^2
res_height.rank <- VarCorr(mod_height.rank)$residual$sd[1]^2
T.var <- var(dat2$Inv)
(res_mod - res_height.rank) / T.var * 100

### R2
R2.height.rank <- bayes_R2(mod_height.rank)
R2.mod <- bayes_R2(mod)

### Save results
save(mod_height.rank, file = paste0(PATH, "brms_fit.", Inv, "_snail_height_rock"))


#### Plot relationship between inversion genotype and shore height ####

### Merge snail shore heights with inversion genotypes
InvHeight <- merge(SnailEnv[,c("site", "snail_ID", "adj_height", "Source")], dat2,
                   by = "snail_ID")

# Create shore position flag from snail ID
InvHeight$height_rank <- factor(substr(InvHeight$snail_ID, 6, 6))


### Plot overall summary
bar1 <- ggplot(InvHeight)+
  geom_bar(aes(x = height_rank, fill = factor(Inv)),
           position = "dodge")+
  labs(x = "", y = Inv)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank())

### Plot by sample
p1 <- ggplot(data = InvHeight, aes(height_rank, adj_height, colour = factor(Inv)))+
  geom_hline(yintercept = 0, lty = 2, colour = "grey40")+
  geom_line(aes(group = site), alpha = 0.9)+
  geom_point(size = 0.5, alpha = 0.9)+
  labs(x = "Shore position", y = "Height on shore (cm)",
       colour = Inv)+
  theme_classic()


ggarrange(bar1, p1, ncol = 1, nrow = 2,
          heights = c(0.4, 1),
          common.legend = TRUE,
          legend = "bottom")
### 370 x 600
