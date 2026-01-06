################################### Variance explained by brms models #################################
### Determining the variance explained by each tested model. For models with a phenotype as a response, 
### variance is calculated from the residuals of the brms model. Variance explained is estimated as by
### calclating Bayes R^2 (Gelman et al., 2018) for a series of submodels.

### James Reeve - DTU Aqua
### 2025-11-20

### Preparation
rm(list=ls())
options(stringsAsFactors = FALSE)
#options(repos=structure(c(CRAN="https://ftp.acc.umu.se/mirror/CRAN/")))
#options(mc.cores = parallel::detectCores())

### Packages
### Only run the first line once. It is a waste of time to reinstall pacakges
#install.packages("brms", dependencies = TRUE, INSTALL_opts = "--no-lock")
#install.packages('RcppEigen')
library(dplyr)
library(rstan)
library(brms)

# R-stan options 
rstan_options(auto_write = TRUE)
#rstan_options(threads_per_chain = 2)

### Filepaths
PATH <- "/path/to/brms/model/fits/"

### Parameters
#Params <- commandArgs(trailingOnly = TRUE)
Analysis <- "Pheno_InvEnv"
# Possible values = "Pheno_Inv", "Pheno_Env", "Inv_Env", and "Pheno_InvEnv"

Var <- "PC1.pheno"
# Possible values = "PC1.pheno", "PC2.pheno", and "shellLength"

Hab <- "mud-sand"
# Possible values = "rock", "boulder", "mud-sand", "hab.only", and "Inv"



#### Function to calculate variance explained: ####

Var.Exp <- function(variable, habitat) {
  
  ## A: Get R^2 from submodels
  # Variance expalined by Pheno~Inv model
  load(paste0(PATH, "Pheno_Inv/brms_fit.",
              ifelse(variable != "shellLength", substr(variable, 1, 3), variable),
              "_Inv_v2"))
  R2.Pheno_Inv <- bayes_R2(mod, summary = FALSE)
  
  # Variance explained by Pheno~Env model
  load(paste0(PATH, "Pheno_Env/brms_fit.", variable, "_", habitat, "_v2"))
  R2.Pheno_Env <- bayes_R2(mod, summary = FALSE)
  
  # Variance explained by full model (Pheno~InvEnv)
  # Note: order of operations is important. This should be run after the last two models.
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_", habitat, "_v2"))
  R2.full <- bayes_R2(mod, summary = FALSE)
  
  
  ## B: Run four new submodels
  # Access components from 'mod'
  dat <- mod$data
  GCov <- mod$data2 # genetic covariance matrix
  fml <- mod$formula
  
  # 1) Variance explained without genetic covariance matrix
  mod_noGBG <- brm(update.formula(fml, .~.-(1|gr(snail_ID, cov = GD))),
                   data = dat, data2 = GCov,
                   iter = 4000, chains = 4, cores = 2L, threads = threading(2))
  mod_noGBG <- add_criterion(mod_noGBG, criterion = "loo")
  R2.noGBG <- bayes_R2(mod_noGBG, summary = FALSE)
  rm(mod_noGBG) # Removing to save memory
  
  # 2) Variance explained without co-factors
  mod_noCF <- brm(update.formula(fml, .~.-sex-adult-(1|Date)),
                  data = dat, data2 = GCov,
                  iter = 4000, chains = 4, cores = 2L, threads = threading(2))
  mod_noCF <- add_criterion(mod_noCF, criterion = "loo")
  R2.noCF <- bayes_R2(mod_noCF, summary = FALSE)
  rm(mod_noCF)
  
  # 3) Total variance explained by environment (variance without confounds)
  mod_onlyEnv <- brm(update.formula(fml, .~.-sex-adult-(1|gr(snail_ID, cov = GD))-(1|Date)
                                    -LGC1.1-LGC1.2-LGC2.1-LGC4.1-LGC6.1.2-LGC6.1.2b-LGC7.1-LGC7.2
                                    -LGC9.1-LGC10.1-LGC10.2-LGC11.1-LGC12.1-LGC12.2-LGC12.3-LGC12.4
                                    -LGC14.1-LGC14.3-LGC17.1),
                     data = dat, data2 = GCov,
                     iter = 4000, chains = 4, cores = 2L, threads = threading(2))
  mod_onlyEnv <- add_criterion(mod_onlyEnv, criterion = "loo")
  R2.onlyEnv <- bayes_R2(mod_onlyEnv, summary = FALSE)
  rm(mod_onlyEnv)
  
  # 4) Total variance explained by inversions (variance without confounds)
  # There's a slightly different formula for hab.only models
  if(habitat == "hab.only"){
    mod_onlyInv <- brm(update.formula(fml, .~.-sex-adult-(1|gr(snail_ID, cov = GD))-(1|Date)
                                      -Habitat),
                       data = dat, data2 = GCov,
                       iter = 4000, chains = 4, cores = 2L, threads = threading(2))
  } else {
    mod_onlyInv <- brm(update.formula(fml, .~.-sex-adult-(1|gr(snail_ID, cov = GD))-(1|Date)
                                      -PC1.env-PC2.env-PC3.env),
                       data = dat, data2 = GCov,
                       iter = 4000, chains = 4, cores = 2L, threads = threading(2))
  }
  mod_onlyInv <- add_criterion(mod_onlyInv, criterion = "loo")
  R2.onlyInv <- bayes_R2(mod_onlyInv, summary = FALSE)
  rm(mod_onlyInv)
  
  
  ## C: Calculate remaining variance components
  # Variance explained by inversions - with confounds
  R2.Inv <- R2.full - R2.Pheno_Env
  
  # Variance explained by environment - with confounds
  R2.Env <- sample(R2.full, length(R2.Pheno_Inv), replace = F) - R2.Pheno_Inv
  
  # Variance explained by co-factors (sex, maturity and Date)
  R2.CF <- R2.full - R2.noCF
  
  # Variance explained by factors without genetic background
  R2.GBG <- R2.full - R2.noGBG
  
  # Confounded variance with inversions
  R2.confInv <- R2.onlyInv - R2.Inv
  
  # Confounded variance with environment
  R2.confEnv <- sample(R2.onlyEnv, length(R2.Env), replace = F) - R2.Env
  
  
  ## D: Write data frame
  res <- data.frame("Analysis" = Analysis,
                    "Variable" = variable,
                    "Habitat" = habitat,
                    "R2full" = mean(R2.full),
                    "R2full.SD" = sd(R2.full),
                    "R2onlyInv" = mean(R2.onlyInv),
                    "R2onlyInv.SD" = sd(R2.onlyInv),
                    "R2onlyEnv" = mean(R2.onlyEnv),
                    "R2onlyEnv.SD" = sd(R2.onlyEnv),
                    "R2inv" = mean(R2.Inv),
                    "R2inv.SD" = sd(R2.Inv),
                    "R2env" = mean(R2.Env),
                    "R2env.SD" = sd(R2.Env),
                    "R2cf" = mean(R2.CF),
                    "R2cf.SD" = sd(R2.CF),
                    "R2gbg" = mean(R2.GBG),
                    "R2gbg.SD" = sd(R2.GBG),
                    "R2confInv" = mean(R2.confInv),
                    "R2confInv.SD" = sd(R2.confInv),
                    "R2confEnv" = mean(R2.confEnv),
                    "R2confEnv.SD" = sd(R2.confEnv))
  return(res)
  
}


### Run function
Vexp <- Var.Exp(variable = Var, habitat = Hab)

# Save results
write.csv(Vexp, paste0(PATH, "variance_explained/", Analysis, "_", Var, "_", Hab, "_var_exp.csv"), row.names = FALSE, quote = FALSE)
