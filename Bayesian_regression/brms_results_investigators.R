################################## Investigating brms model fits ############################
### Viewing the results from brms fits run on a cluster. Each model is summarised as plots of 
### the posterior distribution, track plots and summary table. Finally, the tested model is 
### compared to a null model using the difference in LOOIC.

### James Reeve - University of Gothenburg
### 2024-01-22

### Preparation
rm(list = ls())
dev.off()

### Packages
library(brms)
library(ggplot2)

### Filepath
PATH <- "/path/to/brms_results/"

#### Parameters ####

### Specify which analysis
Analysis <- ...
# Values: "Pheno_Env", "Pheno_Inv", "Inv_Env", or "Pheno_InvEnv

### Specify the focal response variable
Variable <- ...
# Values for Pheno_Env & Pheno_InvEnv: "PC1.pheno", "PC2.pheno", or "PC3.pheno"
# Values for Pheno_Inv: "PC1", "PC2", or "PC3"
# Values for Inv_Env: "LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", "LGC6.1.2b, "LGC7.1", "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3" or "LGC17.1"

### Specify the habitat for environmental variables
Habitat <- ...
# Values: "rock", "boulder", "mud-sand" or "hab.only"
# Note: for the Pheno_Inv model set this to "Inv"


#### View the model fit results ####

### Load fits
if(Analysis == "Pheno_InvEnv"){
  # Final model fits (Pheno~Inv+Env models)
  load(paste0(PATH, Analysis, "/brms_fit.", Variable, "_InvEnv_", Habitat)) # Full model
  load(paste0(PATH, Analysis, "/brms_fit.", Variable, "_InvEnv_", Habitat, "_null")) # Null model
  
} else {
  # First round of model fits
  load(paste0(PATH, Analysis, "/brms_fit.", Variable, "_", Habitat)) # Full model
  load(paste0(PATH, Analysis, "/brms_fit.", Variable, "_", Habitat, "_null")) # Null model
  
}

### Posterior distribution and trackplots
plot(mod)

### Compare looic scores
mod$criteria
mod_null$criteria
loo_compare(mod_null, mod, criterion = "loo")
