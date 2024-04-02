################################### Variance explained by brms models #################################
### Determining the variance explained by each tested model. For models with a phenotype as a response, 
### variance is calculated from the residuals of the brms model. For the models with inversion 
### arrangements as response variables, variance is instead estimated from deviance.

### James Reeve - University of Gothenburg
### 2024-1-22

### Preparation
rm(list = ls())
dev.off()

### Packages
library(brms)
library(ggplot2)

### Filepath
PATH <- "/Volumes/PhD_backup_JR/brms_fits/"

#### Parameters ####
Vars <- c("PC1.pheno", "PC2.pheno", "PC3.pheno")
Habs <- c("rock", "boulder", "mud-sand", "hab.only")
Invs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", "LGC6.1.2b", "LGC7.1", "LGC7.2", "LGC9.1", 
          "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3",
          "LGC17.1")


#### Function to calculate varaince explained: ####

Var.Exp <- function(variable, habitat){
  if(Analysis == "Pheno_InvEnv"){
    # Final model fits (Pheno~Inv+Env models)
    load(paste0(PATH, Analysis, "/brms_fit.", variable, "_InvEnv_", habitat)) # Full model
    load(paste0(PATH, Analysis, "/brms_fit.", variable, "_InvEnv_", habitat, "_null")) # Null model
    
  } else {
    # First round of model fits
    load(paste0(PATH, Analysis, "/brms_fit.", variable, "_", habitat)) # Full model
    load(paste0(PATH, Analysis, "/brms_fit.", variable, "_", habitat, "_null")) # Null model
    
  }
  
  ### LOOIC scores
  looic.mod <- mod$criteria$loo$estimates[3,1]
  looic.null <- mod_null$criteria$loo$estimates[3,1]
  
  ### Delta-LOOIC
  if(looic.mod > looic.null){
    dlooic <- abs(loo_compare(mod, mod_null)[2,1])
  } else {
    dlooic <- loo_compare(mod, mod_null)[2,1]
  }
  dlooic.se <- loo_compare(mod, mod_null)[2,2]
  
  if(Analysis == "Inv_Env"){
    
    ### Get fitted values
    fits <- predict(mod)[,1]
    # Transform to inversion frequencies
    fits <- fits/2
    
    ### Fit glm comparing observations to fitted values to calculate deviance
    devi <- glm(cbind(mod$data$Inv, 2-mod$data$Inv) ~ fits,
                family = binomial(link = "logit"))
    
    ### Calculate deviance explained by model terms
    dev_mod <- devi$deviance
    dev_null <- devi$null.deviance
    dExp <- ( dev_null - dev_mod ) / dev_null
    
    ### Calculate dispersion
    disp <- dev_null / devi$df.residual
    
    ### Create data frame
    res <- data.frame("Analysis" = Analysis,
                      "Variable" = variable,
                      "Habitat" = habitat,
                      dev_mod, dev_null, 
                      "Dexp" = dExp,
                      "disersion" = disp,
                      looic.mod, looic.null, dlooic, dlooic.se)
    
    } else {
    
    ### Get residual variance
    # Tested values
    res_mod <- VarCorr(mod)$residual$sd[1]^2
    # Model with just sex, maturity and genetic covariance
    res_null <- VarCorr(mod_null)$residual$sd[1]^2
    
    ### Get total variance
    Tvar <- var(mod$data$Pheno)
    
    ### % Variance explained
    Vexp <- ( (res_null - res_mod) / Tvar ) * 100
    
    ### Create data frame
    res <- data.frame("Analysis" = Analysis,
                      "Variable" = variable,
                      "Habitat" = habitat,
                      "Vmod" = res_mod,
                      "Vnull" = res_null,
                      "Vtotal" = Tvar, 
                      Vexp, looic.mod, looic.null, dlooic, dlooic.se)
  }
  return(res)
}


### Pheno~Env
Analysis <- "Pheno_Env"
PE.Vexps <- lapply(Habs, function(h){
  tmp <- lapply(Vars, Var.Exp, habitat = h)
  tmp <- do.call(rbind.data.frame, tmp)
  return(tmp)
})
PE.Vexps <- do.call(rbind.data.frame, PE.Vexps)

ggplot(PE.Vexps)+
  geom_point(aes(x = factor(Habitat, levels = Habs), y = Variable, 
                 size = Vexp, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = FALSE)+
  scale_size_binned(breaks = seq(0, 30, 5), limits = c(0, 30))+
  scale_colour_manual(values = c("grey80", "#008080ff"))+
  labs(title = "Phenotype ~ Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 400 X 120


### Pheno~Inv
Analysis <- "Pheno_Inv"
PI.Vexp <- lapply(substr(Vars, 1, 3), Var.Exp, habitat = "Inv")
PI.Vexps <- do.call(rbind.data.frame, PI.Vexp)

ggplot(PI.Vexps)+
  geom_point(aes(x = Habitat, y = Variable, 
                 size = Vexp, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = FALSE)+
  scale_size_binned(breaks = seq(0, 30, 5), limits = c(0, 30))+
  scale_colour_manual(values = c("grey80", "#008080ff"))+
  labs(title = "Phenotype ~ Inversion")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 120 x 120

### Pheno~Inv+Env
Analysis <- "Pheno_InvEnv"
PIE.Vexps <- lapply(Habs[-3], function(h){
  tmp <- lapply(Vars[-2], Var.Exp, habitat = h)
  tmp <- do.call(rbind.data.frame, tmp)
  return(tmp)
})
PIE.Vexps <- do.call(rbind.data.frame, PIE.Vexps)

ggplot(PIE.Vexps)+
  geom_point(aes(x = factor(Habitat, levels = Habs[-3]), y = Variable, 
                 size = Vexp, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = FALSE)+
  scale_size_binned(breaks = seq(0, 30, 5), limits = c(0, 30))+
  scale_colour_manual(values = c("grey80", "#008080ff"))+
  labs(title = "Phenotype ~ Inversion + Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 300 x 100


### Inv~Env
Analysis <- "Inv_Env"
IE.Vexps <- lapply(Habs, function(h){
  tmp <- lapply(Invs, Var.Exp, habitat = h)
  tmp <- do.call(rbind.data.frame, tmp)
  return(tmp)
})
IE.Vexps <- do.call(rbind.data.frame, IE.Vexps)

ggplot(IE.Vexps)+
  geom_point(aes(x = factor(Habitat, levels = Habs), y = factor(Variable, levels = Invs), 
                 size = Dexp, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = FALSE)+
  scale_size_binned(breaks = seq(0, 1, 0.1), limits = c(0, 1))+
  scale_colour_manual(values = c("grey80", "#008080ff"))+
  labs(title = "Inversion ~ Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 400 X 500 (W x H)



#### Variance decomposition: finding variance attributed to inversions and variance from plasticity (or confounds) ####
### Variance explained by inversions = PIE.Vexps - PE.Vexps

# Merge PIE.Vexps and PE.Vexps
tmp <- left_join(PIE.Vexps, PE.Vexps, by = c("Variable", "Habitat"), suffix = c(".PIE", ".PE")) %>%
  #filter(Habitat != "hab.only") %>% 
  summarise(Variable, Habitat,
            "Vexp.Inv" = Vexp.PIE - Vexp.PE,
            "Vexp.Env" = Vexp.PE,
            "Vexp.cf" = ((Vtotal.PIE - Vnull.PIE) / Vtotal.PIE) * 100) %>%
  mutate("Vunexp" = 100 - Vexp.Inv - Vexp.Env - Vexp.cf) %>%
  pivot_longer(cols = c(Vexp.Inv, Vexp.Env, Vexp.cf, Vunexp), names_to = "type", values_to = "Vexp")

# Manual adjustments of negative Vexp for PC3.pheno rock
tmp[tmp$Variable == "PC3.pheno" & tmp$Habitat == "rock" & tmp$type == "Vexp.Env", "Vexp"] <- 0
tmp[tmp$Variable == "PC3.pheno" & tmp$Habitat == "rock" & tmp$type == "Vunexp", "Vexp"] <- 76.7

### Pie charts
ggplot(tmp, aes(x = "", y = Vexp, fill = type))+
  geom_bar(stat = "identity", position = "fill", col = "grey50")+
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("grey70", "#eff7f7ff", "#008080ff", "white"))+
  facet_grid(rows = vars(factor(Variable, levels = c("PC3.pheno", "PC1.pheno"))), 
             cols = vars(factor(Habitat, levels = c("rock", "boulder", "hab.only"))))+
  theme_void()
# 800 x 400
