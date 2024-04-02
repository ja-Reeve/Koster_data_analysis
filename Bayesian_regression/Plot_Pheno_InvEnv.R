####################################### Plot Pheno~InvEnv results ############################################
### Plot summaries for the Pheno~Inv+Env models. This includes a tile plot of all inversions and environmental 
### PCs, along with a dot plot showing variance explained. This was used to plot Fig. 7 and Fig. S6.

### James Reeve - University of Gothenburg
### 2024-02-02

### Preparation
rm(list = ls())
dev.off()
setwd("~")
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(brms)

### Inversions
Invs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", "LGC6.1.2b", "LGC7.1",
          "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2",
          "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3", "LGC17.1")


#### A: Tile plot ####

### Access data
brs <- read.csv("/path/to/brms/summary/brms_result_summary.csv")
# Subset to Pheno~Inv+Env results
PIE <- brs[brs$Analysis == "Pheno-InvEnv",]
# Subset further by habitat
PIE.RS <- PIE[PIE$Habitat == "RS",]
PIE.BS <- PIE[PIE$Habitat == "BS",]
PIE.hab <- PIE[PIE$Habitat == "hab",]


### Plot inversion effects
ggplot(PIE[grepl("LGC", PIE$Parameter),], 
       aes(x = factor(Parameter, levels = Invs), y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0,2.25), breaks = seq(0,2.25, 0.5))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  facet_wrap(vars(factor(Habitat, levels = c("RS", "BS", "hab"))), nrow = 3)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 2000 X 800

### Plot effects of envrionmental PCs
ggplot(PIE[grepl("PC", PIE$Parameter),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0,2.25), breaks = seq(0,2.25, 0.5))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  facet_wrap(vars(factor(Habitat, levels = c("RS", "BS"))), nrow = 2)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 600 x 520

### Plot among habitat effect
ggplot(PIE[PIE$Parameter %in% c("MSS", "RS"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0,2.25), breaks = seq(0,2.25, 0.5))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  facet_wrap(vars(factor(Habitat, levels = c("RS", "BS", "hab"))), nrow = 3)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 450 x 270



#### B: Variance explained dot-plot ####

### Filepath
PATH <- "/path/to/brms/fits/"


### Function to calcaute varaince explained by each model fit
Var.Exp <- function(variable, habitat){
  ### Model fits (Pheno~Inv+Env models)
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_", habitat)) # Full model
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_", habitat, "_null")) # Null model
  
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
  res <- data.frame("Analysis" = "Pheno-InvEnv",
                    "Variable" = variable,
                    "Habitat" = habitat,
                    "Vmod" = res_mod,
                    "Vnull" = res_null,
                    "Vtotal" = Tvar, 
                    Vexp, looic.mod, looic.null, dlooic, dlooic.se)
  
  return(res)
}


### Get variance explained for each model
PIE.Vexps <- lapply(c("rock", "boulder", "hab.only"), function(h){
  tmp <- lapply(c("PC1.pheno", "PC3.pheno"), Var.Exp, habitat = h)
  tmp <- do.call(rbind.data.frame, tmp)
  return(tmp)
})
PIE.Vexps <- do.call(rbind.data.frame, PIE.Vexps)

ggplot(PIE.Vexps)+
  geom_point(aes(x = factor(Habitat, levels = c("rock", "boulder", "hab.only")), y = Variable, 
                 size = Vexp, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = FALSE)+
  scale_size_binned(breaks = seq(0, 30, 5), limits = c(0, 30))+
  scale_colour_manual(values = c("#eff7f7ff", "#008080ff"))+
  labs(title = "Phenotype ~ Inversion + Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 300 x 100
