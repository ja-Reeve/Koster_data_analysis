################### Figure S8 - Posterior distributions Pheno~Env ####################
### Plot posterior distributions from brms model objects for the Pheno~Env models. Tile
### plots are copies from Fig. 5 in the main text.

### James Reeve - DTU Aqua
### 2025-12-16

### Preparation
rm(list = ls())
try(dev.off())
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(ggpubr)
library(brms)

### Parameters
# Phenotypic variables
Vars <- c("PC1.pheno", "PC2.pheno", "shellLength")
# Habitats
Habs <- c("rock", "boulder", "mud-sand", "hab.only")


#### Posterior distribution plots ####

### Filepath
PATH <- "file/path/"


#### 1. Posterior distribution plots ####
### Rocky shores
pSexRS <- annotate_figure(
  ggarrange(
    plotlist = lapply(Vars, function(pheno){
      load(paste0(PATH, "Pheno_Env/brms_fit.", pheno, "_", Habs[1], "_v2"))
      
      p <- mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
        geom_vline(xintercept = 0, col = "red")+
        labs(x = pheno)+
        xlim(c(-3, 0.5))+
        theme_bw()+
        theme(axis.title.x = element_text(size = 14),
              axis.text.y = element_blank())
      
      return(p)
      }),
    nrow = 1),
  top = "Rocky shores")


### Boulder Fields
pSexBS <- annotate_figure(
  ggarrange(
    plotlist = lapply(Vars, function(pheno){
      load(paste0(PATH, "Pheno_Env/brms_fit.", pheno, "_", Habs[2], "_v2"))
      
      p <- mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
        geom_vline(xintercept = 0, col = "red")+
        labs(x = pheno)+
        xlim(c(-3, 0.5))+
        theme_bw()+
        theme(axis.title.x = element_text(size = 14),
              axis.text.y = element_blank())
      
      return(p)
    }),
    nrow = 1),
  top = "Boulder fields")


### Muddy-sandy Shores
pSexMSS <- annotate_figure(
  ggarrange(
    plotlist = lapply(Vars, function(pheno){
      load(paste0(PATH, "Pheno_Env/brms_fit.", pheno, "_", Habs[3], "_v2"))
      
      p <- mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
        geom_vline(xintercept = 0, col = "red")+
        labs(x = pheno)+
        xlim(c(-3, 0.5))+
        theme_bw()+
        theme(axis.title.x = element_text(size = 14),
              axis.text.y = element_blank())
      
      return(p)
    }),
    nrow = 1),
  top = "Muddy-sand shores")


### Among habitat effects
pSexHab <- annotate_figure(
  ggarrange(
    plotlist = lapply(Vars, function(pheno){
      load(paste0(PATH, "Pheno_Env/brms_fit.", pheno, "_", Habs[4], "_v2"))
      
      p <- mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
        geom_vline(xintercept = 0, col = "red")+
        labs(x = pheno)+
        xlim(c(-3, 0.5))+
        theme_bw()+
        theme(axis.title.x = element_text(size = 14),
              axis.text.y = element_blank())
      
      return(p)
    }),
    nrow = 1),
  top = "Among habitats")

load(paste0(PATH, "Pheno_Env/brms_fit.", Vars[3], "_", Habs[4], "_v2"))

mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
  geom_vline(xintercept = 0, col = "red")+
  labs(x = Vars[3])+
  xlim(c(-3, 0.5))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),
        axis.text.y = element_blank())


#### 2. Tile plot ####

brs <- read.csv("file/path/brms_result_summary_v3.2.csv")
# Subset to Pheno~Env results
PE <- brs[brs$Analysis == "Pheno-Env",]


ggplot(PE[PE$Parameter %in% c("sex", "maturity"),], 
       aes(y = Response, x = Parameter))+
  geom_tile(aes(fill = abs(Estimate)), colour = "grey50")+
  facet_wrap(vars(factor(Habitat, levels = c("rock", "boulder", "mud-sand", "hab.only"))))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 4), breaks = seq(0, 4, 0.8))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 550 x 1000 (plot width x height)
