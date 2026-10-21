################### Figure - Posterior distributions Pheno~Inv + Env ####################
### Plot posterior distributions from brms model objects for the Pheno~Inv models. Tile
### plots are copies from Fig. 6 in the main text. Among habitat models were used in
### Fig. S8C. Other posterior distributions are in Fig. S12.

### James Reeve - DTU Aqua
### 2026-02-06

### Preparation
rm(list = ls())
try(dev.off())
setwd("~")
options(stringsAsFactors = FALSE)
mem.maxVSize(vsize = 20000) # Increase vector memory to handle plots being generated

### Packages
library(ggplot2)
library(ggpubr)
library(brms)

### Parameters
# Phenotypic variables
Vars <- c("PC1.pheno", "PC2.pheno", "shellLength")
# Habitats
Habs <- c("rock", "boulder", "mud-sand", "hab.only")
# Inversion arragements
Invs <- c("LGC1.1.RA", "LGC1.1.RR","LGC1.2.RA", "LGC1.2.RR", "LGC2.1.RA", "LGC2.1.RR", "LGC4.1.RA", "LGC4.1.RR",
          "LGC6.1.2.RA", "LGC6.1.2.RR", "LGC6.1.2.RB", "LGC6.1.2.AB", "LGC6.1.2.BB", "LGC7.1.RA", "LGC7.1.RR",
          "LGC7.2.RA", "LGC7.2.RR", "LGC9.1.RA", "LGC9.1.RR", "LGC10.1.RA", "LGC10.1.RR", "LGC10.2.RA", "LGC10.2.RR",
          "LGC11.1.RA", "LGC11.1.RR", "LGC12.1.RA", "LGC12.1.RR","LGC12.2.RA", "LGC12.2.RR", "LGC12.3.RA", "LGC12.3.RR",
          "LGC12.4.RA", "LGC12.4.RR", "LGC14.1.RA", "LGC14.1.RR","LGC14.3.RA", "LGC14.3.RR", "LGC17.1.RA", "LGC17.1.RR")

### Filepath
PATH <- "file/path/"

#### 1. Plot functions ####
### A) Inversion effects
Post.Plots.Inv <- function(variable, habitat, x.limits = c(-0.5, 0.5)){
  ### Model fits 
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_", habitat, "_v3"))
  
  # Plot for each phenotype
  p <- mcmc_plot(mod, type = "areas", pars = "LGC")+
    geom_vline(xintercept = 0, col = "red")+
    labs(x = "")+
    xlim(x.limits)+
    theme_bw()#+
    #theme(axis.title.x = element_text(size = 14))
  
  return(p)
}

### B) Environmental effects
Post.Plots.Env <- function(variable, habitat, x.limits = c(-0.5, 0.5)){
  ### Model fits 
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_", habitat, "_v3"))
  
  # Plot for each phenotype
  p <- mcmc_plot(mod, type = "areas", pars = ".env")+
    geom_vline(xintercept = 0, col = "red")+
    labs(x = "")+
    xlim(x.limits)+
    theme_bw()#+
  #theme(axis.title.x = element_text(size = 14))
  
  return(p)
}

### C) Among habitat effects
Post.Plots.Hab <- function(variable, x.limits = c(-0.5, 0.5)){
  ### Model fits 
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_hab.only_v3"))
  
  # Plot for each phenotype
  p <- mcmc_plot(mod, type = "areas", pars = "Habitat")+
    geom_vline(xintercept = 0, col = "red")+
    labs(x = "")+
    xlim(x.limits)+
    theme_bw()#+
  #theme(axis.title.x = element_text(size = 14))
  
  return(p)
}

### D) Sex and maturity effects
Post.Plots.Sex <- function(variable, habitat, x.limits = c(-0.5, 0.5)){
  ### Model fits 
  load(paste0(PATH, "Pheno_InvEnv/brms_fit.", variable, "_InvEnv_", habitat, "_v3"))
  
  # Plot for each phenotype
  p <- mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
    geom_vline(xintercept = 0, col = "red")+
    labs(x = variable)+
    xlim(x.limits)+
    theme_bw()+
    theme(axis.title.x = element_text(size = 14))
  
  return(p)
}


#### 2. Posterior distribution plots ####
### Note: this will be run in a multipanel plot, but for fine tuning, independent elements
### can be manipulated in a graphics editor. I recommend the following clipboard dimensions.
### Inversion = 1200 x 750
### Environment = 1200 x 250
### Habitat, sex & maturity = 1200 x 250 

### A) Rocky shores
## Inversions
p1A <- ggarrange(Post.Plots.Inv(Vars[1], Habs[1], x.limits = c(-2.5, 2.5))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Inv(Vars[2], Habs[1], x.limits = c(-1.5, 2)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Inv(Vars[3], Habs[1], x.limits = c(-3.5, 3.5)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Environment
p1B <- ggarrange(Post.Plots.Env(Vars[1], Habs[1], x.limits = c(-0.4, 0.3))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Env(Vars[2], Habs[1], x.limits = c(-0.4, 0.3)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Env(Vars[3], Habs[1], x.limits = c(-0.4, 0.3)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Sex and maturity
p1C <- ggarrange(Post.Plots.Sex(Vars[1], Habs[1], x.limits = c(-1.5, 0.5))+
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Sex(Vars[2], Habs[1], x.limits = c(-0.1, 0.4)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Sex(Vars[3], Habs[1], x.limits = c(-3.5, 0.5)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Multipanel plot (500 x 2700 W x H)
Fig.S7.RS <- ggarrange(p1A, p1B, p1C,
          ncol = 1, heights = c(1, 0.1, 0.1))
ggsave("Fig.S7.RS.png", Fig.S7.RS,
       device = "png", path = "file/path/Figures/",
       width = 200, height = 300, units = "mm")
rm(p1A, p1B, p1C, Fig.S7.RS)


### B) Boulder fields
## Inversions
p2A <- ggarrange(Post.Plots.Inv(Vars[1], Habs[2], x.limits = c(-7, 5))+
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Inv(Vars[2], Habs[2], x.limits = c(-4, 4)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Inv(Vars[3], Habs[2], x.limits = c(-5, 4)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Environment
p2B <- ggarrange(Post.Plots.Env(Vars[1], Habs[2], x.limits = c(-0.4, 0.3))+
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Env(Vars[2], Habs[2], x.limits = c(-0.4, 0.3)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Env(Vars[3], Habs[2], x.limits = c(-0.4, 0.3)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Sex and maturity
p2C <- ggarrange(Post.Plots.Sex(Vars[1], Habs[2], x.limits = c(-2, 0.5))+
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Sex(Vars[2], Habs[2], x.limits = c(-0.5, 1)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 Post.Plots.Sex(Vars[3], Habs[2], x.limits = c(-4.5, 0.5)) + 
                   theme(axis.text.y = element_blank(),
                         plot.margin = margin(0,0,0,0)), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Multipanel plot (500 x 2700 W x H)
Fig.S7.BS <- ggarrange(p2A, p2B, p2C,
          ncol = 1, heights = c(1, 0.1, 0.1))
ggsave("Fig.S7.BS.png", Fig.S7.BS,
       device = "png", path = "file/path/Figures/",
       width = 200, height = 300, units = "mm")
rm(p2A, p2B, p2C, Fig.S7.BS)


### C) Muddy-sandy shores
## Inversions
p3A <- ggarrange(Post.Plots.Inv(Vars[1], Habs[3], x.limits = c(-7, 6))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Inv(Vars[2], Habs[3], x.limits = c(-5, 7)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Inv(Vars[3], Habs[3], x.limits = c(-12, 10)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Environment
p3B <- ggarrange(Post.Plots.Env(Vars[1], Habs[3], x.limits = c(-0.4, 0.3))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Env(Vars[2], Habs[3], x.limits = c(-0.4, 0.3)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Env(Vars[3], Habs[3], x.limits = c(-0.4, 0.3)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Sex and maturity
p3C <- ggarrange(Post.Plots.Sex(Vars[1], Habs[3], x.limits = c(-1.5, 0.5))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Sex(Vars[2], Habs[3], x.limits = c(-0.5, 1)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Sex(Vars[3], Habs[3], x.limits = c(-3, 0.5)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Multipanel plot (500 x 2700 W x H)
Fig.S7.MSS <- ggarrange(p3A, p3B, p3C,
          ncol = 1, heights = c(1, 0.1, 0.1))
ggsave("Fig.S7.MSS.png", Fig.S7.MSS,
       device = "png", path = "file/path/Figures/",
       width = 200, height = 300, units = "mm")
rm(p3A, p3B, p3C, Fig.S7.MSS)


### D) Among habitats
## Inversions
p4A <- ggarrange(Post.Plots.Inv(Vars[1], Habs[4], x.limits = c(-2.5, 3))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Inv(Vars[2], Habs[4], x.limits = c(-1, 1)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Inv(Vars[3], Habs[4], x.limits = c(-3, 4)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Environment
p4B <- ggarrange(Post.Plots.Hab(Vars[1], x.limits = c(-1.2, 0.3))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Hab(Vars[2], x.limits = c(-1.2, 0.3)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Hab(Vars[3], x.limits = c(-1.2, 0.3)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Sex and maturity
p4C <- ggarrange(Post.Plots.Sex(Vars[1], Habs[4], x.limits = c(-1.5, 0.5))+
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Sex(Vars[2], Habs[4], x.limits = c(-0.5, 0.5)) + 
                   theme(axis.text.y = element_blank()), 
                 Post.Plots.Sex(Vars[3], Habs[4], x.limits = c(-3.5, 0.5)) + 
                   theme(axis.text.y = element_blank()), 
                 nrow = 1, widths = c(0.3, 0.3, 0.3))

# Multipanel plot (500 x 2700 W x H)
Fig.S7.hab <- ggarrange(p4A, p4B, p4C,
          ncol = 1, heights = c(1, 0.1, 0.1))
ggsave("Fig.S7.hab.png", Fig.S7.hab,
       device = "png", path = "file/path/Figures/",
       width = 200, height = 300, units = "mm")
rm(p4A, p4B, p4C, Fig.S7.hab)
