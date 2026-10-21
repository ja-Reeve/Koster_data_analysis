################### Figure S8 - Posterior distributions Pheno~Inv ####################
### Plot posterior distributions from brms model objects for the Pheno~Inv models. Tile
### plots are copies from Fig. 5 in the main text.

### James Reeve - University of Gothenburg
### 2024-02-07

### Preparation
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(ggpubr)
library(brms)

### Inversion arragements
Invs <- c("LGC1.1.RA", "LGC1.1.RR","LGC1.2.RA", "LGC1.2.RR", "LGC2.1.RA", "LGC2.1.RR", "LGC4.1.RA", "LGC4.1.RR",
          "LGC6.1.2.RA", "LGC6.1.2.RR", "LGC6.1.2.RB", "LGC6.1.2.AB", "LGC6.1.2.BB", "LGC7.1.RA", "LGC7.1.RR",
          "LGC7.2.RA", "LGC7.2.RR", "LGC9.1.RA", "LGC9.1.RR", "LGC10.1.RA", "LGC10.1.RR", "LGC10.2.RA", "LGC10.2.RR",
          "LGC11.1.RA", "LGC11.1.RR", "LGC12.1.RA", "LGC12.1.RR","LGC12.2.RA", "LGC12.2.RR", "LGC12.3.RA", "LGC12.3.RR",
          "LGC12.4.RA", "LGC12.4.RR", "LGC14.1.RA", "LGC14.1.RR","LGC14.3.RA", "LGC14.3.RR", "LGC17.1.RA", "LGC17.1.RR")

#### A: Tile plot ####

### Access data
brs <- read.csv("file/path/brms_result_summary_v3.2.csv")
# Subset to Pheno~Inv+Env results
PI <- brs[brs$Analysis == "Pheno-Inv",]


### Plot
max(abs(PI$Estimate)) # 2.61

ggplot(PI[grep("LGC", PI$Parameter),], 
       aes(x = Response, y = factor(Parameter, levels = rev(Invs))))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative), colour = "grey50")+
  #geom_point(aes(pch = Marginal.case), size = 4)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 2.8), breaks = seq(0, 2.8, 0.4))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 400 x 1200

### Tile plot for sex and maturity
ggplot(PI[PI$Parameter %in% c("sex", "maturity"),], 
       aes(y = Response, x = Parameter))+
  geom_tile(aes(fill = abs(Estimate)), colour = "grey50")+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 3.2), breaks = seq(0, 3.2, 0.8))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 550 x 1000



#### B: Posterior distribution plots ####

### Filepath
PATH <- "file/path/"


### Function to make plot
Post.Plots <- function(variable){
  ### Model fits 
  load(paste0(PATH, "Pheno_Inv/brms_fit.", variable, "_Inv_v2"))
  
  # Plot for each phenotype
   p <- mcmc_plot(mod, type = "areas", pars = "LGC")+
      geom_vline(xintercept = 0, col = "red")+
      labs(x = variable)+
      xlim(c(-2.4, 3.6))+
      theme_bw()+
      theme(axis.title.x = element_text(size = 14))
 
  return(p)
}

Post.Plots.Sex <- function(variable){
  ### Model fits 
  load(paste0(PATH, "Pheno_Inv/brms_fit.", variable, "_Inv_v2"))
  
  # Plot for each phenotype
  p <- mcmc_plot(mod, type = "areas", pars = c("adultjuvenile", "sexM"))+
    geom_vline(xintercept = 0, col = "red")+
    labs(x = variable)+
    xlim(c(-3.5, 0.5))+
    theme_bw()+
    theme(axis.title.x = element_text(size = 14))
  
  return(p)
}


### Plots 
ggarrange(Post.Plots("PC1"), 
          Post.Plots("PC2") + theme(axis.text.y = element_blank()), 
          Post.Plots("shellLength") + theme(axis.text.y = element_blank()), 
          nrow = 1, widths = c(0.4, 0.3, 0.3))
### 740 x 1200

ggarrange(Post.Plots.Sex("PC1"), 
          Post.Plots.Sex("PC2") + theme(axis.text.y = element_blank()), 
          Post.Plots.Sex("shellLength") + theme(axis.text.y = element_blank()), 
          nrow = 1, widths = c(0.4, 0.3, 0.3))
### 1400 x 300
