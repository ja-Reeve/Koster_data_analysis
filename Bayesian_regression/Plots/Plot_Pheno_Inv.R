##################################### Plot summary of Pheno~Inv models #########################################
### Plot summaries for the Pheno~Inv models. This includes a tile plot of all inversions and the corresponding 
### posterior distributions. This is used for Fig. 6C and Fig. S4.

### James Reeve - University of Gothenburg
### 2024-02-07

### Preparation
rm(list = ls())
dev.off()
setwd("~")
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
brs <- read.csv("/path/to/brms/summary/brms_result_summary_v2.csv")
# Subset to Pheno~Inv+Env results
PI <- brs[brs$Analysis == "Pheno-Inv",]


### Plot
max(abs(PI$Estimate)) # 2.61

ggplot(PI[grep("LGC", PI$Parameter),], 
       aes(x = Response, y = factor(Parameter, levels = rev(Invs))))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative), colour = "grey50")+
  #geom_point(aes(pch = Marginal.case), size = 4)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff", limits = c(0, 2.8), breaks = seq(0, 2.8, 0.4))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 400 x 1200



#### B: Posterior distribution plots ####

### Filepath
PATH <- "/path/to/brms/fits/"


### Function to make plot
Post.Plots <- function(variable){
  ### Model fits 
  load(paste0(PATH, "Pheno_Inv/brms_fit.", variable, "_Inv"))
  
  # Plot for each phenotype
   p <- mcmc_plot(mod, type = "areas", pars = "LGC")+
      geom_vline(xintercept = 0, col = "red")+
      labs(x = variable)+
      xlim(c(-2.4, 3.2))+
      theme_bw()+
      theme(axis.title.x = element_text(size = 14))
 
  return(p)
}


### Plots 
ggarrange(Post.Plots("PC1"), Post.Plots("PC2"), Post.Plots("PC3"), nrow = 1)
### 900 x 1200
