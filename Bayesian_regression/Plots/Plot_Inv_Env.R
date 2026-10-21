#################################### Plot summary of Inv~Env models ################################
### Plot summaries for the Inv~Env models. This includes a tile plot of all environmental PCs and the 
### corresponding posterior distributions. This was used to create Fig. 6E and Fig. S5.

### James Reeve - University of Gothenburg
### 2024-02-09

### Preparation
rm(list = ls())
dev.off()
setwd("~")
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(ggpubr)
library(brms)

### Inversions
Invs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", "LGC6.1.2b", "LGC7.1",
          "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2",
          "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3", "LGC17.1")


#### A: Tile plot ####

### Access data
brs <- read.csv("/path/to/brms/summary/brms_result_summary_v2.csv")
# Subset to Pheno~Inv+Env results
IE <- brs[brs$Analysis == "Inv-Env",]
# Subset into envrionmental PCs and habitat effects data
IE.envs <- IE[IE$Habitat != "hab",]
IE.hab <- IE[IE$Habitat == "hab",]


### Plot environmental PCs effects
max(abs(IE.envs$Estimate)) # 38.85 - but this distribution has very long tails due to wide distributions of LGC6.1/2b in MSS and BS
# Set maximum effect size to 1
IE.envs[IE.envs$Estimate >= 1, "Estimate"] <- 1
IE.envs[IE.envs$Estimate <= -1, "Estimate"] <- -1

ggplot(IE.envs, 
       aes(x = Parameter, y = factor(Response, levels = Invs)))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  facet_wrap(vars(factor(Habitat, levels = c("RS", "BS", "MSS"))), nrow = 1)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 800 x 1000


### Plot habitat effects
max(abs(IE.hab$Estimate)) # 20.07 - also bias from LGC6.1/2b
# Set maximum effect size to 1
IE.hab[IE.hab$Estimate >= 3, "Estimate"] <- 3
IE.hab[IE.hab$Estimate <= -3, "Estimate"] <- -3

# Habitat effects
p1 <- ggplot(IE.hab[IE.hab$Parameter %in% c("RS", "MSS"),], 
             aes(x = Parameter, y = factor(Response, levels = Invs)))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 3), breaks = seq(0, 3, 0.6))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

# Sex & maturity
p2 <- ggplot(IE.hab[IE.hab$Parameter %in% c("sex", "maturity"),], 
             aes(x = Parameter, y = factor(Response, levels = Invs)))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 3), breaks = seq(0, 3, 0.6))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_blank())

ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "right", 
          widths = c(1, 0.65))
# 550 x 1000



#### B: Posterior distribution plots ####

### Filepath
PATH <- "/path/to/brms/fits/"


### Function to make plot
Post.Plots <- function(variable, habitat){
  ### Model fits 
  load(paste0(PATH, "Inv_Env/brms_fit.", variable, "_", habitat))
  
  # Plot for Env PCs
  if(habitat %in% c("rock", "boulder", "mud-sand")){
    p <- mcmc_plot(mod, type = "areas", pars = "PC")+
      geom_vline(xintercept = 0, col = "red")+
      xlim(c(-1.5, 1.5))+ # -1.5, 1.5
      theme_bw()+
      theme(axis.title = element_text(size = 8),
            axis.text = element_blank(),
            plot.margin = unit(c(0,0,0,0), 'cm'))
  }
  
  # Plot for sex and maturity
  if(habitat == "hab.only"){
    p.tmp1 <- mcmc_plot(mod, type = "areas", variable = c("b_HabitatmudMsand", "b_Habitatrock"))+
      geom_vline(xintercept = 0, col = "red")+
      xlim(c(-3, 3))+ #c(-3, 3)
      theme_bw()+
      theme(axis.title = element_text(size = 8),
            axis.text = element_blank(),
            plot.margin = unit(c(0.1,0,0.1,0), 'cm'))
    
    p.tmp2 <- mcmc_plot(mod, type = "areas", variable = c("b_adultjuvenile", "b_sexM"))+
      geom_vline(xintercept = 0, col = "red")+
      xlim(c(-3, 3))+
      theme_bw()+
      theme(axis.title = element_text(size = 8),
            axis.text = element_blank(),
            plot.margin = unit(c(0.1,0,0.1,0), 'cm'))
    
    p <- ggarrange(p.tmp1, p.tmp2, nrow = 1)
  }
  
  return(p)
}


### Plot inversions other than LGC6.1/2
### Doing this so I can set wider axes of LGC6.1/2
# Environmental PCs
ggarrange(ggarrange(plotlist = lapply(Invs[-c(5,6)], Post.Plots, habitat = "rock"), ncol = 1),
          ggarrange(plotlist = lapply(Invs[-c(5,6)], Post.Plots, habitat = "boulder"), ncol = 1),
          ggarrange(plotlist = lapply(Invs[-c(5,6)], Post.Plots, habitat = "mud-sand"), ncol = 1),
          nrow = 1, ncol = 3)

# Habitat, maturity and sex
ggarrange(plotlist = lapply(Invs[-c(5,6)], Post.Plots, habitat = "hab.only"), ncol = 1)


### Plot LGC6.1/2
### This is too hard to fit onto a common set of axes, so I'm doing this part manually

# LGC6.1/2
p.A.RS <- Post.Plots("LGC6.1.2", "rock")+xlim(c(-1, 1))+theme(axis.text.x = element_text())
p.A.BS <- Post.Plots("LGC6.1.2", "boulder")+xlim(c(-5, 5))+theme(axis.text.x = element_text())
p.A.MSS <- Post.Plots("LGC6.1.2", "mud-sand")+xlim(c(-80, 80))+theme(axis.text.x = element_text())

ggarrange(p.A.RS, p.A.BS, p.A.MSS, widths = c(1, 2, 4), nrow = 1)
# 1100 x 120
Post.Plots("LGC6.1.2", "hab.only")
# 600 x 100


# LGC6.1/2b
p.B.RS <- Post.Plots("LGC6.1.2b", "rock")+xlim(c(-1, 1))+theme(axis.text.x = element_text())
p.B.BS <- Post.Plots("LGC6.1.2b", "boulder")+xlim(c(-80, 80))+theme(axis.text.x = element_text())
p.B.MSS <- Post.Plots("LGC6.1.2b", "mud-sand")+xlim(c(-120, 120))+theme(axis.text.x = element_text())
ggarrange(p.B.RS, p.B.BS, p.B.MSS, widths = c(1, 6, 10), nrow = 1)
# 1100 x 120

load(paste0(PATH, "Inv_Env/brms_fit.LGC6.1.2b_hab.only"))
p.tmp1 <- mcmc_plot(mod, type = "areas", variable = c("b_HabitatmudMsand", "b_Habitatrock"))+
  geom_vline(xintercept = 0, col = "red")+
  xlim(c(-5, 50))+ #c(-3, 3)
  theme_bw()+
  theme(axis.title = element_text(size = 8),
        axis.text = element_blank(),
        plot.margin = unit(c(0.1,0,0.1,0), 'cm'))

p.tmp2 <- mcmc_plot(mod, type = "areas", variable = c("b_adultjuvenile", "b_sexM"))+
  geom_vline(xintercept = 0, col = "red")+
  xlim(c(-3, 3))+
  theme_bw()+
  theme(axis.title = element_text(size = 8),
        axis.text = element_blank(),
        plot.margin = unit(c(0.1,0,0.1,0), 'cm'))

p <- ggarrange(p.tmp1, p.tmp2, nrow = 1)
# 600 x 100
