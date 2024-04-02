##################################### Plot Pheno~Env results ##################################
### Plot summaries for the Pheno~Env models. This includes a tile plot of all environmental PCs 
### and the corresponding posterior distributions. This was used in Fig. 6B and Fig. S3.

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

#### A: Tile plot ####

### Access data
brs <- read.csv("/path/to/brms/summary/brms_result_summary_v2.csv")
# Subset to Pheno~Inv+Env results
PE <- brs[brs$Analysis == "Pheno-Env",]
# Subset into envrionmental PCs and habitat effects data
PE.envs <- PE[PE$Habitat != "hab",]
PE.hab <- PE[PE$Habitat == "hab",]


### Plot environmental PCs effects
max(abs(PE.envs$Estimate)) # 0.5

ggplot(PE.envs, 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0,0.5), breaks = seq(0,0.5, 0.1))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  facet_wrap(vars(factor(Habitat, levels = c("RS", "BS", "MSS"))), nrow = 3)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"))
# 450 x 700

### Plot habitat effects
max(abs(PE.hab$Estimate)) # 1.88

# Habitat effects
p1 <- ggplot(PE.hab[PE.hab$Parameter %in% c("RS", "MSS"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 2), breaks = seq(0, 2, 0.5))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

# Sex & maturity
p2 <- ggplot(PE.hab[PE.hab$Parameter %in% c("sex", "maturity"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  #geom_point(aes(pch = Marginal.case), size = 6)+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 2), breaks = seq(0, 2, 0.5))+
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(NA, 8))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_blank())

ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "none", widths = c(1, 0.75))
# 450 x 270



#### B: Posterior distribution plots ####

### Filepath
PATH <- "/path/to/brms/fits/"


### Function to make plot
Post.Plots <- function(variable, habitat){
  ### Model fits 
  load(paste0(PATH, "Pheno_Env/brms_fit.", variable, "_", habitat))
  
  # Plot for Env PCs
  if(habitat %in% c("rock", "boulder", "mud-sand")){
    p <- mcmc_plot(mod, type = "areas", pars = "PC")+
      geom_vline(xintercept = 0, col = "red")+
      labs(x = variable)+
      xlim(c(-0.7, 0.4))+
      theme_bw()+
      theme(axis.title.x = element_text(size = 14),
            axis.text.y = element_blank())
  }
  
  # Plot for sex and maturity
  if(habitat == "hab.only"){
    p.tmp1 <- mcmc_plot(mod, type = "areas", variable = c("b_HabitatmudMsand", "b_Habitatrock"))+
      geom_vline(xintercept = 0, col = "red")+
      labs(title = "Habitat effects", x = variable)+
      xlim(c(-2.3, 0.7))+
      theme_bw()+
      theme(axis.title.x = element_text(size = 14),
            axis.text.y = element_blank())
    
    p.tmp2 <- mcmc_plot(mod, type = "areas", variable = c("b_adultjuvenile", "b_sexM"))+
      geom_vline(xintercept = 0, col = "red")+
      labs(title = "Maturity and sex", x = variable)+
      xlim(c(-2.3, 0.7))+
      theme_bw()+
      theme(axis.title.x = element_text(size = 14),
            axis.text.y = element_blank())
    
    p <- ggarrange(p.tmp1, p.tmp2, nrow = 2)
  }
  
  return(p)
}


### Plots of environmental PCs
ggarrange(Post.Plots("PC1.pheno", "rock"), Post.Plots("PC2.pheno", "rock"), Post.Plots("PC3.pheno", "rock"),
          Post.Plots("PC1.pheno", "boulder"), Post.Plots("PC2.pheno", "boulder"), Post.Plots("PC3.pheno", "boulder"),
          Post.Plots("PC1.pheno", "mud-sand"), Post.Plots("PC2.pheno", "mud-sand"), Post.Plots("PC3.pheno", "mud-sand"),
          nrow = 3, ncol = 3)
### 800 x 600


### Plots of habitat and sex effects
ggarrange(Post.Plots("PC1.pheno", "hab.only"), Post.Plots("PC2.pheno", "hab.only"), Post.Plots("PC3.pheno", "hab.only"),
          nrow = 1)
### 800 x 400
