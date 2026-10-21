#################### Figure 5: Bayesian regression summary plots ##########################
### Summary plots for Baysian regression models (Pheno~Env, Pheno~Inv, & Inv~Env). Plots
### inlcude: A, C, & E tile plots of effect sizes, and B, D, & F dot plots showing variance
### explained per model (or deviance explained in the logistic fits for Inv~Env).

### James Reeve - University of Gothenburg
### 2025-12-15

### Preparation
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(forcats)
library(brms)

### Inversions
Invs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", "LGC6.1.2b", "LGC7.1",
          "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2",
          "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3", "LGC17.1")

### Habitats
Habs <- c("rock", "boulder", "mud-sand", "hab.only")


#### A: access data ####

### 1) Parameter effects summary
brs.sum <- read.csv("/Users/james/Documents/snail-scape/brms_result_summary_v3.2.csv")

## Subset to each type of model
PE.sum <- brs.sum[brs.sum$Analysis == "Pheno-Env",] # Pheno~Env models
PI.sum <- brs.sum[brs.sum$Analysis == "Pheno-Inv",] # Pheno~Inv models
IE.sum <- brs.sum[brs.sum$Analysis == "Inv-Env",] # Inv~Env models


### 2) Model reliability
## Reliability is summarised with two scores: the variance explained (R^2) and the 
## difference in leave-one-out confidence intervals (∆LOOIC).
brs.sup <- read.csv("/file/path/brms_LOO_summary_v3.csv")

## Subset to each type of model
PE.sup <- brs.sup[brs.sup$Analysis == "Pheno_Env",]
PI.sup <- brs.sup[brs.sup$Analysis == "Pheno_Inv",]
IE.sup <- brs.sup[brs.sup$Analysis == "Inv_Env",]


#### B: Pheno~Env plots ####

### 1) Tile plots
## Per habitat effects
ggplot(PE.sum[grepl(".env", PE.sum$Parameter),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff")+
  scale_alpha_manual(values = c(0.3, 1))+
  facet_wrap(vars(factor(Habitat, levels = Habs[1:3])), nrow = 3)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"),
        strip.text = element_text(face = "bold", colour = "white", size = 8))
# 420 x 600 (W x H)

## Among habitat effects
ggplot(PE.sum[PE.sum$Habitat == "hab.only" & PE.sum$Parameter %in% c("MSS", "RS"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5))+
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
# 520 x 400

## Co-factors
# Note: while co-factors were included in all models, I am only presenting their effects
# on the among-habitat models, as they include all samples.
ggplot(PE.sum[PE.sum$Habitat == "hab.only" & PE.sum$Parameter %in% c("sex", "maturity"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5),
                      na.value = "#008080ff")+
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
# 500 x 400


### 2) Dot plot
ggplot(PE.sup)+
  geom_point(aes(x = factor(Habitat, levels = Habs), y = Variable,
                 size = Vexp, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = FALSE)+
  scale_size_binned(breaks = seq(0, 0.08, 0.02), limits = c(0, 0.08))+
  scale_colour_manual(values = c("#b1d6d5ff", "#008080ff"))+
  labs(title = "Phenotype ~ Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 500 x 150  
  


#### C: Pheno~Inv plots ####

### 1) Tile plots
## Vector to store order of inversions in Parameter names
ord <- match(substr(PI.sum$Parameter, 1, (nchar(PI.sum$Parameter)-2)), Invs)
# Reverse the order
ord <- 20 - ord

## Plot
ggplot(PI.sum, 
       aes(x = Response, y = fct_reorder(Parameter, ord) ))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 2.9), breaks = seq(0, 2.9, 0.4))+
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
# 350 x 1100 (W x H)


### 2) Dot plot
ggplot(PI.sup)+
  geom_point(aes(x = Habitat, y = Variable,
                 size = Vexp), #colour = dlooic < -4 | dlooic < -2*dlooic.se), Note: all models were informative!
             show.legend = FALSE, colour = "#008080ff")+
  scale_size_binned(breaks = seq(0, 0.08, 0.02), limits = c(0, 0.08))+
  #scale_colour_manual(values = c("#b1d6d5ff", "#008080ff"))+
  labs(title = "Phenotype ~ Inversion")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 200 x 150  



#### D: Inv~Env plots ####

### 1) Tile plots
## Per habitat effects
ggplot(IE.sum[grepl(".env", IE.sum$Parameter),], 
       aes(x = Parameter, y = factor(Response, levels = rev(Invs)) ))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 1), breaks = seq(0, 1, 0.2),
                      na.value = "#008080ff")+  # Setting NAs created by limits to fill of max
  scale_alpha_manual(values = c(0.3, 1))+
  facet_wrap(vars(factor(Habitat, levels = Habs[1:3])), ncol = 3)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"),
        strip.text = element_text(face = "bold", colour = "white", size = 8))
# 800 x 1000 (W x H)

## Among habitat effects
ggplot(IE.sum[IE.sum$Parameter %in% c("MSS", "RS"),], 
       aes(x = Parameter, y = factor(Response, levels = rev(Invs)) ))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 3), breaks = seq(0, 3.6, 0.6),
                      na.value = "#008080ff")+  # Setting NAs created by limits to fill of max
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"),
        strip.text = element_text(face = "bold", colour = "white", size = 8))
# 400 x 1000 (W x H)

## Co-factors
ggplot(IE.sum[IE.sum$Parameter %in% c("sex", "maturity"),], 
       aes(x = Parameter, y = factor(Response, levels = rev(Invs)) ))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 3), breaks = seq(0, 3, 0.6),
                      na.value = "#008080ff")+  # Setting NAs created by limits to fill of max
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"),
        strip.text = element_text(face = "bold", colour = "white", size = 8))
# 400 x 1000 (W x H)


### 2) Dot plot
ggplot(IE.sup)+
  geom_point(aes(x = factor(Habitat, levels = Habs), y = factor(Variable, levels = rev(Invs)),
                 size = Vexp*100, colour = dlooic < -4 | dlooic < -2*dlooic.se), show.legend = T)+
  scale_size_binned(breaks = seq(0, 80, 10), limits = c(0, 80))+
  scale_colour_manual(values = c("#b1d6d5ff", "#008080ff"))+
  labs(title = "Inversion ~ Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 500 x 620  
