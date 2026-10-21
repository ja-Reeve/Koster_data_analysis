################### Figure 6 - Pheno~Inv+Env model summary plots ####################
### Plot summaries of the Pheno~Inv+Env models. Three plots are generated per model; 
### i) a tile plot showing the effect size of parameters and informativeness,
### ii) a dot plot indicating the overall reliability of the models,
### iii) pie charts showing varaince explained per sets of parameters.

### James Reeve - DTU Aqua
### 2025-12-15

### Preparation
rm(list = ls())
try(dev.off())
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(brms)

### Inversions
Invs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1.2", #"LGC6.1.2b", "LGC7.1",
          "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2",
          "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.3", "LGC17.1")

### Habitats
Habs <- c("rock", "boulder", "mud-sand", "hab.only")


#### A: access data ####

### 1) Parameter effects summary
brs.sum <- read.csv("file/path/brms_result_summary_v3.2.csv")
## Subset to Pheno~Inv+Env models
PIE.sum <- brs.sum[brs.sum$Analysis == "Pheno-InvEnv",] 
# Subset again to remove muddy-sandy shore models, due to poor model fit
PIE.sum <- PIE.sum[PIE.sum$Habitat != "mud-sand",]

### 2) Model reliability
## Reliability is summarised with two scores: the variance explained (R^2) and the 
## difference in leave-one-out confidence intervals (∆LOOIC).
brs.sup <- read.csv("file/path/brms_LOO_summary_v3.csv")
## Subset to Pheno~Inv+Env models
PIE.sup <- brs.sup[brs.sup$Analysis == "Pheno_InvEnv",]
# Subset again to remove muddy-sandy shore models, due to poor model fit
PIE.sup <- PIE.sup[PIE.sup$Habitat != "mud-sand",]

### Variance decomposition
Vexp <- read.csv("file/path/Variance_explained_summary_Pheno-InvEnv_v4.csv")
## Subset to remove muddy-sandy shore models, due to poor model fit
Vexp <- Vexp[Vexp$Habitat != "mud-sand",]

#### B: Tile plot ####

### 1) Plot inversion effects
# Factorize inversion frequencies
inv.ord <- unique(PIE.sum$Parameter)[grep("LGC", unique(PIE.sum$Parameter))]
# Plot
ggplot(PIE.sum[grepl("LGC", PIE.sum$Parameter),], 
       aes(x = factor(Parameter, levels = inv.ord), y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 1.6), na.value = "#008080ff")+ # Setting NAs created by limits to fill of max
  scale_alpha_manual(values = c(0.3, 1))+
  facet_wrap(vars(factor(Habitat, levels = Habs)), nrow = 4)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"),
        strip.text = element_text(face = "bold", colour = "white", size = 8))
# 1600 x 1000

### 2) Plot environmental effects
## Within each habitat
ggplot(PIE.sum[grepl(".env", PIE.sum$Parameter),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 1.6), na.value = "#008080ff")+ # Setting NAs created by limits to fill of max
  scale_alpha_manual(values = c(0.3, 1))+
  facet_wrap(vars(factor(Habitat, levels = Habs)), nrow = 4)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "#008080ff"),
        strip.text = element_text(face = "bold", colour = "white", size = 8))
# 420 x 650

## Among habitats
ggplot(PIE.sum[PIE.sum$Habitat == "hab.only" & PIE.sum$Parameter %in% c("MSS", "RS"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 1.6), na.value = "#008080ff")+ # Setting NAs created by limits to fill of max
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
# 420 x 300

## Sex and maturity
ggplot(PIE.sum[PIE.sum$Parameter %in% c("sex", "maturity"),], 
       aes(x = Parameter, y = Response))+
  geom_tile(aes(fill = abs(Estimate), alpha = Informative))+
  scale_fill_gradient(low = "#d3eaeaff", high = "#008080ff",
                      limits = c(0, 1.6), na.value = "#008080ff")+ # Setting NAs created by limits to fill of max
  facet_wrap(vars(factor(Habitat, levels = Habs[-3])))+
  scale_alpha_manual(values = c(0.3, 1))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))
# 720 x 300



#### C: Dot plot ####
ggplot(PIE.sup)+
  geom_point(aes(x = factor(Habitat, levels = c("rock", "boulder", "hab.only")), y = Variable,
                 size = Vexp*100, colour = dlooic < -4 | dlooic < -2*dlooic.se), 
             show.legend = T)+
  scale_size_binned(breaks = seq(0, 21, 3), limits = c(0, 21))+
  scale_colour_manual(values = c("#b1d6d5ff", "#008080ff"))+
  labs(title = "Phenotype ~ Inversion + Environment")+
  theme_void()+
  theme(axis.text = element_text(size = 14),
        panel.grid = element_line(linewidth = 0.4, colour = "grey95"))
# 400 x 130


#### D. Variance decomposition pie chart ####

### Reformat as long format
## While doing this, I'm also calculating the unexplained variance
tmp <- Vexp %>% 
  mutate("R2unexp" = 1 - R2full,
         "R2conf" = R2full - (R2env + R2inv + R2cf + R2gbg)) %>%
  select(Variable, Habitat, R2inv, R2env, R2conf, #R2confInv, R2confEnv, 
         R2cf, R2gbg, R2unexp) %>%
  pivot_longer(cols = starts_with("R2"), names_to = "type", values_to = "Vexp")

## Set Vexp < 0 to 0
tmp$Vexp[tmp$Vexp < 0] <- 0


### Pie charts
ggplot(tmp, aes(x = "", y = Vexp, 
                fill = factor(type, 
                              levels = c("R2env", #"R2confEnv", 
                                         "R2inv", "R2conf", "R2gbg", "R2cf", "R2unexp"))))+
  geom_bar(stat = "identity", position = "fill", col = "grey35", linewidth = 0.1)+
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("#00FF15", #"#d3eaeaff", 
                               "#3182bd", "#d3eaeaff", #"#8dbfe2",
                               "#fdae6b", "grey80", "white"))+
  labs(fill = "Component")+
  facet_grid(rows = vars(factor(Variable, levels = c("shellLength", "PC2.pheno", "PC1.pheno"))),
             cols = vars(factor(Habitat, levels = Habs[-3])))+
  theme_void()
# 1200 x 900
