######################## Figure 4: phenotypic variation and ecotypes #############################
### Ploting scatterplots of phenotypic PCs and shell length vs. PC1, using points that are coloured
### by ecotypes. Sample ecotypes are either inffered through habitat characteristics or assinged
### from previous studies at the same location (Westram et al. 2018; Koch et al., 2022, 2023).
### Ecotype categories and colours:
###                    - Known Wave    = dark red ("#d73027")
###                    - Inferred Wave = light red ("#fc8d59")
###                    - Known Crab    = saturated blue ("#4585d9")
###                    - Inferred Crab = desaturated blue ("#91bfdb")
###                    - Intermediate phenotypes = "grey50"

### James Reeve - DTU Aqua
### 2025-09-18

### Preparation
rm(list = ls())
try(dev.off())
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

##### 1. Access data ####
PATH1 <- "file/path/PCAmix/"
PATH2 <- "file/path/Phenotypes/"

### A: Access data
# Phenotypic PCA
PCs <- read.csv(paste0(PATH1, "principal_components/Pheno_PC1-8.v3.csv"))
Pvar <- read.csv(paste0(PATH1, "eigenvalues/Pheno_variance_explained.csv"))

# Infered ecotype labels
CW <- read.csv(paste0(PATH2, "Inferring_ecotype_from_Amin.csv"), sep = ";")

# Shell lengths
SL <- read.csv(paste0(PATH2, "KT_ShellShaper_parameters.v2.csv"))[,c("snailID", "shellLength")]

# Known Crab and Wave sites
Wave <- "KT007"
Crab <- "KT009"


### B: Wrangle the data
# Add 'Site' to ecotype labels
CW$Site <- substr(CW$Snail_ID, 1, 5)

# List of Wave sites
wave.sites <- CW[CW$ecotype == "Wave", "Site"]
# List of Crab sites
crab.sites <- CW[CW$ecotype == "Crab", "Site"]

# Add 'Site' to PCs
colnames(PCs)[1] <- "Snail_ID"
PCs$Site <- substr(PCs$Snail_ID, 1, 5)

# Remove "_ac" from the snail IDs for shellLength
SL$snailID <- substr(SL$snailID, 1, 6)


### C: Merge data

# Add Ecotype
PCs$Ecotype <- "Other"

PCs[PCs$Site %in% crab.sites, "Ecotype"] <- "Crab"
PCs[PCs$Site %in% wave.sites, "Ecotype"] <- "Wave"

# Add shellLength
PCs$shellLength <- SL[SL$snailID %in% PCs$Snail_ID, "shellLength"]



#### 2. Plot known and inferrd sites ####
### A: Shell length vs. PC1
p1a <- ggplot(PCs, aes(x = shellLength, y = PC1))+
  geom_point(colour = "grey50", alpha = 0.4, size = 0.8)+
  geom_point(data = PCs[PCs$Ecotype == "Wave", ], colour = "#fc8d59", alpha = 0.3)+
  geom_point(data = PCs[PCs$Ecotype == "Crab", ], colour = "#91bfdb", alpha = 0.3)+
  geom_point(data = PCs[PCs$Site == Wave, ], colour = "#d73027", pch = 15)+
  geom_point(data = PCs[PCs$Site == Crab, ], colour = "#4585d9", pch = 15)+
  labs(x = "Shell length (cm)",
       y = paste0("PC1 (", round(Pvar[Pvar$PC == "PC1", "p.var"], 2), "%)"))+
  scale_y_continuous(breaks = 0)+
  theme_classic()+
  theme(panel.grid.major.y = element_line(colour = "grey50", linewidth = 0.1),
        plot.margin = unit(c(1,1,1,1), "mm"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, angle = 90),
        plot.title = element_text(size = 18, face = "bold"))


### B: PC1 vs. PC2
p1b <- ggplot(PCs, aes(x = PC1, y = PC2))+
  geom_point(colour = "grey50", alpha = 0.4, size = 0.8)+
  geom_point(data = PCs[PCs$Ecotype == "Wave", ], colour = "#fc8d59", alpha = 0.3)+
  geom_point(data = PCs[PCs$Ecotype == "Crab", ], colour = "#91bfdb", alpha = 0.3)+
  geom_point(data = PCs[PCs$Site == Wave, ], colour = "#d73027", pch = 15)+
  geom_point(data = PCs[PCs$Site == Crab, ], colour = "#4585d9", pch = 15)+
  labs(x = paste0("PC1 (", round(Pvar[Pvar$PC == "PC1", "p.var"], 2), "%)"),
       y = paste0("PC2 (", round(Pvar[Pvar$PC == "PC2", "p.var"], 2), "%)"))+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "grey50", linewidth = 0.1),
        plot.margin = unit(c(1,1,1,1), "mm"),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, angle = 90),
        plot.title = element_text(size = 18, face = "bold"))

### C: PC1 vs. PC3
p1c <- ggplot(PCs, aes(x = PC1, y = PC3))+
  geom_point(colour = "grey50", alpha = 0.4, size = 0.8)+
  geom_point(data = PCs[PCs$Ecotype == "Wave", ], colour = "#fc8d59", alpha = 0.3)+
  geom_point(data = PCs[PCs$Ecotype == "Crab", ], colour = "#91bfdb", alpha = 0.3)+
  geom_point(data = PCs[PCs$Site == Wave, ], colour = "#d73027", pch = 15)+
  geom_point(data = PCs[PCs$Site == Crab, ], colour = "#4585d9", pch = 15)+
  labs(x = paste0("PC1 (", round(Pvar[Pvar$PC == "PC1", "p.var"], 2), "%)"),
       y = paste0("PC3 (", round(Pvar[Pvar$PC == "PC3", "p.var"], 2), "%)"))+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "grey50", linewidth = 0.1),
        plot.margin = unit(c(1,1,1,1), "mm"),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, angle = 90),
        plot.title = element_text(size = 18, face = "bold"))

### D: PC2 vs. PC3
p1d <- ggplot(PCs, aes(x = PC2, y = PC3))+
  geom_point(colour = "grey50", alpha = 0.4, size = 0.8)+
  geom_point(data = PCs[PCs$Ecotype == "Wave", ], colour = "#fc8d59", alpha = 0.3)+
  geom_point(data = PCs[PCs$Ecotype == "Crab", ], colour = "#91bfdb", alpha = 0.3)+
  geom_point(data = PCs[PCs$Site == Wave, ], colour = "#d73027", pch = 15)+
  geom_point(data = PCs[PCs$Site == Crab, ], colour = "#4585d9", pch = 15)+
  labs(x = paste0("PC2 (", round(Pvar[Pvar$PC == "PC2", "p.var"], 2), "%)"),
       y = paste0("PC3 (", round(Pvar[Pvar$PC == "PC3", "p.var"], 2), "%)"))+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "grey50", linewidth = 0.1),
        plot.margin = unit(c(1,1,1,1), "mm"),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, angle = 90),
        plot.title = element_text(size = 18, face = "bold"))

### Multipanel plot
# 1200 x 300
ggarrange(p1a, p1b, p1c, p1d)

