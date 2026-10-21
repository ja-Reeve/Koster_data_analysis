####################### Environmental variance among habitats ######################
### Running PCAmix on envrionmental data, excluding any measurements which were 
### specific to a single habitat. The variance explained by each PC will then be 
### contrasted among habitats. Observations in the field, and plotting varaibles 
### independently, suggest that rocky shores should be more varaible.

### James Reeve - University of Gothenburg
### 2025-04-25

### Preparation
rm(list = ls())
try(dev.off())
setwd("/Users/james/")

### Pacakges
library(PCAmixdata)
library(ggplot2)
library(dplyr)

### Filepath
PATH <- "file/path/"


#### 1. Access data ####

### A: Field notes for each site - after adding sea level data
FB_Site <- read.csv(paste0(PATH, "KT_fieldnotes.v2.csv"))

### B: Drone photo measurements
DPh <- read.csv(paste0(PATH, "KT_drone_photos_v2.csv"))

### C: GiS data
GisD <- read.csv(paste0(PATH, "KT_GIS_data_v2.csv"))

### D: Phenotypic PCs
Pheno <- read.csv("file/path/Phenotypes/Pheno_PC1-8.v3.csv",
                  col.names = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8"))



#### 2. Data Wrangling ####

### A: Remove unused variables
FB_Site <- FB_Site[,c(1,7,14:20,24)]
DPh <- DPh[,c(1:ncol(DPh)-1)]

### B: Merge date
tmp <- merge(FB_Site, DPh, by= "Site") # Note: "no snail" sites are removed here
Env <- merge(tmp, GisD, by = "Site")
rm(tmp)

### C: Calculate shore slope
# shore height = field notes
# shore length = drone photo measurement
# slope = arc tanget (shore height / shore length) [radians]
Env$slope <- atan(Env$Transect_height / Env$Shore_length)


### D: Adjust drone photo measurements by shore slope
slope.adjust <- function(X){X/cos(Env$slope)}
Env$Adj.FW <- slope.adjust(Env$FW)
Env$Adj.BZW <- slope.adjust(Env$BZW)
Env$Adj.TZW <- slope.adjust(Env$TZW)
Env$Adj.LVD <- slope.adjust(Env$LVD)

# Manual adjustment: revert change for cliff sites (KT195 & KT196), as a slope
# of 90° gives ridiculously high adjusted scores.
Env$Adj.FW[195:196] <- Env$FW[195:196]
Env$Adj.BZW[195:196] <- Env$BZW[195:196] 
Env$Adj.TZW[195:196] <- Env$TZW[195:196]
Env$Adj.LVD[195:196] <- Env$LVD[195:196]


### E: Adjust aspect to make 0 and 359 close
# Roger suggested to adjust aspect as the sin and cos of aspect (02/03/2023)
Env$cos_aspect <- cos(pi*Env$Aspect_shore/180)
Env$sin_aspect <- sin(pi*Env$Aspect_shore/180)


### F: Infer missing field note measurements of Ascophyllum and Barnacles using drone photos
# Replace '?' with 'Un'
Env[Env$Ascophyllum %in% c("???", "????", "?????"), "Ascophyllum"] <- "Un"
Env[Env$Barnacles == "????", "Barnacles"] <- "Un"

# Replace "Un" with values from drone photos
Env$Ascophyllum <- ifelse(Env$Ascophyllum == "Un", Env$Asc, Env$Ascophyllum)
Env$Barnacles <- ifelse(Env$Barnacles == "Un", Env$Barn, Env$Barnacles)


### G: Select variables to keep
keeps <- c("Site", "Habitat.y", "Exposure_FETCH", "sin_aspect", "cos_aspect", "slope", 
           "Shore_length", "Fucus", "Ascophyllum", "Barnacles", "Boulders_rocks_stones", 
           "HI", "Adj.FW", "FP", "Adj.LVD", "HTD")
Env <- Env[, keeps]

### Fix name of "Habitat"
colnames(Env)[2] <- "Habitat"


#### 3. Run PCAmix ####

### A: Function to run PCAmix
ShorePCA <- function(dat){
  # Convert site names to rownames 
  rownames(dat) <- dat[,"Site"]
  dat$Site <- NULL
  
  # Drop "habitat while running PCA
  dat$Habitat <- NULL
  
  # Separate out continuous and categorical variables
  tmp <- splitmix(dat)
  qual <- tmp$X.quali
  quant <- tmp$X.quanti
  
  # Rescale and recentre quantitative variables
  quant <- scale(quant, center = TRUE, scale = TRUE)
  
  # Run the PCA
  tmp_PCA <- PCAmix(X.quanti = quant, X.quali = qual,
                    rename.level = TRUE, graph = FALSE)
  
  return(tmp_PCA)
}


### B: Run the PCA
EnvPCA <- ShorePCA(Env)


### C: Default PCA plots
plot(EnvPCA, axes = c(1,2), choice = "ind", label = TRUE)
plot(EnvPCA, choice = "levels")
plot(EnvPCA, choice = "cor")
plot(EnvPCA, axes = c(1,2), choice = "sqload", coloring.var = T, 
               leg=TRUE, posleg="none")


### D: Custom PCA plots with ggplot() 
# Copy to clipboard: 800x400 (w x h)

# Create dataframe
PCs <- as.data.frame(EnvPCA$ind$coord)
PCs$Site <- rownames(PCs)

# Add habitat to scores
PCs$Habitat <- Env$Habitat[Env$Site %in% PCs$Site]

# Plot
p.PCA <- ggplot()+
  geom_point(data = PCs, aes(`dim 1`, `dim 2`, colour = Habitat))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  geom_segment(data = Loads, aes(x = 0, xend = `dim 1`*8,
                   y = 0, yend = `dim 2`*8), colour = "grey25")+
  geom_label(data = Loads, aes(x = `dim 1`*8, y = `dim 2`*8, label = Var, fill = Class), 
             size = 3, colour = "white", show.legend = FALSE)+
  scale_fill_manual(values = c("firebrick", "grey50"))+
  labs(x = paste0("Dim 1 (", round(EnvPCA$eig[1,"Proportion"], 2), "%)"), 
       y = paste0("Dim 2 (", round(EnvPCA$eig[2,"Proportion"], 2), "%)"), 
       title = "All habitats")+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom")
p.PCA


### E: Custom factor loading plot with ggplot

Loads <- as.data.frame(EnvPCA$sqload)
Loads$Var <- rownames(Loads)
Loads$Class <- "continous"
Loads$Class[Loads$Var %in% c("Fucus", "Ascophyllum", "Barnacles")] <- "catagorical"

p.Loads <- ggplot(Loads)+
  geom_segment(aes(x = 0, xend = `dim 1`,
                   y = 0, yend = `dim 2`,
                   colour = Class))+
  geom_label(aes(x = `dim 1`, y = `dim 2`, label = Var, colour = Class),
             size = 3)+
  labs(x = paste0("Dim 1 (", round(EnvPCA$eig[1,"Proportion"], 2), "%)"), 
       y = paste0("Dim 2 (", round(EnvPCA$eig[2,"Proportion"], 2), "%)"), 
       title = "Factor loadings")+
  scale_colour_manual(values = c("firebrick", "dodgerblue"))+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"), 
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.85))



#### 4. Plot varaince among habitats ####

### A: Calcaulte variance of each habitat
vars <- PCs %>% pivot_longer(cols = -c(Site, Habitat), 
                             names_to = "dimension", 
                             values_to = "score") %>%
  group_by(Habitat, dimension) %>%
  summarise("var" = var(score),
            "mean" = mean(score))


### B: Simple column plot
# Copy to clipboard: 1000x250 (w x h)
p.var <- ggplot(vars)+
  geom_col(aes(x = Habitat, y = var, fill = Habitat), show.legend = FALSE)+
  facet_wrap(vars(dimension), nrow = 1)+
  labs(x = "", y = "Variance in Envrionment")+
  scale_fill_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30))


### C: More detailed violin plots
# Copy to clipboard: 1000x300 (w x h)
p1 <- ggplot(PCs)+
  geom_violin(aes(x = Habitat, y = `dim 1`, fill = Habitat, colour = Habitat),
              alpha = 0.4)+
  geom_jitter(aes(x = Habitat, y = `dim 1`, colour = Habitat), 
              width = 0.2, height = 0, size = 0.6, alpha = 0.7)+
  geom_point(data = vars[vars$dimension == "dim 1",], 
             aes(x = Habitat, y = mean, colour = Habitat))+
  labs(x = "", y = "Envrionmental PCs", title = "Dim 1")+
  scale_fill_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "none")

p2 <- ggplot(PCs)+
  geom_violin(aes(x = Habitat, y = `dim 2`, fill = Habitat, colour = Habitat),
              alpha = 0.4)+
  geom_jitter(aes(x = Habitat, y = `dim 2`, colour = Habitat), 
              width = 0.2, height = 0, size = 0.6, alpha = 0.7)+
  geom_point(data = vars[vars$dimension == "dim 2",], 
             aes(x = Habitat, y = mean, colour = Habitat))+
  labs(x = "", y = "", title = "Dim 2")+
  scale_fill_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "none")

p3 <- ggplot(PCs)+
  geom_violin(aes(x = Habitat, y = `dim 3`, fill = Habitat, colour = Habitat),
              alpha = 0.4)+
  geom_jitter(aes(x = Habitat, y = `dim 3`, colour = Habitat), 
              width = 0.2, height = 0, size = 0.6, alpha = 0.7)+
  geom_point(data = vars[vars$dimension == "dim 3",], 
             aes(x = Habitat, y = mean, colour = Habitat))+
  labs(x = "", y = "", title = "Dim 3")+
  scale_fill_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "none")

ggarrange(p1, p2, p3, nrow = 1)


### D: Multipanel plot
ggarrange(
  ggarrange(p.PCA, p.Loads, nrow = 1),
  p.var,
  ggarrange(p1, p2, p3, nrow = 1),
  nrow = 3, heights = c(1, 0.4, 0.4)
)



#### 5. Plot phenotypic variance among habitats ####

### A: Data wrangle
# Add site ID to phenotypic data
Pheno$Site <- substr(Pheno$ID, 1, 5)

# Merge Phenotypic and envrionmental PCA
PhenoEnv <- merge(Pheno, PCs, by = "Site")

# Mean and Variance of each phenotypic PC
vars.pheno <- PhenoEnv %>% 
  pivot_longer(cols = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9),
               values_to = "score", names_to = "PC") %>%
  group_by(Habitat, PC) %>%
  summarise("var" = var(score),
            "mean" = mean(score))


### B: Violin plots of Phenotypic PCs among habitats
p1 <- ggplot(PhenoEnv)+
  geom_violin(aes(x = Habitat, y = PC1, fill = Habitat, colour = Habitat),
              alpha = 0.4)+
  geom_jitter(aes(x = Habitat, y = PC1, colour = Habitat), 
              width = 0.2, height = 0, size = 0.6, alpha = 0.6)+
  geom_point(data = vars.pheno[vars.pheno$PC == "PC1",], 
             aes(x = Habitat, y = mean, colour = Habitat), size = 3)+
  labs(x = "", y = "Phenotypic PCs", title = "PC1")+
  scale_fill_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "none")


### C: Scatterplots comparing environmental and phenotypic PCs
p.E1P1 <- ggplot(PhenoEnv, aes(x = `dim 1`, y = PC1, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 1", y = expression(PC1[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E1P2 <- ggplot(PhenoEnv, aes(x = `dim 1`, y = PC2, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 1", y = expression(PC2[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E1P3 <- ggplot(PhenoEnv, aes(x = `dim 1`, y = PC3, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 1", y = expression(PC3[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()
  
p.E2P1 <- ggplot(PhenoEnv, aes(x = `dim 2`, y = PC1, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 2", y = expression(PC1[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E2P2 <- ggplot(PhenoEnv, aes(x = `dim 2`, y = PC2, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 2", y = expression(PC2[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E2P3 <- ggplot(PhenoEnv, aes(x = `dim 2`, y = PC3, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 2", y = expression(PC3[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E3P1 <- ggplot(PhenoEnv, aes(x = `dim 3`, y = PC1, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 3", y = expression(PC1[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E3P2 <- ggplot(PhenoEnv, aes(x = `dim 3`, y = PC2, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 3", y = expression(PC2[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

p.E3P3 <- ggplot(PhenoEnv, aes(x = `dim 3`, y = PC3, group = Habitat, colour = Habitat))+
  geom_point(size = 0.8, alpha = 0.6)+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Environmental Dimension 3", y = expression(PC3[Pheno]))+
  scale_colour_manual(values = c("#a2c4c9ff", "#b4a7d6ff", "#ea9999ff"))+
  theme_classic()

ggarrange(p.E1P1, p.E2P1, p.E3P1,
          p.E1P2, p.E2P2, p.E3P2,
          p.E1P3, p.E2P3, p.E3P3,
          nrow = 3, ncol = 3, common.legend = TRUE)
