############# Perform PCA on Phynotypic data ######################
### This script performs a simple principle component analysis on
### phenotypic measurements from 1,488 Littorina saxatilis that 
### were collected around Tjärnö, Sweden in the summer of 2020.
### These measurement come from two different sources, i) dissection
### records and ii) photos of snail shells that were analysed by 
### Jenny Larsson's ShellShaper program in MatLab (https://github.com/jslarsson/ShellShaper)
### Version 5 of this scripts separates shell length from shape measurements
### and adjust the ShellShaper paramteres by shell lengths.
### James Reeve - DTU Aqua
### 16/09/2025

### Measurements ####
### Dissection records
#   i. "site" = name of sampling site
#   ii. "snail_ID" = ID for individual snail; a=lowshore to h=highshore
#   iii. "avg_thick_mm" = average shell thickness (mm)
#   iv. "sex" = Trinary score; M=Male, F=Female and J=Juvenile (i.e. undeterminable)
#   v. "adult" = Life history stage (TRUE/FASLE)
#   vi. "parasite = presence of parasites (TRUE/FASLE)

### ShellShaper parameters
#   i. "snail_ID" = as above
#   ii. "gw" = Growth factor for shell width
#   iii. "gh" = Growth factor for shell height
#   iv. "r0" = Sprial radius (mm) - centre of apeture to centre of spiral
#   v. "z0" = Spiral height (mm) - apeture centre to apex
#   vi. "a0" = Apeture size (mm) - centre to right edge
#   vii. "eccentricity" = how circular the apeture is 0=circle; >0 ovoid
#   viii. "apAngle" = Angle of apeture opening
#   ix. "shellLength" = Length from apeture end to apex (mm)
#   x. "apex_x" = horizontal position of the apex (pixels) [ignored]
#   xi. "apex_y" = verticale position of the apex (pixels) [ignored]
#   xii. "scaleFactor" = conversion facotr to turn pixels to mm [ignored]


### Preparation ####
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)
if(Sys.info()['sysname'] == "Windows") {
  setwd("C:/Users/Windows/Folder")
} else {
  setwd("/Users/OSX/directory")
}


### Packages
library(tidyverse)
library(ggpubr)
library(ggcorrplot)


#### A: Access and wrangle data ####

### Access data
PATH <- "/path/to/data/"
# Dissection records
Diss <- read.csv(paste0(PATH, "KT_dissection_log.csv"))
# ShellShaper parameters
ShS <- read.csv(paste0(PATH, "ShellShaper/KT_ShellShaper_parameters.v2.csv"))

### Remove "_ac" from ShS$snailID
ShS$snailID <- gsub("_ac.*", "", ShS$snailID)

### Remove unnecessary rows
Diss <- Diss[,-which(colnames(Diss) %in% c("site", "date", "thick1", "thick2", "thick3"))]
ShS <- ShS[,-which(colnames(ShS) %in% c("apex_x", "apex_y", "scaleFactor"))]

### Convert sex == 'J' to 'F'
# We are assuming all juvenile snails were female, since we identified juvenile male
Diss[Diss$sex == "J", "sex"] <- "F"

### Re-score adult as c("adult", "juvenile")
Diss$adult <- ifelse(Diss$adult, "adult", "juvenile")

### Log transform growth parameters (gw & gh)
ShS$gw <- log(ShS$gw)
ShS$gh <- log(ShS$gh)

### Merge data
Pheno <- merge(ShS, Diss, by.x = "snailID", by.y = "snail_ID")
# Note: several snails were photographed, but then found dead upon dissection

### Adjust length measurements by shellLength
### This accounts for the correlations between these measurements, and should
### decouple the effects of shell shape from snail size.
Pheno[, which(colnames(Pheno) %in% c("r0", "z0", "a0", "avg_thick_mm"))] <- 
  apply(Pheno[, which(colnames(Pheno) %in% c("r0", "z0", "a0", "avg_thick_mm"))], 
        2, function(X){X / Pheno$shellLength})

# Check adjustment with Corrplot
ggcorrplot(cor(Pheno %>% select(where(is.numeric))), type = "upper",
           p.mat = cor_pmat(Pheno %>% select(where(is.numeric))),
           title = "Littorina saxatilis phenotypes")

### Adjust eccentricity by aperture size (a0)
### This is a correction recommended by Jenny Larsson
### It scales the eccentricity measurement by the radius of the aperture
Pheno$eccentricity <- Pheno$eccentricity / Pheno$a0  



#### B: Run the PCA ####
### Drop ShellLength
Pheno2 <- Pheno[,which(colnames(Pheno) != "shellLength")]

### Scale and centre continuous data
tmp <- scale(Pheno2[ ,which(sapply(Pheno2, class) == "numeric")], center = TRUE, scale = TRUE)

# Add sample as rownames
rownames(tmp) <- Pheno2$snailID

# Run the PCA
Pheno.PCA <- prcomp(tmp, center = FALSE)

# Get PC values
PCs <- as.data.frame("Snail_ID" = rownames(Pheno.PCA$x),
                     Pheno.PCA$x)

# Get variance explained per axis
Pvar <- data.frame("PC" = colnames(Pheno.PCA$rotation),
                   "p.var" = (Pheno.PCA$sdev^2) / sum(Pheno.PCA$sdev^2) * 100)

# Extract variable loadings
PCAloadings <- data.frame("Vars" = rownames(Pheno.PCA$rotation), 
                          Pheno.PCA$rotation)



#### C: Plots ####

### PC1 vs PC2 scatterplot
# Copy dimensions 460 x 400
ggplot() +
  geom_point(data = PCs, aes(x = PC1, y = PC2, colour = Pheno2$sex), 
             size = 2, alpha = 0.4, show.legend = FALSE)+
  geom_segment(data = PCAloadings, 
               aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "grey30")+
  geom_label(data = PCAloadings, aes(x = PC1*5, y = PC2*5, label = Vars),
             size = 2.6)+
  labs(x = paste0("PC1 (", round(Pvar[Pvar$PC == "PC1", "p.var"], 2), "%)"),
       y = paste0("PC2 (", round(Pvar[Pvar$PC == "PC2", "p.var"], 2), "%)"),
       title = "Littorina saxatilis phenotypes")+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  scale_colour_manual(values = c("#FF55DD", "#2A7FFF"))+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "grey50", linewidth = 0.1),
        plot.margin = unit(c(1,1,1,1), "mm"),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, angle = 90),
        plot.title = element_text(size = 18, face = "bold"))


### Screeplot 
# Copy ratio: 400 x 400
ggplot(Pvar, aes(x = as.numeric(gsub("PC", "", PC)), y = p.var))+
  geom_line()+
  geom_point()+
  labs(x = "Principal component", y = "Percentage varaince explained")+
  ylim(c(0, 100))+
  scale_x_continuous(breaks = 1:nrow(Pvar))+
  theme_classic()


### Histogram of PC1-5
# Copy ratio: 300 X 500
PCA.hist <- function(data, principal.compent, binwidth = 0.2,
                     highlight.colour = "grey50"){
  Xlim <- c(floor(min(sapply(data, range))),
            ceiling(max(sapply(data, range)))) # Set plot limits based on max and min PC scores
  Xbreaks <- seq(Xlim[1], Xlim[2], 1)               # Set breaks based on 'Xlim'
  
  p <- ggplot(data, aes(x = data[, principal.compent]))+
    geom_histogram(binwidth = binwidth, colour = "grey50", fill = NA, alpha = 0.4)+
    geom_density(aes(y = binwidth * ..count..), alpha = 0.4, col = highlight.colour)+
    labs(y = paste0("PC", principal.compent))+
    scale_x_continuous(breaks = Xbreaks, limits = Xlim)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return(p)
}

ggarrange(PCA.hist(PCs, principal.compent = 1, binwidth = 0.4,
                   highlight.colour = "#888EEB"),
          PCA.hist(PCs, principal.compent = 2, binwidth = 0.4,
                   highlight.colour = "#74EBA3"),
          PCA.hist(PCs, principal.compent = 3, binwidth = 0.4,
                   highlight.colour = "#D48CE6"),
          PCA.hist(PCs, principal.compent = 4, binwidth = 0.4,
                   highlight.colour = "#E6D38A"),
          PCA.hist(PCs, principal.compent = 5, binwidth = 0.4,
                   highlight.colour = "#EB817A"),
          ncol = 1)


### PCA scatter plot PC1-8
# Copy ratio: 560 x 560
PCA.Scatterplot <- function(i, j){
  ggplot(PCs, aes(x = PCs[,i], y = PCs[,j], colour = Pheno2$sex))+
    geom_point(size = 0.2, alpha = 0.4, stroke = 1.2, show.legend = FALSE)+
    labs(x = paste0("PC", i,"(", round(Pvar$p.var[i], 2), "%)"),
         y = paste0("PC", j, "(", round(Pvar$p.var[j], 2), "%)"),
         colour = "Sex")+
    scale_x_continuous(breaks = 0)+
    scale_y_continuous(breaks = 0)+
    scale_colour_manual(values = c("#FF55DD", "#2A7FFF"))+ # Colour palette for sex
    #scale_colour_manual(values = c("#7570b3", "#1b9e77"))+ # Colour palette for maturity
    #scale_colour_manual(values = c("grey50", "red"))+ # Colour palette for parasites
    theme_void()+
    theme(panel.grid.major = element_line(colour = "grey50", linewidth = 0.1),
          plot.margin = unit(c(1,1,1,1), "mm"),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6, angle = 90))
}

p_ <- ggplot()+theme_void()
ggarrange(
          # Row1: PC1 vs. PC2-8
          ggarrange(PCA.Scatterplot(1,2), PCA.Scatterplot(1,3), PCA.Scatterplot(1,4),
                    PCA.Scatterplot(1,5), PCA.Scatterplot(1,6), PCA.Scatterplot(1,7),
                    PCA.Scatterplot(1,8), nrow = 1),
          # Row2: PC2 vs. PC3-8
          ggarrange(p_, PCA.Scatterplot(2,3), PCA.Scatterplot(2,4), 
                    PCA.Scatterplot(2,5), PCA.Scatterplot(2,6), PCA.Scatterplot(2,7),
                    PCA.Scatterplot(2,8), nrow = 1),
          # Row3: PC3 vs. PC4-8
          ggarrange(p_, p_, PCA.Scatterplot(3,4), 
                    PCA.Scatterplot(3,5), PCA.Scatterplot(3,6), PCA.Scatterplot(3,7),
                    PCA.Scatterplot(3,8), nrow = 1),
          # Row4: PC4 vs. PC5-8
          ggarrange(p_, p_, p_, 
                    PCA.Scatterplot(4,5), PCA.Scatterplot(4,6), PCA.Scatterplot(4,7),
                    PCA.Scatterplot(4,8), nrow = 1),
          # Row5: PC5 vs. PC6-8
          ggarrange(p_, p_, p_, 
                    p_, PCA.Scatterplot(5,6), PCA.Scatterplot(5,7),
                    PCA.Scatterplot(5,8), nrow = 1),
          # Row6: PC6 vs. PC7-8
          ggarrange(p_, p_, p_, p_, p_, PCA.Scatterplot(6,7), 
                    PCA.Scatterplot(6,8), nrow = 1),
          # Row7: PC7 vs. PC8
          ggarrange(p_, p_, p_, p_, p_, p_, 
                    PCA.Scatterplot(7,8), nrow = 1),
          ncol = 1)


### Factor loading plots
# Copy ratio: 1000 x 250
PCA.load <- function(data, principal.component, highlight.colour = "grey50"){
  ggplot(data)+
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.1)+
    geom_segment(aes(x = 0, y = Vars, xend = {{principal.component}}, yend = Vars), 
                 arrow = arrow(length = unit(0.8, "picas")),
                 color = highlight.colour, linewidth = 1.2)+
    xlim(c(-1, 1))+
    labs(x = deparse(substitute(principal.component)))+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())
}

ggarrange(PCA.load(PCAloadings, principal.component = PC1, highlight.colour = "#888EEB"),
          PCA.load(PCAloadings, principal.component = PC2, highlight.colour = "#74EBA3"),
          PCA.load(PCAloadings, principal.component = PC3, highlight.colour = "#D48CE6"),
          PCA.load(PCAloadings, principal.component = PC4, highlight.colour = "#E6D38A"),
          PCA.load(PCAloadings, principal.component = PC5, highlight.colour = "#EB817A"),
          nrow = 1)



#### D: Determine number of PCs to retain ####

### Broken stick model
Pvar$BSexp <- sapply(1:nrow(Pvar), function(i){ sum(1/i:nrow(Pvar))/nrow(Pvar) * 100})

ggplot(Pvar, aes(x = PC))+
  geom_point(aes(y = p.var), col = "black")+
  geom_line(aes(y = p.var, group = 1), col = "black")+
  geom_point(aes(y = BSexp), col = "grey80")+
  geom_line(aes(y = BSexp, group = 1), col = "grey80")+
  labs(x = "Principal component", y = "Variance explained (%)")+
  theme_bw()

### Inspect hard cuts
Pvar$cumm.var <- Pvar$p.var[1]
for(i in 2:nrow(Pvar)){
  Pvar$cumm.var[i] <- Pvar$p.var[i] + Pvar$cumm.var[i-1]
}
# 75% = PC1-3
# 80% = PC1-4
# 90% = PC1-5



#### E: Save PCA ####
OUT <- "/output/filepath/"

# Save principal components
write.csv(PCs, paste0(OUT, "principal_components/Pheno_PC1-8.v3.csv"), 
          quote = FALSE)

# Save eigenvalues and percent of variance explained
write.csv(Pvar, paste0(OUT, "eigenvalues/Pheno_variance_explained.v3.csv"), quote = FALSE)
