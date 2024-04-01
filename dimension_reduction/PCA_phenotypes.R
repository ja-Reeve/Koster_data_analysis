############# Perform PCA on Phynotypic data ######################
### This script performs a simple principle component analysis on
### phenotypic measurements from 1,488 Littorina saxatilis that 
### were collected around Tjärnö, Sweden in the summer of 2020.
### These measurement come from two different sources, i) dissection
### records and ii) photos of snail shells that were analysed by 
### Jenny Larsson's ShellShaper program in MatLab (https://github.com/jslarsson/ShellShaper)
### Version 4 of this scripts drops PCAmix to focus just on the phenotypic PCs
### James Reeve - University of Gothenburg
### 01/09/2023

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


### Packages
library(tidyverse)
#library(ggfortify)
library(ggpubr)
library(ggcorrplot)
#library(PCAmixdata)
library(cluster)


### Access data
# Filepath
PATH <- "/input/dir/"
# Dissection records
Diss <- read.csv(paste0(PATH, "KT_dissection_log.csv"))
# ShellShaper parameters
ShS <- read.csv(paste0(PATH, "KT_ShellShaper_parameters.v2.csv"))


#### A: Format and merge data ####

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


#### B: Explore the data ####

### B1: Categorical variables
# Counts
table(Pheno$sex) # Sex ratio
table(Pheno$adult) # Maturity
table(Pheno$parasite) # Presence of parasites

# Plots
p1 <- ggplot(Pheno)+
  geom_bar(aes(x = sex), fill = "#9dbfc4ff")+
  labs(title = "Sex")+
  ylim(c(0,1410))+
  theme_bw()+
  theme(axis.title = element_blank())

p2 <- ggplot(Pheno)+
  geom_bar(aes(x = adult), fill = "#9dbfc4ff")+
  labs(title = "Adult")+
  ylim(c(0,1410))+
  theme_bw()+
  theme(axis.title = element_blank())

p3 <- ggplot(Pheno)+
  geom_bar(aes(x = parasite), fill = "#9dbfc4ff")+
  labs(title = "Paracitised")+
  ylim(c(0,1410))+
  theme_bw()+
  theme(axis.title = element_blank())

ggarrange(p1,p2,p3, nrow = 1)


### B2: Continuous variables
# Distribution summaries
sapply(Pheno %>% select(where(is.numeric)), summary)

# Violin plots function
Violin.Plot <- function(x, y, axis.labels){
  if(!(class({{x}}) %in% c("character", "factor"))){stop("Error: 'x' must be a catagorical variable or factor")}
  if(!(class({{y}}) %in% c("integer", "numeric"))){stop("Error: 'y' must be numeric or an iteger vector")}
  if(length(axis.labels) != 2){stop("Error: 'axis.labels' must be a vector of two characters c(x, y)")}
  
  # Calculate mean
  tmp <- Pheno %>% group_by({{x}}) %>% summarise("Var" = unique({{x}}), "Avg" = mean(y))
  
  # Violin plot
  p <- ggplot(Pheno, aes({{x}}, {{y}}, colour = {{x}}))+
    geom_jitter(size = 0.2, width = 0.1, height = 0)+
    geom_violin(aes(fill = {{x}}), alpha = 0.2, colour = "black", linewidth = 0.1)+
    geom_point(data = tmp, aes(Var, Avg), colour = "grey20", pch = 18, size = 4)+
    labs(x = axis.labels[1], y = axis.labels[2])+
    scale_colour_manual(values = c("#FF55DD", "#2A7FFF"))+
    theme_classic()+
    theme(legend.position = "none")
  
  return(p)
}

# Violin plots by sex
ggarrange(plotlist = list(
  Violin.Plot(sex, shellLength, c("", "Shell length (mm)")),
  Violin.Plot(sex, avg_thick_mm, c("", "Shell thickness (mm)")),
  Violin.Plot(sex, apAngle, c("", "Apeture angle (°)")),
  Violin.Plot(sex, eccentricity, c("", "Eccentricity")),
  Violin.Plot(sex, gh, c("", "Growth height")),
  Violin.Plot(sex, gw, c("", "Growth width")),
  Violin.Plot(sex, r0, c("", "Spiral radius (mm)")), 
  Violin.Plot(sex, z0, c("", "Spiral height (mm)")),
  Violin.Plot(sex, a0, c("", "Apeture size (mm)"))),
  nrow = 3, ncol = 3)

# Violin plots by maturity
ggarrange(plotlist = list(
  Violin.Plot(adult, shellLength, c("", "Shell length (mm)")),
  Violin.Plot(adult, avg_thick_mm, c("", "Shell thickness (mm)")),
  Violin.Plot(adult, apAngle, c("", "Apeture angle (°)")),
  Violin.Plot(adult, eccentricity, c("", "Eccentricity")),
  Violin.Plot(adult, gh, c("", "Growth height")),
  Violin.Plot(adult, gw, c("", "Growth width")),
  Violin.Plot(adult, r0, c("", "Spiral radius (mm)")), 
  Violin.Plot(adult, z0, c("", "Spiral height (mm)")),
  Violin.Plot(adult, a0, c("", "Apeture size (mm)"))),
  nrow = 3, ncol = 3)


# Pairwise scatterplots
Pairwise.Scatterplot <- function(x, y, axis.labels){
  if(!(class({{x}}) %in% c("integer", "numeric"))){stop("Error: 'x' must be numeric or an iteger vector")}
  if(!(class({{y}}) %in% c("integer", "numeric"))){stop("Error: 'y' must be numeric or an iteger vector")}
  if(length(axis.labels) != 2){stop("Error: 'axis.labels' must be a vector of two characters c(x, y)")}
  
  ggplot(Pheno, aes({{x}}, {{y}}))+
    geom_point(size = 0.6, alpha = 0.4)+
    labs(x = axis.labels[1], y = axis.labels[2])+
    theme_classic()+
    theme(axis.title = element_text(size = 6))
}

ggarrange(
  ggarrange(plotlist = list(Pairwise.Scatterplot(shellLength, avg_thick_mm, c("Shell length (mm)", "Shell thickness (mm)")),
                            Pairwise.Scatterplot(shellLength, apAngle, c("Shell length (mm)", "Apeture angle (°)")),
                            Pairwise.Scatterplot(shellLength, eccentricity, c("Shell length (mm)", "Eccentricity")),
                            Pairwise.Scatterplot(shellLength, gh, c("Shell length (mm)", "Growth height")),
                            Pairwise.Scatterplot(shellLength, gw, c("Shell length (mm)", "Growth width")),
                            Pairwise.Scatterplot(shellLength, r0, c("Shell length (mm)", "Spiral radius (mm)")),
                            Pairwise.Scatterplot(shellLength, z0, c("Shell length (mm)", "Spiral height (mm)")),
                            Pairwise.Scatterplot(shellLength, a0, c("Shell length (mm)", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            Pairwise.Scatterplot(avg_thick_mm, apAngle, c("Shell thickness (mm)", "Apeture angle (°)")),
                            Pairwise.Scatterplot(avg_thick_mm, eccentricity, c("Shell thickness (mm)", "Eccentricity")),
                            Pairwise.Scatterplot(avg_thick_mm, gh, c("Shell thickness (mm)", "Growth height")),
                            Pairwise.Scatterplot(avg_thick_mm, gw, c("Shell thickness (mm)", "Growth width")),
                            Pairwise.Scatterplot(avg_thick_mm, r0, c("Shell thickness (mm)", "Spiral radius (mm)")),
                            Pairwise.Scatterplot(avg_thick_mm, z0, c("Shell thickness (mm)", "Spiral height (mm)")),
                            Pairwise.Scatterplot(avg_thick_mm, a0, c("Shell thickness (mm)", "Apeture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            Pairwise.Scatterplot(apAngle, eccentricity, c("Apeture angle (°)", "Eccentricity")),
                            Pairwise.Scatterplot(apAngle, gh, c("Apeture angle (°)", "Growth height")),
                            Pairwise.Scatterplot(apAngle, gw, c("Apeture angle (°)", "Growth width")),
                            Pairwise.Scatterplot(apAngle, r0, c("Apeture angle (°)", "Spiral radius (mm)")),
                            Pairwise.Scatterplot(apAngle, z0, c("Apeture angle (°)", "Spiral height (mm)")),
                            Pairwise.Scatterplot(apAngle, a0, c("Apeture angle (°)", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            Pairwise.Scatterplot(eccentricity, gh, c("Eccentricity", "Growth height")),
                            Pairwise.Scatterplot(eccentricity, gw, c("Eccentricity", "Growth width")),
                            Pairwise.Scatterplot(eccentricity, r0, c("Eccentricity", "Spiral radius (mm)")),
                            Pairwise.Scatterplot(eccentricity, z0, c("Eccentricity", "Spiral height (mm)")),
                            Pairwise.Scatterplot(eccentricity, a0, c("Eccentricity", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            Pairwise.Scatterplot(gh, gw, c("Growth height", "Growth width")),
                            Pairwise.Scatterplot(gh, r0, c("Growth height", "Spiral radius (mm)")),
                            Pairwise.Scatterplot(gh, z0, c("Growth height", "Spiral height (mm)")),
                            Pairwise.Scatterplot(gh, a0, c("Growth height", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            Pairwise.Scatterplot(gw, r0, c("Growth width", "Spiral radius (mm)")),
                            Pairwise.Scatterplot(gw, z0, c("Growth width", "Spiral height (mm)")),
                            Pairwise.Scatterplot(gw, a0, c("Growth width", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            Pairwise.Scatterplot(r0, z0, c("Spiral radius (mm)", "Spiral height (mm)")),
                            Pairwise.Scatterplot(r0, a0, c("Spiral radius (mm)", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
  ggarrange(plotlist = list(ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            ggplot()+theme_void(),
                            Pairwise.Scatterplot(z0, a0, c("Spiral height (mm)", "Aperture size (mm)"))), 
            nrow = 1, ncol = 8),
nrow = 8, ncol = 1)

# Correlation heatmap among variables
ggcorrplot(cor(Pheno %>% select(where(is.numeric))), type = "upper",
           p.mat = cor_pmat(Pheno %>% select(where(is.numeric))),
           title = "Littorina saxatilis phenotypes")


#### C: Adjustments based on correlations ####

Pheno2 <- Pheno

### All measurements in mm are strongly correlated with shell length, with more variation
### from larger shells. To account for this heteroscedaticity I will scale dividing r0, z0, 
### a0 and shell thickness by shell length.
#Pheno2[, which(colnames(Pheno2) %in% c("r0", "z0", "a0", "avg_thick_mm"))] <- 
#  apply(Pheno2[, which(colnames(Pheno2) %in% c("r0", "z0", "a0", "avg_thick_mm"))], 
#        2, function(X){X / Pheno2$shellLength})

### Adjust eccentricity by aperture size (a0)
### This is a correction recommended by Jenny Larsson
### It scales the eccentricity measurement by the radius of the aperture
Pheno2$eccentricity <- Pheno2$eccentricity / Pheno2$a0

### Remove apAngle due to weak correlation with other variables
#Pheno2$apAngle <-  NULL

### Merge gw and gh to reduce total axes considered in PCA
### gh - gw = a measurement of the convexity of whorls along the shell
### This is a correction recommended by Jenny Larsson
#Pheno2$convexity <- Pheno2$gh- Pheno2$gw
#Pheno2$gh <- NULL
#Pheno2$gw <- NULL


#### D: PCA ####

### Scale and centre categorical data
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


#### E: Plot PCA ####

### E1: PC1 vs PC2 scatterplot
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


### E2: Screeplot 
# Copy ratio: 400 x 400
ggplot(Pvar, aes(x = as.numeric(gsub("PC", "", PC)), y = p.var))+
  geom_line()+
  geom_point()+
  labs(x = "Principal component", y = "Percentage varaince explained")+
  ylim(c(0, 100))+
  scale_x_continuous(breaks = 1:nrow(Pvar2))+
  theme_classic()

### E3: Histogram of PC1-5
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


### E4: PCA scatter plot PC1-9
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
ggarrange(ggarrange(PCA.Scatterplot(1,2), PCA.Scatterplot(1,3), PCA.Scatterplot(1,4),
                    PCA.Scatterplot(1,5), PCA.Scatterplot(1,6), PCA.Scatterplot(1,7),
                    PCA.Scatterplot(1,8), PCA.Scatterplot(1,9), nrow = 1),
          ggarrange(p_, PCA.Scatterplot(2,3), PCA.Scatterplot(2,4), 
                    PCA.Scatterplot(2,5), PCA.Scatterplot(2,6), PCA.Scatterplot(2,7),
                    PCA.Scatterplot(2,8), PCA.Scatterplot(2,9), nrow = 1),
          ggarrange(p_, p_, PCA.Scatterplot(3,4), 
                    PCA.Scatterplot(3,5), PCA.Scatterplot(3,6), PCA.Scatterplot(3,7),
                    PCA.Scatterplot(3,8), PCA.Scatterplot(3,9), nrow = 1),
          ggarrange(p_, p_, p_, 
                    PCA.Scatterplot(4,5), PCA.Scatterplot(4,6), PCA.Scatterplot(4,7),
                    PCA.Scatterplot(4,8), PCA.Scatterplot(4,9), nrow = 1),
          ggarrange(p_, p_, p_, 
                    p_, PCA.Scatterplot(5,6), PCA.Scatterplot(5,7),
                    PCA.Scatterplot(5,8), PCA.Scatterplot(5,9), nrow = 1),
          ggarrange(p_, p_, p_, p_, p_, PCA.Scatterplot(6,7), 
                    PCA.Scatterplot(6,8), PCA.Scatterplot(6,9), nrow = 1),
          ggarrange(p_, p_, p_, p_, p_, p_, 
                    PCA.Scatterplot(7,8), PCA.Scatterplot(7,9), nrow = 1),
          ggarrange(p_, p_, p_, p_, p_, p_, p_, PCA.Scatterplot(8,9), nrow = 1),
          ncol = 1)


### E5: Factor loading plots
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


#### F: Determine number of PCs to retain ####

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

# 75% = PC1-2 (almost)
# 80% = PC1-3
# 90% = PC1-4


#### G: Save PCA ####
OUT <- "/output/dir/"

# Save principal components
write.csv(PCs, paste0(OUT, "principal_components/Pheno_PC1-9.v2.csv"), quote = FALSE)

# Save eigenvalues and percent of variance explained
write.csv(Pvar, paste0(OUT, "eigenvalues/Pheno_variance_explained.csv"), quote = FALSE)
