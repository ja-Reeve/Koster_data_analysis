########## Isolation by distance analysis for Koster data ############
### Identifying if there is isolation by distance among sampling sites in the Koster dataset. IBD is not a focus 
### of this project, but something we want to understand, so we can account for it in later analysis. This script 
### was modified from Roger Butlin's original, which was based on the adegenet tutorials.

### James Reeve - University of Gothenburg
### 2024-02-02
### Edit 2025-09-16: drop SNPs in LGC5.1 from genetic analysis

### Preparation
rm(list = ls())
dev.off()
setwd("~")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(reshape2)

# Mapping
library(terra)
library(tidyterra)
library(raster)
library(sf)
library(rgdal)
library(geosphere)

# Genetic analysis
library(adegenet)
library(MASS)



#### A: Access data ####

# Filepath
PATH<- "/path/to/working/directory/"


### 1: SNP genotypes
SNP <- read.csv(paste0(PATH, "Koster_SNP_imputed_20221031.csv"))
# Subset to Koster data
SNP <- SNP[grep("KT", SNP$snail),]
# Remove duplicated snails
SNP <- SNP[!duplicated(SNP$snail),]


### 2: SNP information
SNP_info <- read.csv(paste0(PATH, "Koster_SNP_info.csv"))


### 3: Site details
Env <- read.csv(paste0(PATH, "KT_fieldnotes.v2.csv"))


### 4: Map data
tmp_shp <- vect("/path/to/map/data/coastline_clipped.shp")
# Subset to Kosterhavets
KT_shp <- terra::crop(tmp_shp, ext(10.95, 11.16, 58.82, 58.93))
rm(tmp_shp)



#### B: Data wrangling ####

### 1: Separate SNP and inversion genotypes
GT <- SNP[,c(1, grep("Contig", colnames(SNP)))]

# Further sample 'GT' to just genetic background
SNPs <- SNP_info %>% 
  filter(inv == "background") %>% # Only keep genetic background markers
  filter(!(LG == 5 & between(mp, 15, 47))) %>% # Drop SNPs in LGC5.1 (LG5: 15 - 47 cM)
  filter(!(str_detect(NAME, "Contig70288_47522"))) %>% # Finally, drop outlier which impacts PCA pattern
  reframe("SNP" = paste(CHROM, POS, sep = "_"))
GT <- GT[,SNPs]

# Add rownames
rownames(GT) <- SNP$snail

# Clean-up files
rm(SNPs, SNP_info)


### 2: Convert genotypes into adgenet format
GT2 <- GT
GT2[GT2 == 0] <- "00"
GT2[GT2 == 1] <- "01"
GT2[GT2 == 2] <- "11"
GT2[is.na(GT2)] <- "99"

# Reorder by snail name
GT2 <- GT2[order(rownames(GT2)),]


### 3: Convert into adgenet formats
# Genind (genotypes by individual)
gen <- df2genind(GT2, ploidy = 2, type = "codom", NA.char = "9", sep = "", pop = substr(rownames(GT2), 1, 5)) # Beware, removes snails where all values are missing
# Genpop (genotypes by population)
genp <- genind2genpop(gen)



#### C: PCA of genetic background ####

### Interpolate missing genotypes using mean
gen_scale <- scaleGen(gen, scale = FALSE, NA.method = "mean")

### PCA
pca1 <- dudi.pca(gen_scale, scale = FALSE, center = FALSE, scannf = FALSE, nf = 2)

# Find number of SNPs
Nsnp <- length(levels(gen$loc.fac))

# Variance explained (%)
Vexp.PC1 <- round(pca1$eig[1] / sum(pca1$eig) * 100, digits = 2)
Vexp.PC2 <- round(pca1$eig[2] / sum(pca1$eig) * 100, digits = 2)

# Plot
tmp <- as.data.frame(pca1$li)
tmp$Snail <- rownames(tmp)
tmp$Site <- substr(tmp$Snail, 1, 5)

ggplot(tmp)+
  geom_point(aes(x = Axis1, y = Axis2), alpha = 0.4)+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"),
       title = paste(Nsnp, "SNPs"))+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "grey50", linewidth = 0.1),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, angle = 90),
        plot.title = element_text(size = 18, face = "bold"))
  
# 600 x 450



#### D: Mantel test for Isolation by Distance ####
### Mantel tests look for significant differences among two correlation matrices.
### Comparing a spatial distance and genetic distance matrix can be used to 
### test for the presence of isolation by distance.

### 1: Genetic distance matrix
Dgen <- dist.genpop(genp, method = 2, diag = TRUE, upper = TRUE) # Method 2 = Edward's distance


### 2: Spatial distance matrix
Env <- Env[order(Env$Site),]
# Only keep sites with snails
Env2 <- Env[Env$Site %in% rownames(genp$tab),]

# Correction: replace bottom coordinates with top at KT148
Env2[Env2$Site == "KT148", "Long_bottom"] <- Env2[Env2$Site == "KT148", "Long_top"]
Env2[Env2$Site == "KT148", "Lat_bottom"] <- Env2[Env2$Site == "KT148", "Lat_top"]

# Spatial distances (km)
Dspat <- spDists(x = as.matrix(Env2[,c("Long_bottom", "Lat_bottom")]), longlat = TRUE)
colnames(Dspat) <- Env2$Site
rownames(Dspat) <- Env2$Site


### 3: Mantel test
mantel.randtest(as.dist(Dspat), as.dist(Dgen))


### 4: Convert matrices to a data frame
## Custom function
mat.long <- function(x){
  tmp <- melt(x)
  colnames(tmp) <- c("site.x", "site.y", "r")
  # Drop autocorrelations
  tmp <- tmp[tmp$site.x != tmp$site.y,]
  
  return(tmp)
}

refmt.Dspat <- mat.long(Dspat)
refmt.Dgen <- mat.long(Dgen)

### Merge together
# Add a comparison flag between sites
refmt.Dspat$comp <- paste(refmt.Dspat$site.x, refmt.Dspat$site.y, sep = "_")
refmt.Dgen$comp <- paste(refmt.Dgen$site.x, refmt.Dgen$site.y, sep = "_")

# The merger
tmp <- merge(refmt.Dspat, refmt.Dgen[,c("r", "comp")], by = "comp")
colnames(tmp)[4:5] <- c("SD", "GD")


### 5: Plot comparison among matrices

mod <- lm(tmp$Dgen ~ log(tmp$Dspat))
# Distances > 2km
tmp.2km <- tmp[tmp$Dspat > 2,]
mod.2km <- lm(tmp.2km$GD ~ log(tmp.2km$SD))

ggplot(tmp, aes(x = SD, y = GD))+ 
  geom_density_2d_filled(n = 300, alpha = 1.0, show.legend = FALSE, bins = 9)+ # Runs MASS:kde2d in ggplot
  scale_fill_brewer(palette = "Greens")+
  geom_point(size = 0.2, alpha = 0.1)+
  geom_smooth(method = "lm", formula = y~log(x), colour = "red", linewidth = 0.4)+
  labs(x = "Spatial distances (km)", y = "Genetic distances (E)")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axix.text = element_text(size = 14),
        axis.title = element_text(size = 18))
# 900 x 600



#### E: Spatial PCA ####

### 1: Site coordiantes
Coords <- Env2[,c("Long_bottom","Lat_bottom")]


### 2: Run spatial PCA
# Prompts: connection network = type 1 (Delaunay triangulation)
#         Number of axes with positive spatial autocorrelation = 1
#         Number of axes with negative spatial autocorrelation = 0
Spca <- spca(genp, Coords)


### 3: Map
Coords$Similarity <- Spca$li$`Axis 1`
  
ggplot()+
  geom_spatvector(data = KT_shp, colour = "grey70")+
  geom_point(data = Coords, aes(Long_bottom, Lat_bottom, fill = Similarity), size = 2, pch = 21)+
  scale_fill_gradient2(low = "navy", mid = "#ebffff", high = "firebrick", midpoint = 0)+
  labs(fill = "sPC1")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        legend.position = "bottom",
        legend.key.width = unit(1.6, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# 450 x 450
