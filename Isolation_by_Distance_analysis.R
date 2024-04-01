########## Isolation by distance analysis for Koster data ############
### Identifying if there is isolation by distance among sampling sites in the Koster dataset. IBD is not a focus 
### of this project, but something we want to understand, so we can account for it in later analysis. This script 
### was modified from Roger Butlin's original, which was based on the adegenet tutorials.

### James Reeve - University of Gothenburg
### 2024-02-02

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
tmp_shp <- vect("/path/to/map/data/coastline_line.shp")
# Subset to Kosterhavets
KT_shp <- terra::crop(tmp_shp, ext(10.95, 11.16, 58.82, 58.93))
rm(tmp_shp)



#### B: Data wrangling ####

### 1: Separate SNP and inversion genotypes
GT <- SNP[,c(1, grep("Contig", colnames(SNP)))]

# Further sample 'GT' to just genetic background
SNPs <- SNP_info %>% 
  filter(inv == "background") %>%
  summarise("SNP" = paste(CHROM, POS, sep = "_")) %>% 
  filter(SNP != "Contig70288_47522") %>% # Removing SNP which contributed strongly to PCA
  unlist()
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
GT2 <- GT2[order(GT2$snail),]

# Drop 'snail' column
GT2 <- GT2[,-1]


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
  annotate("text", x = 4.5, y = -3, size = 6, label = paste(Nsnp, "SNPs"))+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"), y = paste0("PC2 (", Vexp.PC2, "%)") )+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))
# 600 x 450


#### D: Mantel test for Isolation by Distance ####
### Mantel tests look for significant differences among two correlation matrices.
### Comparing a spatial distance and genetic distance matrix can be used to 
### test for the presence of isolation by distance.

### 1: Genetic distance matrix
Dgen <- dist.genpop(genp, method = 2) # Method 2 = Edward's distance


### 2: Spatial distance matrix
Env <- Env[order(Env$Site),]
# Only keep sites with snails
Env2 <- Env[Env$Site %in% rownames(genp$tab),]

# Correction: replace bottom coordinates with top at KT148
Env2[Env2$Site == "KT148", "Long_bottom"] <- Env2[Env2$Site == "KT148", "Long_top"]
Env2[Env2$Site == "KT148", "Lat_bottom"] <- Env2[Env2$Site == "KT148", "Lat_top"]

# Spatial distances (km)
Dspat <- distm(Env2[,c("Long_bottom","Lat_bottom")])/1000


### 3: Mantel test
mantel.randtest(as.dist(Dspat), Dgen)


### 4: Convert matrices to a data frame

# Turn genetic distance matrix to a data frame
tmp <- Dgen %>% as.matrix() %>% melt()
colnames(tmp) <- c("Site.x", "Site.y", "Dgen")

# Convert spatial distances to a single vector
Dspat.l <- melt(Dspat)$value
# Add spatial distances to data frame
tmp$Dspat <- Dspat.l

# Remove autocorrelations
tmp <- tmp[tmp$Site.x != tmp$Site.y,]


### 4: Plot comparison among matrices
# Convert into densities using MASS
#dens <- kde2d(tmp$Dspat, tmp$Dgen, n = 300)

mod <- lm(tmp$Dgen ~ log(tmp$Dspat))
# Distances > 2km
tmp.2km <- tmp[tmp$Dspat > 2,]
mod.2km <- lm(tmp.2km$Dgen ~ log(tmp.2km$Dspat))

ggplot(tmp, aes(x = Dspat, y = Dgen))+ 
  geom_point(size = 0.2, alpha = 0.4)+
  geom_density_2d_filled(n = 300, alpha = 0.5, show.legend = FALSE)+ # Runs MASS:kde2d in ggplot
  scale_fill_brewer(palette = "Greens")+
  geom_smooth(method = "lm", formula = y~log(x), colour = "red", linewidth = 0.4)+
  labs(x = "Spatial distances (km)", y = "Genetic distances (E)")+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
# 800 x 400



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
  geom_spatvector(data = KT_shp)+
  geom_point(data = Coords, aes(Long_bottom, Lat_bottom, fill = Similarity), size = 2, pch = 21)+
  scale_fill_gradient(low = "#eff7f7ff", high = "#008080ff")+
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1.6, "cm"))
# 450 x 450
