######### Snail-scape: Choosing sampling sites ##########
### This script shows how I picked sampling sites for a PhD project
### exploring the genetic diversity present in sea snail (Littorina saxatilis)
### from the region around Tjärnö, Sweden.
### James Reeve 
### Göteborgs Universitet
### 24/03/2020

#### Preparation ####
rm(list = ls())
dev.off()
setwd("..") # PC
#setwd("~") #Mac

### Packages
library(SpatialTools)
library(png)
library(grid)
library(ggplot2)

### Data
# Image of field site - screenshot from Google Maps
img <- png::readPNG("input/dir/Tjarno_Koster_google_map.png")

#### A. Create grid (KTR = [K]osterhavets - [T]järnö [R]egion) ####

### Define parameters
nx <- 100 # Number of rows
ny <- 100 # Number of columns

### Create a coordinates matrix
col.num <- rep(1:nx, each = ny)
row.num <- rep(1:ny, nx)
KTR <- cbind(row.num, col.num)

### Convert img file into a plot object (i.e. grob)
KTR_map <- rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

### reorder matrix to dataframe
KTR.df <- as.data.frame(KTR)
# add a cell number vector
KTR.df$cell <- 1:10000

### Plot region map with overlaid grid
ggplot(KTR.df, aes(x = long, y = lat))+
  annotation_custom(KTR_map,  -Inf, Inf, -Inf, Inf)+ # Adds image to background!
  geom_tile(colour = "lightgrey", alpha = 0.0)+
  geom_text(aes(label = cell), size = 0.8, colour = "white")+
  scale_y_continuous(expand = c(0,0))+ # Removes margin between plot and axes
  scale_x_continuous(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none", 
        axis.title = element_blank())

#### B. Create a probability matrix for sampling ####
### The chance a site is pick is determined by the probability matrix

### Define centre of probability matrix
# Ängklovebukten (ANG)
ANG <- rbind(c(x = 76, y = 58)) 

### Calcualte distances from each gird location to ANG - using SpatialData
# This approach uses Euclidian distances among cells
distances <- dist2(KTR,ANG)

### Calculate probabilities
# create function that transforms distances into proportions
scale_prop <- function(x, m = 0) {
  (x - min(x)) / (m + max(x) - min(x))}

# Take the inverse proportion as a probability
## linear ### probs <- 1 - scale_prop(distances, m = 1) 
## log-linear
probs <- 0.9*exp(-5*(distances/max(distances)))+0.1

### Test it worked by visualising distances
image(matrix(distances, nrow = nx, byrow = TRUE))

### Test 2; visualising probabilities
# It should be inverse of previous (large dist = low prob)
image(matrix(probs, nrow = nx, byrow = TRUE))

#### C. exclude unsuitable sites ####
# Load vector of cells to keep
# These sites were determined manually by removing any cells which did not have a coastline (i.e., on land or open sea).
keeps <- read.table("inpiut/dir/keeps.txt")
tmp <- strsplit(as.character(keeps[1,1]), split = ",")
tmp2 <- sapply(unlist(tmp), function(X){eval(parse(text = X))})
keeps <- unlist(tmp2)
rm(tmp, tmp2)

### Subset probabilites vector
keeps.probs <- probs[keeps]

### Create subset of 200 sites
Sites <- sample(keeps, size = 200, prob = keeps.probs)

### Visualize on map
# Add new columns to data frame
KTR.df$site <- KTR.df$cell %in% Sites
KTR.df$coastal <- KTR.df$cell %in% keeps

ggplot(KTR.df, aes(x = long, y = lat))+
  annotation_custom(KTR_map,  -Inf, Inf, -Inf, Inf)+ ### Adds image to background!
  geom_tile(aes(fill = site, alpha = coastal), colour = "lightgrey")+
  scale_alpha_manual(values = c(0.5,0.3))+
  geom_text(aes(label = cell), colour = "white", size = 0.8)+
  scale_fill_manual(values = c("black", "yellow"))+
  scale_y_continuous(expand = c(0,0))+ 
  scale_x_continuous(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none", 
        axis.title = element_blank())

### Add coordinates
N <- 58.908934 # Northern boundary
S <- 58.820627 # Southern boundary
E <- 11.159007 # Eastern boundary
W <- 10.988648 # Western boundary

Dist_lat <- (N - S) / 100
Dist_long <- (E - W) / 100

lat_calculator <- function(cell.number){S + Dist_lat * cell.number - (Dist_lat / 2)}
long_calculator <- function(cell.number){W + Dist_long * cell.number - (Dist_long / 2)}

# Adding coordinates to KTR.df usign the row and column numbers
KTR.df$lat <- lat_calculator(KTR.df$row.num)
KTR.df$long <- long_calculator(KTR.df$col.num)

#### D. Save the data ####
write.csv(KTR.df[KTR.df$site == TRUE,], 
          "out/dir/sampling_sites.csv", quote = FALSE)
