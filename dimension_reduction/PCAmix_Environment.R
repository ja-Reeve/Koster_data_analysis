############# Perform PCA on Environmental data ######################
### This script performs a simple principle component analysis on
### environmental measurements from 200 Littorina saxatilis sampling 
### sites around Tjärnö, Sweden that were collected in the summer of 2020.
### These measurement come from drone photos taking over each site
### James Reeve - University of Gothenburg
### 12/04/2021

### Measurements ####
### Fieldnotes
#   i. "site" = name of sampling site
#   ii. "Site_name" = name of location or island
#   iii. "Lat_bottom" = latitude at the bottom of the transect [ignored]
#   iv. "Long_bottom" = longitude at the bottom of the transect, and a character in the Harry Potter series [ignored]
#   v. "Lat_top" = latitude at the top of the transect [ignored]
#   vi. "Long_top" = longitude at the top of the transect [ignored]
#   vii. "Date" = date of collection (DD/MM/YYYY)
#   viii. "Time" = time of collection (Swedish time - GMT+1 & GMT+2 in Summer) [ignored]
#   ix. "Height_bottom" = height from sea level at the bottom of the transect (m) [ignored]
#   x. "Height_top" = Height from sea level at the top of the transect (m) [ignored]
#   xi. "Sea_level" = SMHI sea level record from bouy at Kungsvik (cm) [https://www.smhi.se/data/hydrologi/vattenstand - in Swedish] [ignored]
#   xii. "Adj_height_bottom" = Height_bottom + Sea_level/100 (m) [ignored]
#   xiii. "Adj_height_top" = Height_top + Sea_level/100 (m) [ignored]
#   xiv. "Transect_height" = Adj_height_top - Adj_height_bottom (m)
#   xV. "slope" = acos(shore_length/Transect_height) (radians)
#   xvi. "Habitat" = substrate category of each site; R=rock, B=boulder, M=mud, S=sand & H..=hybrid zone
#   xvii. "Exposure_notes" = exposure to waves for site; S=sheltered, E=exposed, X.=exterme site
#   xviii. "Fucus" = presence of Fucus vesiculousus, 'PATCH' sites had only small patches
#   xix. "Ascophyllum" = presence of Ascphyllum nodosum, 'PATCH' sites only had small patches
#   xx. "Barnacles" = presence of barnacle bed, 'PATCH' sites only had small patches
#   xxi. "Date_redo_photo" = date of second flight of drone if needed [ignored]
#   xxii. "Time_redo_photo" = time of second flight of drone if needed [ignored]
#   xxiii. "TimeStamp" = 'Date' + 'Time', used for adjusting heights - see "Add_sea_level_to"fieldnotes.R" [ignored]
#   xxiv. "Time.UTC" = time of collection in UTC

### Drone photo - width measurements are mean of 5 lines along shore
#   i. "Site" = name of sampling site
#   ii. "Habitat" = rock | boulder | mud-sand
#   iii. "Boulder_rock_stones" = presence of patch of alternative habitat (TRUE|FALSE)
#   iv. "HI" = presence of human infrastructure (TRUE|FALSE)
#   v. "Barn" = presence of barnacles (TRUE|FALSE)
#   vi. "FW" = Fucus belt width (m)
#   vii. "FP" = Fucus patchiness (%)
#   viii. "BZW" = black zone width (m)
#   ix. "TZW" = tidal zone width (m)
#   x. "GZW" = grey zone width (m)
#   xi. "LVD" = land vegetation distance (m)
#   xii. "HTD" = distance to habitat transition (m)
#   xiii. "GS" = grain size of boulders (integer)
#   xiv. "RF" = rock fissure count (integer)
# For a more detailed description of each of the measurements please see the supplementry
# method: https://docs.google.com/document/d/1cCgqbYA1x18d0wYcO4LOfPqBRq25GXWH/edit?usp=sharing&ouid=108508049506954195112&rtpof=true&sd=true

### GIS data - from Per Bergström's fetch model
#   i. "Site" = name of sampling site
#   ii. "Depth" = Depth drop-off near site - not measured [ignored]
#   iii. "Exposure_FETCH_m.2.s" = wave exposure estimate from fetch model (m^2/s)
#   iv. "Aspect_shore_ºN" = orientation of the shore in decimal degrees clockwise from North
#   v. "Reach_shore" = distance to next shore line - not measured [ignored]

### Edits
#   14/03/2023: drop GZW and transect height due to correlation with other variables
#               and remove exposure notes, due to strong influence on PCAmix
#   17/03/2023: select PCs to keep using the package PCAtest()
#   24/09/2025: Clean up code and added broken stick model to scree plot

### Preparation ####
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)
if(Sys.info()['sysname'] == "Windows") {
  setwd("C:/Users/jalre")
} else {
  setwd("/Users/james/")
}

### Pacakges
library(tidyverse)
library(ggfortify)
library(ggpubr)
library(ggcorrplot)
library(PCAmixdata)
#devtools::install_github("arleyc/PCAtest")
library(PCAtest)


#### A: Access and wrangle data ####

PATH <- "/file/path/"

### Access data
# Field notes for each site - after adding sea level data
FB_Site <- read.csv(paste0(PATH, "Environmental_data/KT_fieldnotes.v2.csv"))
# Drone photo measurements
DPh <- read.csv(paste0(PATH, "Environmental_data/KT_drone_photos_v2.csv"))
# GIS data - Per Bergström
GisD <- read.csv(paste0(PATH, "Environmental_data/KT_GIS_data_v2.csv"))

### Remove unused variables
FB_Site <- FB_Site[,c(1,7,14:20,24)]
DPh <- DPh[,c(1:ncol(DPh)-1)]

### Merge date
Env <- merge(GisD, merge(FB_Site, DPh, by = "Site"), by = "Site")

### Replace '?' with 'Un'
Env[Env$Ascophyllum %in% c("???", "????", "?????"), "Ascophyllum"] <- "Un"
Env[Env$Barnacles == "????", "Barnacles"] <- "Un"

### Calculate shore slope
# shore height = field notes
# shore length = drone photo measurement
# slope = arc tanget (shore height / shore length) [radians]
Env$slope <- atan(Env$Transect_height / Env$Shore_length)


### Adjust drone photo measurements by shore slope
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

### Adjust aspect to make 0 and 359 close
# Roger suggested to adjust aspect as the sin and cos of aspect (02/03/2023)
Env$cos_aspect <- cos(pi*Env$Aspect_shore/180)
Env$sin_aspect <- sin(pi*Env$Aspect_shore/180)

### Infer missing field note measurements of Ascophyllum and Barnacles using drone photos
Env$Ascophyllum <- ifelse(Env$Ascophyllum == "Un", Env$Asc, Env$Ascophyllum)
Env$Barnacles <- ifelse(Env$Barnacles == "Un", Env$Barn, Env$Barnacles)

### Select variables to keep
keeps <- c("Site", "Habitat.y", "Exposure_FETCH", "sin_aspect", "cos_aspect", 
           "slope", "Shore_length", "Fucus", "Ascophyllum", "Barnacles", "Boulders_rocks_stones",
           "HI", "Adj.FW", "FP", "Adj.BZW", "Adj.TZW", "Adj.LVD", "HTD", "GS", "RF")
Env <- Env[, keeps]

### Fix name of "Habitat"
colnames(Env)[2] <- "Habitat"

### Split data into three different habitat categories
### Doing this to avoid NA's in correlations
RS <- Env[Env$Habitat == "rock", !colnames(Env) %in% c("Habitat", "GS", "Adj.TZW")]
BS <- Env[Env$Habitat == "boulder", !colnames(Env) %in% c("Habitat", "RF", "Adj.TZW")]
MSS <- Env[Env$Habitat == "mud-sand", !colnames(Env) %in% c("Habitat", "GS", "RF", "Adj.BZW")]

# Drop sites without snails
RS <- RS[complete.cases(RS),]
BS <- BS[complete.cases(BS),]
MSS <- MSS[complete.cases(MSS),]



#### B: Investigate data ####

### B1: Investigate categorical data

# Identify categorical variables
CatVar <- which(sapply(Env, is.character) | sapply(Env, is.logical))

# View counts of each variable
sapply(Env[,CatVar], table)

# Plot each variable
cat.plot <- function(X, title){
  ggplot(Env)+
    geom_bar(aes(X, fill = Habitat), 
             position = position_dodge2(width = 0.9, preserve = "single"))+
    labs(x = title)+
    scale_fill_manual(values = c("#a2c4c9", "#b4a7d6", "#ea9999"))+
    theme_bw()+
    theme(axis.title.y = element_blank())
}

p1 <- ggplot(Env)+
  geom_bar(aes(Habitat, fill = Habitat))+
  labs(x = "Habitats")+
  scale_fill_manual(values = c("#a2c4c9", "#b4a7d6", "#ea9999"))+
  theme_bw()+
  theme(axis.title.y = element_blank())

p2 <- cat.plot(Env$Fucus, title = "Fucus present")
p3 <- cat.plot(Env$Ascophyllum, title = "Ascophyllum present")
p4 <- cat.plot(Env$Barnacles, title = "Barnacles present") 
p5 <- cat.plot(Env$Boulders_rocks_stones, title = "Microhabitat")
p6 <- cat.plot(Env$HI, title = "Human infrastructure")

ggarrange(p1,p2,p3,p4,p5,p6, nrow = 2, ncol = 3, 
          common.legend = TRUE, legend = "bottom")


###  Alternative plot of categorical variables
tmp <- Env %>% select(Habitat, Ascophyllum,
                      Barnacles, HI, Boulders_rocks_stones) %>%
  gather("variable", "score", -Habitat) %>%
  group_by(Habitat, variable) %>% 
  count(score)

ggplot(tmp)+
  geom_bar(aes(x = variable, y = n, alpha = score, 
               fill = Habitat, col = Habitat),
           stat = "identity", position = "dodge", show.legend = FALSE)+
  facet_wrap(vars(Habitat), strip.position = "top", nrow = 3, ncol = 1)+
  scale_colour_manual(values = c("#9dbfc4ff", "#ad9ed1ff", "#e69193ff"))+
  scale_fill_manual(values = c("#9dbfc4ff", "#ad9ed1ff", "#e69193ff"))+
  scale_alpha_manual(values = c(0.2, 0.6, 1))+
  scale_x_discrete(labels = c("Ascophyllum", "Barnacles", 
                              "Human infrastructure", "Microhabitat"))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


### B2: Investigate continuous variables

# Identify continuous variables
ConVar <- which(sapply(Env, is.numeric) | sapply(Env, is.integer))

# View summary stats for each variable
sapply(Env[,ConVar], summary)

# Plot each variable
cont.plot <- function(X, title){
  tmp <- data.frame("Habitat" = Env$Habitat, 
                    "Var" = X)
  tmp2 <- tmp %>% group_by(Habitat) %>% 
                summarise("Avg" = mean(Var, na.rm = TRUE),
                          "SD" = sd(Var, na.rm = TRUE))
  ggplot(Env, aes(x = Habitat))+
    geom_violin(aes(y = X, fill = Habitat), alpha = 0.4)+
    geom_jitter(aes(y = X), height = 0, width = 0.1, size = 0.3)+
    geom_point(data = tmp2, aes(y = Avg, x = Habitat), 
               colour = "grey50", size = 2, pch = 15)+
    geom_errorbar(data = tmp2, aes(ymin = Avg - SD, ymax = Avg + SD), 
                  width = 0.3, colour = "grey50")+
    scale_fill_manual(values = c("#a2c4c9", "#b4a7d6", "#ea9999"))+
    labs(x = "Habitat", y = title)+
    theme_bw()+
    theme(axis.title = element_text(size = 8))
}

p1 <- cont.plot(log10(Env$Exposure_FETCH), title = expression(Log[10](wave~fetch)))
p2 <- cont.plot(Env$sin_aspect, title = "Sin(Shore aspect)")
p3 <- cont.plot(Env$slope, title = "Shore angle (radians)")
p4 <- cont.plot(Env$FP, title = "Fucus patchiness")
p5 <- cont.plot(Env$HTD, title = "Habitat transition distance (m)")
p6 <- cont.plot(Env$GS, title = "Boulder numbers")
p7 <- cont.plot(Env$RF, title = "Number of rock fisures") + scale_fill_manual(values = c("#ea9999"))
p8 <- cont.plot(Env$Adj.FW, title = "FW (m)")
p9 <- cont.plot(Env$Adj.BZW, title = "BZW (m)") + scale_fill_manual(values = c("#a2c4c9", "#ea9999"))
p10 <- cont.plot(Env$Adj.TZW, title = "TZW (m)") + scale_fill_manual(values = c("#b4a7d6"))
p11 <- cont.plot(Env$Adj.LVD, title = "LVD (m)")

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, ncol = 6, nrow = 2,
          common.legend = TRUE, legend = "bottom")


#### B3: Correlations ####
### This section will generate multiple plots to see the two-way correlations between environmental variables.
### I advise running it line by line.

### 2-way correlation plots plots
scat.plots <- function(x, y, data, title.x, title.y){
  ggplot(data, aes(x, y))+
    geom_point()+
    geom_smooth()+
    labs(x = title.x, y = title.y, title = paste("rho =", round(cor(x,y), 2)))+
    theme_bw()
}

### Rocky shore correlations ####
scat.plots(log10(RS$Exposure_FETCH), RS$sin_aspect, RS, 
           expression(Log[10](fetch)), "Shore aspect")
scat.plots(log10(RS$Exposure_FETCH), RS$Shore_length, RS, 
           expression(Log[10](fetch)), "Shore length (m)")
scat.plots(log10(RS$Exposure_FETCH), RS$slope, RS, 
           expression(Log[10](fetch)), "Transect slope (radians)")
scat.plots(log10(RS$Exposure_FETCH), RS$HTD, RS, 
           expression(Log[10](fetch)), "Habitat transition distance (m)")
scat.plots(log10(RS$Exposure_FETCH), RS$FP, RS, 
           expression(Log[10](fetch)), "Fucus patchiness (%)")
scat.plots(log10(RS$Exposure_FETCH), RS$RF, RS, 
           expression(Log[10](fetch)), "Number of rock fissures")
scat.plots(log10(RS$Exposure_FETCH), RS$Adj.FW, RS, 
           expression(Log[10](fetch)), "Fucus belt width (m)")
scat.plots(log10(RS$Exposure_FETCH), RS$Adj.BZW, RS, 
           expression(Log[10](fetch)), "Black zone width (m)")
scat.plots(log10(RS$Exposure_FETCH), RS$Adj.LVD, RS, 
           expression(Log[10](fetch)), "Land vegetation distance (m)")

scat.plots(RS$sin_aspect, RS$Shore_length, RS, 
           "Aspect shore", "Shore length (m)")
scat.plots(RS$sin_aspect, RS$slope, RS, 
           "Aspect shore", "Transect slope (radians)")
scat.plots(RS$sin_aspect, RS$HTD, RS, 
           "Aspect shore", "Habitat transition distance (m)")
scat.plots(RS$sin_aspect, RS$FP, RS, 
           "Aspect shore", "Fucus patchiness (%)")
scat.plots(RS$sin_aspect, RS$RF, RS, 
           "Aspect shore", "Number of rock fissures")
scat.plots(RS$sin_aspect, RS$Adj.FW, RS, 
           "Aspect shore", "Fucus width (m)")
scat.plots(RS$sin_aspect, RS$Adj.BZW, RS, 
           "Aspect shore", "Black zone width (m)")
scat.plots(RS$sin_aspect, RS$Adj.LVD, RS, 
           "Aspect shore", "Land vegetation distance (m)")

scat.plots(RS$Shore_length, RS$slope, RS, 
           "Shore length (m)", "Transect slope (radians)")
scat.plots(RS$Shore_length, RS$HTD, RS, 
           "Shore length (m)", "Habitat transition distance (m)")
scat.plots(RS$Shore_length, RS$FP, RS, 
           "Shore length (m)", "Fucus patchiness (%)")
scat.plots(RS$Shore_length, RS$RF, RS, 
           "Shore length (m)", "Number of rock fissures")
scat.plots(RS$Shore_length, RS$Adj.FW, RS, 
           "Shore length (m)", "Fucus width (m)")
scat.plots(RS$Shore_length, RS$Adj.BZW, RS, 
           "Shore length (m)", "Black zone width (m)")
scat.plots(RS$Shore_length, RS$Adj.LVD, RS, 
           "Shore length (m)", "Land vegetation distance (m)")

scat.plots(RS$slope, RS$HTD, RS, 
           "Transect slope (radians)", "Habitat transition distance (m)")
scat.plots(RS$slope, RS$FP, RS, 
           "Transect slope (radians)", "Fucus patchiness (%)")
scat.plots(RS$slope, RS$RF, RS, 
           "Transect slope (radians)", "Number of rock fissures")
scat.plots(RS$slope, RS$Adj.FW, RS, 
           "Transect slope (radians)", "Fucus width (m)")
scat.plots(RS$slope, RS$Adj.BZW, RS, 
           "Transect slope (radians)", "Black zone width (m)")
scat.plots(RS$slope, RS$Adj.LVD, RS, 
           "Transect slope (radians)", "Land vegetation distance (m)")

scat.plots(RS$HTD, RS$FP, RS, 
           "Habitat transition distance (m)", "Fucus patchiness (%)")
scat.plots(RS$HTD, RS$RF, RS, 
           "Habitat transition distance (m)", "Number of rock fissures")
scat.plots(RS$HTD, RS$Adj.FW, RS, 
           "Habitat transition distance (m)", "Fucus width (m)")
scat.plots(RS$HTD, RS$Adj.BZW, RS, 
           "Habitat transition distance (m)", "Black zone width (m)")
scat.plots(RS$HTD, RS$Adj.LVD, RS, 
           "Habitat transition distance (m)", "Land vegetation distance (m)")

scat.plots(RS$FP,RS$RF, RS, 
           "Fucus patchiness (%)", "Number of rock fissures")
scat.plots(RS$FP, RS$Adj.FW, RS, 
           "Fucus patchiness (%)", "Fucus width (m)")
scat.plots(RS$FP, RS$Adj.BZW, RS, 
           "Fucus patchiness (%)", "Black zone width (m)")
scat.plots(RS$FP, RS$Adj.LVD, RS, 
           "Fucus patchiness (%)", "Land vegetation distance (m)")

scat.plots(RS$RF, RS$Adj.FW, RS, 
           "Number of rock fissures", "Fucus width (m)")
scat.plots(RS$RF, RS$Adj.BZW, RS, 
           "Number of rock fissures", "Black zone width (m)")
scat.plots(RS$RF, RS$Adj.LVD, RS, 
           "Number of rock fissures", "Land vegetation distance (m)")

scat.plots(RS$Adj.FW, RS$Adj.BZW, RS, 
           "Fucus width (m)", "Black zone width (m)")
scat.plots(RS$Adj.FW, RS$Adj.LVD, RS, 
           "Fucus width (m)", "Land vegetation distance (m)")

scat.plots(RS$Adj.FW, RS$Adj.LVD, RS, 
           "Black zone width (m)", "Land vegetation distance (m)")


rscor <- ggcorrplot(cor(RS %>% select(where(is.numeric))), type = "upper",
           p.mat = cor_pmat(RS %>% select(where(is.numeric))),
           title = "Rocky shores")


### Boulder shore correlations ####
scat.plots(log10(BS$Exposure_FETCH), BS$sin_aspect, BS, 
           expression(Log[10](fetch)), "Shore aspect")
scat.plots(log10(BS$Exposure_FETCH), BS$Shore_length, BS, 
           expression(Log[10](fetch)), "Shore length (m)")
scat.plots(log10(BS$Exposure_FETCH), BS$slope, BS, 
           expression(Log[10](fetch)), "Transect slope (radians)")
scat.plots(log10(BS$Exposure_FETCH), BS$HTD, BS, 
           expression(Log[10](fetch)), "Habitat transition distance (m)")
scat.plots(log10(BS$Exposure_FETCH), BS$FP, BS, 
           expression(Log[10](fetch)), "Fucus patchiness (%)")
scat.plots(log10(BS$Exposure_FETCH), BS$GS, BS, 
           expression(Log[10](fetch)), "Number of boulders")
scat.plots(log10(BS$Exposure_FETCH), BS$Adj.FW, BS, 
           expression(Log[10](fetch)), "Fucus belt width (m)")
scat.plots(log10(BS$Exposure_FETCH), BS$Adj.BZW, BS, 
           expression(Log[10](fetch)), "Black zone width (m)")
scat.plots(log10(BS$Exposure_FETCH), BS$Adj.LVD, BS, 
           expression(Log[10](fetch)), "Land vegetation distance (m)")

scat.plots(BS$sin_aspect, BS$Shore_length, BS, 
           "Aspect shore", "Shore length (m)")
scat.plots(BS$sin_aspect, BS$slope, BS, 
           "Aspect shore", "Transect slope (radians)")
scat.plots(BS$sin_aspect, BS$HTD, BS, 
           "Aspect shore", "Habitat transition distance (m)")
scat.plots(BS$sin_aspect, BS$FP, BS, 
           "Aspect shore", "Fucus patchiness (%)")
scat.plots(BS$sin_aspect, BS$GS, BS, 
           "Aspect shore", "Number of boulders")
scat.plots(BS$sin_aspect, BS$Adj.FW, BS, 
           "Aspect shore", "Fucus width (m)")
scat.plots(BS$sin_aspect, BS$Adj.BZW, BS, 
           "Aspect shore", "Black zone width (m)")
scat.plots(BS$sin_aspect, BS$Adj.LVD, BS, 
           "Aspect shore", "Land vegetation distance (m)")

scat.plots(BS$Shore_length, BS$slope, BS, 
           "Shore length (m)", "Transect slope (radians)")
scat.plots(BS$Shore_length, BS$HTD, BS, 
           "Shore length (m)", "Habitat transition distance (m)")
scat.plots(BS$Shore_length, BS$FP, BS, 
           "Shore length (m)", "Fucus patchiness (%)")
scat.plots(BS$Shore_length, BS$GS, BS, 
           "Shore length (m)", "Number of boulders")
scat.plots(BS$Shore_length, BS$Adj.FW, BS, 
           "Shore length (m)", "Fucus width (m)")
scat.plots(BS$Shore_length, BS$Adj.BZW, BS, 
           "Shore length (m)", "Black zone width (m)")
scat.plots(BS$Shore_length, BS$Adj.LVD, BS, 
           "Shore length (m)", "Land vegetation distance (m)")

scat.plots(BS$slope, BS$HTD, BS, 
           "Transect slope (radians)", "Habitat transition distance (m)")
scat.plots(BS$slope, BS$FP, BS, 
           "Transect slope (radians)", "Fucus patchiness (%)")
scat.plots(BS$slope, BS$GS, BS, 
           "Transect slope (radians)", "Number of boulders")
scat.plots(BS$slope, BS$Adj.FW, BS, 
           "Transect slope (radians)", "Fucus width (m)")
scat.plots(BS$slope, BS$Adj.BZW, BS, 
           "Transect slope (radians)", "Black zone width (m)")
scat.plots(BS$slope, BS$Adj.LVD, BS, 
           "Transect slope (radians)", "Land vegetation distance (m)")

scat.plots(BS$HTD, BS$FP, BS, 
           "Habitat transition distance (m)", "Fucus patchiness (%)")
scat.plots(BS$HTD, BS$GS, BS, 
           "Habitat transition distance (m)", "Number of boulders")
scat.plots(BS$HTD, BS$Adj.FW, BS, 
           "Habitat transition distance (m)", "Fucus width (m)")
scat.plots(BS$HTD, BS$Adj.BZW, BS, 
           "Habitat transition distance (m)", "Black zone width (m)")
scat.plots(BS$HTD, BS$Adj.LVD, BS, 
           "Habitat transition distance (m)", "Land vegetation distance (m)")

scat.plots(BS$FP, BS$GS, BS, 
           "Fucus patchiness (%)", "Number of boulders")
scat.plots(BS$FP, BS$Adj.FW, BS, 
           "Fucus patchiness (%)", "Fucus width (m)")
scat.plots(BS$FP, BS$Adj.BZW, BS, 
           "Fucus patchiness (%)", "Black zone width (m)")
scat.plots(BS$FP, BS$Adj.LVD, BS, 
           "Fucus patchiness (%)", "Land vegetation distance (m)")

scat.plots(BS$GS, BS$Adj.FW, BS, 
           "Number of boulders", "Fucus width (m)")
scat.plots(BS$GS, BS$Adj.BZW, BS, 
           "Number of boulder", "Black zone width (m)")
scat.plots(BS$GS, BS$Adj.LVD, BS, 
           "Number of boulder", "Land vegetation distance (m)")

scat.plots(BS$Adj.FW, BS$Adj.BZW, BS, 
           "Fucus width (m)", "Black zone width (m)")
scat.plots(BS$Adj.FW, BS$Adj.LVD, BS, 
           "Fucus width (m)", "Land vegetation distance (m)")

scat.plots(BS$Adj.FW, BS$Adj.LVD, BS, 
           "Black zone width (m)", "Land vegetation distance (m)")


bscor <- ggcorrplot(cor(BS %>% select(where(is.numeric))), type = "upper",
           p.mat = cor_pmat(BS %>% select(where(is.numeric))),
           title = "Boulder shores")


### Muddy-sandy shore correlations ####
scat.plots(log10(MSS$Exposure_FETCH), MSS$sin_aspect, MSS, 
           expression(Log[10](fetch)), "Shore aspect")
scat.plots(log10(MSS$Exposure_FETCH), MSS$Shore_length, MSS, 
           expression(Log[10](fetch)), "Shore length (m)")
scat.plots(log10(MSS$Exposure_FETCH), MSS$slope, MSS, 
           expression(Log[10](fetch)), "Transect slope (radians)")
scat.plots(log10(MSS$Exposure_FETCH), MSS$HTD, MSS, 
           expression(Log[10](fetch)), "Habitat transition distance (m)")
scat.plots(log10(MSS$Exposure_FETCH), MSS$FP, MSS, 
           expression(Log[10](fetch)), "Fucus patchiness (%)")
scat.plots(log10(MSS$Exposure_FETCH), MSS$Adj.FW, MSS, 
           expression(Log[10](fetch)), "Fucus belt width (m)")
scat.plots(log10(MSS$Exposure_FETCH), MSS$Adj.TZW, MSS, 
           expression(Log[10](fetch)), "Tide zone width (m)")
scat.plots(log10(MSS$Exposure_FETCH), MSS$Adj.LVD, MSS, 
           expression(Log[10](fetch)), "Land vegetation distance (m)")

scat.plots(MSS$sin_aspect, MSS$Shore_length, MSS, 
           "Shore aspect", "Shore length (m)")
scat.plots(MSS$sin_aspect, MSS$slope, MSS, 
           "Shore aspect", "Transect slope (radians)")
scat.plots(MSS$sin_aspect, MSS$HTD, MSS, 
           "Shore aspect", "Habitat transition distance (m)")
scat.plots(MSS$sin_aspect, MSS$FP, MSS, 
           "Shore aspect", "Fucus patchiness (%)")
scat.plots(MSS$sin_aspect, MSS$Adj.FW, MSS, 
           "Shore aspect", "Fucus belt width (m)")
scat.plots(MSS$sin_aspect, MSS$Adj.TZW, MSS, 
           "Shore aspect", "Tide zone width (m)")
scat.plots(MSS$sin_aspect, MSS$Adj.LVD, MSS, 
           "Shore aspect", "Land vegetation distance (m)")

scat.plots(MSS$Shore_length, MSS$slope, MSS, 
           "Shore length (m)", "Transect slope (radians)")
scat.plots(MSS$Shore_length, MSS$HTD, MSS, 
           "Shore length (m)", "Habitat transition distance (m)")
scat.plots(MSS$Shore_length, MSS$FP, MSS, 
           "Shore length (m)", "Fucus patchiness (%)")
scat.plots(MSS$Shore_length, MSS$Adj.FW, MSS, 
           "Shore length (m)", "Fucus belt width (m)")
scat.plots(MSS$Shore_length, MSS$Adj.TZW, MSS, 
           "Shore length (m)", "Tide zone width (m)")
scat.plots(MSS$Shore_length, MSS$Adj.LVD, MSS, 
           "Shore length (m)", "Land vegetation distance (m)")

scat.plots(MSS$slope, MSS$HTD, MSS, 
           "Transect slope (radians)", "Habitat transition distance (m)")
scat.plots(MSS$slope, MSS$FP, MSS, 
           "Transect slope (radians)", "Fucus patchiness (%)")
scat.plots(MSS$slope, MSS$Adj.FW, MSS, 
           "Transect slope (radians)", "Fucus belt width (m)")
scat.plots(MSS$slope, MSS$Adj.TZW, MSS, 
           "Transect slope (radians)", "Tide zone width (m)")
scat.plots(MSS$slope, MSS$Adj.LVD, MSS, 
           "Transect slope (radians)", "Land vegetation distance (m)")

scat.plots(MSS$HTD, MSS$FP, MSS, 
           "Habitat transition distance (m)", "Fucus patchiness (%)")
scat.plots(MSS$HTD, MSS$Adj.FW, MSS, 
           "Habitat transition distance (m)", "Fucus belt width (m)")
scat.plots(MSS$HTD, MSS$Adj.TZW, MSS, 
           "Habitat transition distance (m)", "Tide zone width (m)")
scat.plots(MSS$HTD, MSS$Adj.LVD, MSS, 
           "Habitat transition distance (m)", "Land vegetation distance (m)")

scat.plots(MSS$FP, MSS$Adj.FW, MSS, 
           "Fucus patchiness (%)", "Fucus belt width (m)")
scat.plots(MSS$FP, MSS$Adj.TZW, MSS, 
           "Fucus patchiness (%)", "Tide zone width (m)")
scat.plots(MSS$FP, MSS$Adj.LVD, MSS, 
           "Fucus patchiness (%)", "Land vegetation distance (m)")

scat.plots(MSS$Adj.FW, MSS$Adj.TZW, MSS, 
           "Fucus belt width (m)", "Tide zone width (m)")
scat.plots(MSS$Adj.FW, MSS$Adj.LVD, MSS, 
           "Fucus belt width (m)", "Land vegetation distance (m)")

scat.plots(MSS$Adj.TZW, MSS$Adj.LVD, MSS, 
           "Tide zone width (m)", "Land vegetation distance (m)")

msscor <- ggcorrplot(cor(MSS %>% select(where(is.numeric))), type = "upper",
           p.mat = cor_pmat(MSS %>% select(where(is.numeric))),
          title = "Muddy-sandy shores")

ggarrange(rscor, bscor, msscor, nrow = 1, common.legend = TRUE)


#### D: Run PCAmix ####

ShorePCA <- function(dat){
  # Convert site names to rownames
  rownames(dat) <- dat[,"Site"]
  tmp <- dat[,-c(1)]
  
  # Drop 'Exposure_notes'
  tmp$Exposure_notes <- NULL
  
  # Separate out continuous and categorical variables
  tmp2 <- splitmix(tmp)
  qual <- tmp2$X.quali
  quant <- tmp2$X.quanti
  
  # Rescale and recentre quantitative variables
  quant <- scale(quant, center = TRUE, scale = TRUE)
  
  # Run the PCA
  tmp_PCA <- PCAmix(X.quanti = quant, X.quali = qual,
                    rename.level = TRUE, graph = FALSE)
  
  return(tmp_PCA)
}

RS_PCA <- ShorePCA(RS)
BS_PCA <- ShorePCA(BS)
MSS_PCA <- ShorePCA(MSS)

### Format as dataframe
# Eigenvalues
RS_PCA_dat <- as.data.frame(RS_PCA$ind$coord)
BS_PCA_dat <- as.data.frame(BS_PCA$ind$coord)
MSS_PCA_dat <- as.data.frame(MSS_PCA$ind$coord)

# Variance explained
RS_Pvar <- as.data.frame(RS_PCA$eig)
RS_Pvar$dimension <- rownames(RS_PCA$eig)

BS_Pvar <- as.data.frame(BS_PCA$eig)
BS_Pvar$dimension <- rownames(BS_PCA$eig)

MSS_Pvar <- as.data.frame(MSS_PCA$eig)
MSS_Pvar$dimension <- rownames(MSS_PCA$eig)



#### E: Plot the results ####
# Default PCA plots
plot(RS_PCA, axes = c(1,2), choice = "ind", label = TRUE)
plot(RS_PCA, choice = "levels")
plot(RS_PCA, choice = "cor")
plot(RS_PCA, axes = c(1,2), choice = "sqload", coloring.var = T, leg=TRUE,
     posleg="none")

### Custom PCA plots with ggplot() ####
### Copy to clipboard: 1000x360 (w x h)
# Rocky shores
RS_loads <- RS_PCA$sqload
rownames(RS_loads) <- gsub("Adj.", "", rownames(RS_loads))

ggplot(RS_PCA_dat)+
  geom_point(aes(`dim 1`, `dim 2`, colour = log(GisD[GisD$Site %in% RS$Site, "Exposure_FETCH"])), 
             show.legend = TRUE, size = 3)+
  geom_segment(data = RS_loads, 
               aes(x = 0, y = 0, xend = `dim 1`*10, yend = `dim 2`*10), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "grey30")+
  geom_label(data = RS_loads, aes(x = `dim 1`*10, y = `dim 2`*10, label = rownames(RS_loads)),
             size = 2.6)+
  labs(x = paste0("Dim 1 (", round(RS_PCA$eig[1,"Proportion"], 2), "%)"), 
       y = paste0("Dim 2 (", round(RS_PCA$eig[2,"Proportion"], 2), "%)"), 
       title = "Rocky shores", 
       colour = expression(Log[e]~"(Fetch)"))+
  scale_colour_gradient(low = "#d3eaeaff", high = "#008080ff", 
                        limits = c(1.5, 13.5), breaks = seq(1.5, 13.5, 4))+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.margin = margin(l = 0.4, unit = "npc"))

# Boulder shores
BS_loads <- BS_PCA$sqload
rownames(BS_loads) <- gsub("Adj.", "", rownames(BS_loads))

ggplot(BS_PCA_dat)+
  geom_point(aes(`dim 1`, `dim 2`, colour = log(GisD[GisD$Site %in% BS$Site, "Exposure_FETCH"])), 
             show.legend = TRUE, size = 3)+
  geom_segment(data = BS_loads, 
               aes(x = 0, y = 0, xend = `dim 1`*10, yend = `dim 2`*10), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "grey30")+
  geom_label(data = BS_loads, aes(x = `dim 1`*10, y = `dim 2`*10, label = rownames(BS_loads)),
             size = 2.6)+
  labs(x = paste0("Dim 1 (", round(BS_PCA$eig[1,"Proportion"], 2), "%)"), 
       y = paste0("Dim 2 (", round(BS_PCA$eig[2,"Proportion"], 2), "%)"), 
       title = "Boulder fields", 
       colour = expression(Log[e]~"(Fetch)"))+
  scale_colour_gradient(low = "#d3eaeaff", high = "#008080ff", 
                         limits = c(1.5, 13.5), breaks = seq(1.5, 13.5, 4))+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.margin = margin(l = 0.4, unit = "npc"))

# Muddy and sandy shores
MSS_loads <- MSS_PCA$sqload
rownames(MSS_loads) <- gsub("Adj.", "", rownames(MSS_loads))

ggplot(MSS_PCA_dat)+
  geom_point(aes(`dim 1`, `dim 2`, colour = log(GisD[GisD$Site %in% MSS$Site, "Exposure_FETCH"])), 
             show.legend = TRUE, size = 3)+
  geom_segment(data = MSS_loads, 
               aes(x = 0, y = 0, xend = `dim 1`*8, yend = `dim 2`*8), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "grey30")+
  geom_label(data = MSS_loads, aes(x = `dim 1`*8, y = `dim 2`*8, label = rownames(MSS_loads)),
             size = 2.6)+
  labs(x = paste0("Dim 1 (", round(MSS_PCA$eig[1,"Proportion"], 2), "%)"), 
       y = paste0("Dim 2 (", round(MSS_PCA$eig[2,"Proportion"], 2), "%)"), 
       title = "Muddy-Sandy shores", 
       colour = expression(Log[e]~"(Fetch)"))+
  scale_colour_gradient(low = "#d3eaeaff", high = "#008080ff", 
                        limits = c(1.5, 13.5), breaks = seq(1.5, 13.5, 4))+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.margin = margin(l = 0.4, unit = "npc"))



### F: Scree plot of variance explained ####
### Copy to clipboard: 500x500 (w x h)

### Create broken stick estimate
RS_Pvar$BSexp <- sapply(1:nrow(RS_Pvar), function(i){ 
  sum(1/i:nrow(RS_Pvar))/nrow(RS_Pvar) * 100})
BS_Pvar$BSexp <- sapply(1:nrow(BS_Pvar), function(i){ 
  sum(1/i:nrow(BS_Pvar))/nrow(BS_Pvar) * 100})
MSS_Pvar$BSexp <- sapply(1:nrow(MSS_Pvar), function(i){ 
  sum(1/i:nrow(MSS_Pvar))/nrow(MSS_Pvar) * 100})


# Rocky shores
ggplot(RS_Pvar, aes(x = dimension))+
  geom_line(aes(y = Proportion, group = ""), col = "black")+
  geom_point(aes(y = Proportion), col = "black")+
  geom_point(aes(y = BSexp), col = "grey80")+
  geom_line(aes(y = BSexp, group = 1), col = "grey80")+
  labs(x = "Principal component", y = "Varaince explained (%)")+#, title = "Rocky shores")+
  ylim(c(0, 30))+
  scale_x_discrete(limits = paste0("dim ", 1:nrow(RS_PCA$eig)))+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 16, face = "bold"))

# Boulder shores
ggplot(BS_Pvar, aes(x = dimension))+
  geom_line(aes(y = Proportion, group = ""), col = "black")+
  geom_point(aes(y = Proportion), col = "black")+
  geom_point(aes(y = BSexp), col = "grey80")+
  geom_line(aes(y = BSexp, group = 1), col = "grey80")+
  labs(x = "Principal component", y = "Varaince explained (%)")+#, title = "Boulder fields")+
  ylim(c(0, 30))+
  scale_x_discrete(limits = paste0("dim ", 1:nrow(BS_PCA$eig)))+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 16, face = "bold"))

# Muddy and sandy shores
ggplot(MSS_Pvar, aes(x = dimension))+
  geom_line(aes(y = Proportion, group = ""), col = "black")+
  geom_point(aes(y = Proportion), col = "black")+
  geom_point(aes(y = BSexp), col = "grey80")+
  geom_line(aes(y = BSexp, group = 1), col = "grey80")+
  labs(x = "Principal component", y = "Varaince explained (%)")+#, title = "Muddy & Sandy shores")+
  ylim(c(0, 30))+
  scale_x_discrete(limits = paste0("dim ", 1:nrow(MSS_PCA$eig)))+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 16, face = "bold"))


### G: Histogram of first 5 dimensions of PCAmix ####
### Copy to clipboard: 360 x 500 (w x h)
PCA.hist <- function(data, principal.compent, binwidth = 0.2,
                     highlight.colour = "grey50"){
  
  Xlim <- c(floor(min(sapply(data, range))),
            ceiling(max(sapply(data, range)))) # Set plot limits based on max and min PC scores
  Xbreaks <- seq(Xlim[1], Xlim[2], 1)               # Set breaks based on 'Xlim'
  
  ggplot(data, aes(x = data[, principal.compent]))+
    geom_histogram(binwidth = binwidth, colour = "grey50", fill = NA, alpha = 0.4)+
    geom_density(aes(y = binwidth * ..count..), alpha = 0.4, col = highlight.colour)+
    labs(y = paste("Dim", principal.compent))+
    scale_x_continuous(breaks = Xbreaks, limits = Xlim)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text.y = element_blank())
}

# Rocky shore
p1a <- PCA.hist(RS_PCA_dat, principal.compent = 1, binwidth = 0.4,
         highlight.colour = "#888EEB")
p1b <- PCA.hist(RS_PCA_dat, principal.compent = 2, binwidth = 0.4,
         highlight.colour = "#74EBA3")
p1c <- PCA.hist(RS_PCA_dat, principal.compent = 3, binwidth = 0.4,
         highlight.colour = "#D48CE6")
p1d <- PCA.hist(RS_PCA_dat, principal.compent = 4, binwidth = 0.3,
         highlight.colour = "#E6D38A")
p1e <- PCA.hist(RS_PCA_dat, principal.compent = 5, binwidth = 0.3,
         highlight.colour = "#EB817A")
ggarrange(p1a, p1b, p1c, p1d, p1e, ncol = 1)

# Boulder shores
p1a <- PCA.hist(BS_PCA_dat, principal.compent = 1, binwidth = 0.4,
                highlight.colour = "#888EEB")
p1b <- PCA.hist(BS_PCA_dat, principal.compent = 2, binwidth = 0.4,
                highlight.colour = "#74EBA3")
p1c <- PCA.hist(BS_PCA_dat, principal.compent = 3, binwidth = 0.4,
                highlight.colour = "#D48CE6")
p1d <- PCA.hist(BS_PCA_dat, principal.compent = 4, binwidth = 0.3,
                highlight.colour = "#E6D38A")
p1e <- PCA.hist(BS_PCA_dat, principal.compent = 5, binwidth = 0.3,
                highlight.colour = "#EB817A")
ggarrange(p1a, p1b, p1c, p1d, p1e, ncol = 1)

# Muddy and sandy shores
p1a <- PCA.hist(MSS_PCA_dat, principal.compent = 1, binwidth = 0.4,
                highlight.colour = "#888EEB")
p1b <- PCA.hist(MSS_PCA_dat, principal.compent = 2, binwidth = 0.4,
                highlight.colour = "#74EBA3")
p1c <- PCA.hist(MSS_PCA_dat, principal.compent = 3, binwidth = 0.4,
                highlight.colour = "#D48CE6")
p1d <- PCA.hist(MSS_PCA_dat, principal.compent = 4, binwidth = 0.3,
                highlight.colour = "#E6D38A")
p1e <- PCA.hist(MSS_PCA_dat, principal.compent = 5, binwidth = 0.3,
                highlight.colour = "#EB817A")
ggarrange(p1a, p1b, p1c, p1d, p1e, ncol = 1)



### H: Plot top 5 dimensions in large grid ####
### Copy to clipboard: 500 x 500 (w x h)
PCA.scat <- function(data, x, y){
  # Select exposure data to colour plot by
  #if(identical(data, RS_PCA$ind$coord)){
    #Exposure <-  FB_Site[FB_Site$Site %in% RS$Site, "Exposure_notes"]
    #colPal <- c("#7EB6E0", "#EBB76A", "#025EA3", "#A36302")}
  #if(identical(data, BS_PCA$ind$coord)){
    #Exposure <-  FB_Site[FB_Site$Site %in% BS$Site, "Exposure_notes"]
    #colPal <- c("#7EB6E0", "#EBB76A", "#A36302")}
  #if(identical(data, MSS_PCA$ind$coord)){
    #Exposure <-  FB_Site[FB_Site$Site %in% MSS$Site, "Exposure_notes"]
    #colPal <- c("#EBB76A", "#A36302")}
  
  ggplot(as.data.frame(data))+
    geom_point(aes(x = data[, x], y = data[, y],
                   colour = log(GisD[GisD$Site %in% rownames(data), "Exposure_FETCH"])), 
               show.legend = FALSE)+
    scale_x_continuous(breaks = 0)+
    scale_y_continuous(breaks = 0)+
    #scale_colour_manual(values = colPal)+
    scale_colour_gradient(low = "#d3eaeaff", high = "#008080ff")+
    labs(x = paste("Dim", x), y = paste("Dim", y))+
    theme_classic()+
    theme(panel.grid.major = element_line(colour = "lightgrey"),
          axis.title = element_text(size = 8))
}

# Rocky shores
p2a <- PCA.scat(RS_PCA$ind$coord, 1, 2)
p2b <- PCA.scat(RS_PCA$ind$coord, 1, 3)
p2c <- PCA.scat(RS_PCA$ind$coord, 1, 4)
p2d <- PCA.scat(RS_PCA$ind$coord, 1, 5)
p2e <- PCA.scat(RS_PCA$ind$coord, 2, 3)
p2f <- PCA.scat(RS_PCA$ind$coord, 2, 4)
p2g <- PCA.scat(RS_PCA$ind$coord, 2, 5)
p2h <- PCA.scat(RS_PCA$ind$coord, 3, 4)
p2i <- PCA.scat(RS_PCA$ind$coord, 3, 5)
p2j <- PCA.scat(RS_PCA$ind$coord, 4, 5)
p2_ <- ggplot()+theme_void()

ggarrange(ggarrange(p2a, p2b, p2c, p2d, nrow = 1),
          ggarrange(p2_, p2e, p2f, p2g, nrow = 1),
          ggarrange(p2_, p2_, p2h, p2i, nrow = 1),
          ggarrange(p2_, p2_, p2_, p2j, nrow = 1),
          nrow = 4, ncol = 1)

# Boulder shores
p2a <- PCA.scat(BS_PCA$ind$coord, 1, 2)
p2b <- PCA.scat(BS_PCA$ind$coord, 1, 3)
p2c <- PCA.scat(BS_PCA$ind$coord, 1, 4)
p2d <- PCA.scat(BS_PCA$ind$coord, 1, 5)
p2e <- PCA.scat(BS_PCA$ind$coord, 2, 3)
p2f <- PCA.scat(BS_PCA$ind$coord, 2, 4)
p2g <- PCA.scat(BS_PCA$ind$coord, 2, 5)
p2h <- PCA.scat(BS_PCA$ind$coord, 3, 4)
p2i <- PCA.scat(BS_PCA$ind$coord, 3, 5)
p2j <- PCA.scat(BS_PCA$ind$coord, 4, 5)
p2_ <- ggplot()+theme_void()

ggarrange(ggarrange(p2a, p2b, p2c, p2d, nrow = 1),
          ggarrange(p2_, p2e, p2f, p2g, nrow = 1),
          ggarrange(p2_, p2_, p2h, p2i, nrow = 1),
          ggarrange(p2_, p2_, p2_, p2j, nrow = 1),
          nrow = 4, ncol = 1)

# Muddy and sandy shores
p2a <- PCA.scat(MSS_PCA$ind$coord, 1, 2)
p2b <- PCA.scat(MSS_PCA$ind$coord, 1, 3)
p2c <- PCA.scat(MSS_PCA$ind$coord, 1, 4)
p2d <- PCA.scat(MSS_PCA$ind$coord, 1, 5)
p2e <- PCA.scat(MSS_PCA$ind$coord, 2, 3)
p2f <- PCA.scat(MSS_PCA$ind$coord, 2, 4)
p2g <- PCA.scat(MSS_PCA$ind$coord, 2, 5)
p2h <- PCA.scat(MSS_PCA$ind$coord, 3, 4)
p2i <- PCA.scat(MSS_PCA$ind$coord, 3, 5)
p2j <- PCA.scat(MSS_PCA$ind$coord, 4, 5)
p2_ <- ggplot()+theme_void()

ggarrange(ggarrange(p2a, p2b, p2c, p2d, nrow = 1),
          ggarrange(p2_, p2e, p2f, p2g, nrow = 1),
          ggarrange(p2_, p2_, p2h, p2i, nrow = 1),
          ggarrange(p2_, p2_, p2_, p2j, nrow = 1),
          nrow = 4, ncol = 1)



### I: Plot variable loadings for first 5 dimensions ####
### Copy to clipboard: 1000 x 250 (w x h)
PCA.load <- function(data, principal.component, highlight.colour = "grey50"){
  tmp <- as.data.frame(data$sqload)
  tmp$var <- rownames(tmp)
  
  ggplot(tmp)+
    #geom_point(aes(x = tmp[, principal.component], y = var), 
     #          pch = '>', size = 3, col = highlight.colour)+
    #geom_col(aes(x = tmp[, principal.component], y = var), 
     #        width = 0.05, col = highlight.colour)+
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.1)+
    geom_segment(aes(x = 0, y = var, xend = tmp[, principal.component], yend = var), 
                 arrow = arrow(length = unit(0.8, "picas")),
                 color = highlight.colour, linewidth = 1.2)+
    labs(x = paste("Dim", principal.component))+
    lims(x = c(0,1))+
    theme_classic()+
    theme(axis.title.y = element_blank())
}


# Rocky shores
p3a <- PCA.load(RS_PCA, 1, highlight.colour = "#888EEB")
p3b <- PCA.load(RS_PCA, 2, highlight.colour = "#74EBA3")+theme(axis.text.y = element_blank())
p3c <- PCA.load(RS_PCA, 3, highlight.colour = "#D48CE6")+theme(axis.text.y = element_blank())
p3d <- PCA.load(RS_PCA, 4, highlight.colour = "#E6D38A")+theme(axis.text.y = element_blank())
p3e <- PCA.load(RS_PCA, 5, highlight.colour = "#EB817A")+theme(axis.text.y = element_blank())

ggarrange(p3a, p3b, p3c, #p3d, p3e, 
          nrow = 1, widths = c(1.5,1,1))

# Boulder shores
p3a <- PCA.load(BS_PCA, 1, highlight.colour = "#888EEB")
p3b <- PCA.load(BS_PCA, 2, highlight.colour = "#74EBA3")+theme(axis.text.y = element_blank())
p3c <- PCA.load(BS_PCA, 3, highlight.colour = "#D48CE6")+theme(axis.text.y = element_blank())
p3d <- PCA.load(BS_PCA, 4, highlight.colour = "#E6D38A")+theme(axis.text.y = element_blank())
p3e <- PCA.load(BS_PCA, 5, highlight.colour = "#EB817A")+theme(axis.text.y = element_blank())

ggarrange(p3a, p3b, p3c, #p3d, p3e, 
          nrow = 1, widths = c(1.5,1,1))

# Muddy and sandy shores
p3a <- PCA.load(MSS_PCA, 1, highlight.colour = "#888EEB")
p3b <- PCA.load(MSS_PCA, 2, highlight.colour = "#74EBA3")+theme(axis.text.y = element_blank())
p3c <- PCA.load(MSS_PCA, 3, highlight.colour = "#D48CE6")+theme(axis.text.y = element_blank())
p3d <- PCA.load(MSS_PCA, 4, highlight.colour = "#E6D38A")+theme(axis.text.y = element_blank())
p3e <- PCA.load(MSS_PCA, 5, highlight.colour = "#EB817A")+theme(axis.text.y = element_blank())

ggarrange(p3a, p3b, p3c, #p3d, p3e, 
          nrow = 1, widths = c(1.5,1,1))


### F: Evaluate number of dimensions to save using PCAtest ####

### Note: this package calculates statistics to evaluate observed PCA
### against a null distribution of permutations. I am not focusing on the 
### overall "significance" of the PCA, instead I am using this to look at
### which dimensions explain more variation than expect by chance.

## Rocky shores
PCAtest(RS_PCA$Z)


#### Currently not working!!! I think this is a scalling problem with the 
### package. I've made an error request on the GitHub.


#### G: Save PCAmix results ####

### Saving principal components (top 5 dimensions)
PATH.1 <- "/output/dir/principal_components/"

write.csv(RS_PCA$ind$coord, paste0(PATH.1, "Env_RS_PC1-5.csv"), quote = FALSE)
write.csv(BS_PCA$ind$coord, paste0(PATH.1, "Env_BS_PC1-5.csv"), quote = FALSE)
write.csv(MSS_PCA$ind$coord, paste0(PATH.1, "Env_MSS_PC1-5.csv"), quote = FALSE)

### Save eigenvalues and variance explained
PATH.2 <- "/output/dir/eigenvalues/"

write.csv(RS_PCA$eig, paste0(PATH.2, "Env_RS_eig.csv"), quote = FALSE)
write.csv(BS_PCA$eig, paste0(PATH.2, "Env_BS_eig.csv"), quote = FALSE)
write.csv(MSS_PCA$eig, paste0(PATH.2, "Env_MSS_eig.csv"), quote = FALSE)



#### H: Addendum: testing skew and modality of data ####
### Dim 1 in rocky shores looks to be bi-modal. While other dimensions look to be non-normally
### distributed. These are investigated by looking at the skewness, kurtosis (testing for normality), 
### and unimodality of the first three dimensions.

### Packages
library(moments) # to calcaulte skewness and kurtosis
library(diptest) # Runs Hartigan's dip test for unimodality

### Skewness ( > 0 = postive skew: < 0 = negative skew)
skewness(RS_PCA$dim.1)
skewness(RS_PCA$dim.2)
skewness(RS_PCA$dim.3)

skewness(BS_PCA$dim.1)
skewness(BS_PCA$dim.2)
skewness(BS_PCA$dim.3)

skewness(MSS_PCA$dim.1)
skewness(MSS_PCA$dim.2)
skewness(MSS_PCA$dim.3)

### Kurtosis ( > 0 wider tails than normal distribution; < 0  narrower tails)
kurtosis(RS_PCA$dim.1)
kurtosis(RS_PCA$dim.2)
kurtosis(RS_PCA$dim.3)

kurtosis(BS_PCA$dim.1)
kurtosis(BS_PCA$dim.2)
kurtosis(BS_PCA$dim.3)

kurtosis(MSS_PCA$dim.1)
kurtosis(MSS_PCA$dim.2)
kurtosis(MSS_PCA$dim.3)

### Hartigan's test for unimodality
dip.test(RS_PCA_dat$`dim 1`)
dip.test(RS_PCA_dat$`dim 2`)
dip.test(RS_PCA_dat$`dim 3`)

dip.test(BS_PCA_dat$`dim 1`)
dip.test(BS_PCA_dat$`dim 2`)
dip.test(BS_PCA_dat$`dim 3`)

dip.test(MSS_PCA_dat$`dim 1`)
dip.test(MSS_PCA_dat$`dim 2`)
dip.test(MSS_PCA_dat$`dim 3`)
