########### Select SNP markers for genetic background ##################
### This script selects ~600 markers which represent the genetic background
### of Littorina saxatilis in Sweden. These SNPs were selected by randomly
### picking map positions and then randomly selecting one SNP from each
### map position.

### James Reeve - University of Gothenburg
### 09/12/2021

### Preparation
rm(list = ls())
dev.off()
setwd("~")
PATH <- "/path/to/working/directory/"

### Packages
library(tidyverse)
library(ggpubr)

### Parameters
N <- 449 # Number of SNPs to return!
buffer <- 2 # Buffer around putative breakpoints (cM)

### Access data
# Linkage map (crab map)
LMap <- read.table("/Users/james/Documents/Inversion_detection/Map_filtrated.txt", header = TRUE)
LMap$LG <- paste0("LG", LMap$LG)
# Inversion breakpoints
InvPos <- read.csv("/Users/james/Dropbox (Personal)/PhD/Inversion_positions_on_old_map.csv", header = T)

### A: Filter out inversions from linkage map ####
# Add contig column to LMap
#LMap$Contig <- gsub("_.*", "", LMap$Marker) ### Redundant for crab map

# Make list of linkage groups
LGs <- unique(LMap$LG)

# Filter inversion positions to just consensus positions
InvPos <- InvPos[InvPos$Site == "consensus",]

# Filter out markers inside inversion
LMap.CL <- lapply(LGs, function(LG){
  # Subset inversion positions to given linkage group
  LGinvs <- InvPos[InvPos$LG == LG,]
  # Number of inversions on LG
  Ninv <- nrow(LGinvs)/2
  # Subset the linkage map
  LMap.LG <- LMap[LMap$LG == LG,]
  
  # Return unfiltered map if no inversions are present
  if(Ninv == 0) return(LMap.LG)
  
  # 1 inversion present
  if(Ninv == 1) {
    Spos <- LGinvs[LGinvs$Start == T, "Pos"][1] - buffer
    Epos <- LGinvs[LGinvs$Start == F, "Pos"][1] + buffer
    subMap <- LMap.LG[!(LMap.LG$av >= Spos & LMap.LG$av <= Epos), ]
    return(subMap)
  }
  
  # 2 inversion present
  if(Ninv == 2) {
    Spos.1 <- LGinvs[LGinvs$Start == T, "Pos"][1] - buffer
    Epos.1 <- LGinvs[LGinvs$Start == F, "Pos"][1] + buffer
    Spos.2 <- LGinvs[LGinvs$Start == T, "Pos"][2] - buffer
    Epos.2 <- LGinvs[LGinvs$Start == F, "Pos"][2] + buffer
    subMap <- LMap.LG[!(LMap.LG$av >= Spos.1 & LMap.LG$av <= Epos.1 | LMap.LG$av >= Spos.2 & LMap.LG$av <= Epos.2), ]
    return(subMap)
  }
  
  # 3 inversions present
  if(Ninv == 3) {
    Spos.1 <- LGinvs[LGinvs$Start == T, "Pos"][1] - buffer
    Epos.1 <- LGinvs[LGinvs$Start == F, "Pos"][1] + buffer
    Spos.2 <- LGinvs[LGinvs$Start == T, "Pos"][2] - buffer
    Epos.2 <- LGinvs[LGinvs$Start == F, "Pos"][2] + buffer
    Spos.3 <- LGinvs[LGinvs$Start == T, "Pos"][3] - buffer
    Epos.3 <- LGinvs[LGinvs$Start == F, "Pos"][3] + buffer
    subMap <- LMap.LG[!(LMap.LG$av >= Spos.1 & LMap.LG$av <= Epos.1 | LMap.LG$av >= Spos.2 & LMap.LG$av <= Epos.2 | 
                        LMap.LG$av >= Spos.3 & LMap.LG$av <= Epos.3), ]
    return(subMap)
  }
  
  # 4 inversions present
  if(Ninv == 4) {
    Spos.1 <- LGinvs[LGinvs$Start == T, "Pos"][1] - buffer
    Epos.1 <- LGinvs[LGinvs$Start == F, "Pos"][1] + buffer
    Spos.2 <- LGinvs[LGinvs$Start == T, "Pos"][2] - buffer
    Epos.2 <- LGinvs[LGinvs$Start == F, "Pos"][2] + buffer
    Spos.3 <- LGinvs[LGinvs$Start == T, "Pos"][3] - buffer
    Epos.3 <- LGinvs[LGinvs$Start == F, "Pos"][3] + buffer
    Spos.4 <- LGinvs[LGinvs$Start == T, "Pos"][4] - buffer
    Epos.4 <- LGinvs[LGinvs$Start == F, "Pos"][4] + buffer
    subMap <- LMap.LG[!(LMap.LG$av >= Spos.1 & LMap.LG$av <= Epos.1 | LMap.LG$av >= Spos.2 & LMap.LG$av <= Epos.2 |
                         LMap.LG$av >= Spos.3 & LMap.LG$av <= Epos.3 | LMap.LG$av >= Spos.4 & LMap.LG$av <= Epos.4), ]
    return(subMap)
  }
  # No linkage group has more than 4 inversions in L.saxatilis
  # however, I will need to get back to this code and make it robust
  # to an indefinitely large number of inversions!
})
LMap.CL <- do.call(rbind.data.frame, LMap.CL)

# Check that inversion regions are missing with ggplot
Logi <- LMap$av %in% LMap.CL$av
ggplot(LMap, aes(x = av, colour = Logi))+
  geom_bar()+
  facet_wrap(vars(LG), nrow = 3, ncol = 6)+
  labs(x = "Map position (cM)")+
  theme(legend.position = "none")

### B: Filter SNPs in WGS data to those on the collinear linkage map ####
SNPs <- lapply(LGs, function(LG){
  # Access [filtered] SNP list for specific linkage group
  tmp <-   read.table(paste0(PATH, "SNP_list/", LG, "_SNPs.txt"), col.names = c("contig", "pos"))
  # Merge contig and position to get SNP ID
  tmp$SNP <- paste(tmp$contig, tmp$pos, sep = "_")
  return(tmp)
})
SNPs <- do.call(rbind.data.frame, SNPs)  

# Filtering step
SNPs <- SNPs[SNPs$contig %in% LMap.CL$contig,]

### C: Select SNPs randomly ####
# Add map position to SNPs
GBsnps <- SNPs %>% group_by(contig) %>%
  left_join(LMap.CL, by = "contig") %>% # Add map positions to SNPs
  summarise(SNP, LG, av) %>%
  # Select 1 random SNP per map positions
  group_by(LG, av) %>%
  slice_sample(n = 1) %>%
  # Select N random SNPs
  ungroup() %>%
  slice_sample(n = N)

# Save data
write.table(GBsnps, paste0(PATH, "Genetic_background_SNPs_N", N,  "_", buffer, "cM_buffer", ".txt"), row.names = FALSE, quote = FALSE)
