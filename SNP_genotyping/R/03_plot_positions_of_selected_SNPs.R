########## Plot first round SNP panel onto the linkage map ##############
### Plotting the position of SNP panel candidates onto the Littorina 
### saxatilis linakge map. This is a sanity check to make sure SNP panel
### candidates are even spread across the reference genome, and that 
### genetic background markers are not landing in inversions.

### James Reeve - University of Gothenburg
### 10/12/2021

### Preparation
rm(list=ls())
dev.off()
options(stringsAsFactors = FALSE)

#### 1: Access data ####
# Linkage map
LMap <- read.table("/path/to/linakge/map_v11.txt", header = TRUE)
LMap$LG <-  paste0("LG", LMap$LG)

# SNP panel 
SNP_panel <- read.table("/path/to/Koster_SNP_panel_2022-02-11.bed", 
                        col.names = c("CONTIG", "START", "END", "NAME"))
SNP_panel <- SNP_panel[complete.cases(SNP_panel),]

# Inversion positions
Invs <- read.csv("/path/to/Inversion_positions_on_old_map.csv")


#### 2: Data wrangeling ####
# Split up SNP names to get more information
tmp <- strsplit(SNP_panel$NAME, split = "_")
tmp <- data.frame("LG" = sapply(tmp, `[[`, 3),
                  "map.pos" = as.numeric(gsub("cM", "", sapply(tmp, `[[`, 4))),
                  "dataset" = paste0(sapply(tmp, `[[`, 5),"_", sapply(tmp, `[[`, 6)))

# Add new column for dataset
tmp[tmp$dataset == "Paralleogram_inversion", "dataset"] <- "ipr"
tmp[tmp$dataset == "colinear_outliers", "dataset"] <- "cot"
tmp[tmp$dataset == "colinear_genetic", "dataset"] <- "gbg"
tmp[tmp$dataset == "Koster_inversion", "dataset"] <- "ikt"
tmp[tmp$dataset == "Koster_background", "dataset"] <- "ikt" # This is LGC14.3 markers

# Subset inversion positions
Invs <- Invs[Invs$Site == "consensus",]
Invs$Start_pos <- Invs[Invs$Start == TRUE, "Pos"]
Invs$End_pos <- Invs[Invs$Start == FALSE, "Pos"]
Invs <- Invs %>% select(-Pos, -Start) %>% distinct()


#### 3: Plot ####
ggplot()+
  geom_rect(data = Invs, aes(ymin = as.numeric(gsub("LG", "", LG))-0.3, 
                             ymax = as.numeric(gsub("LG", "", LG))+0.3, 
                             xmin = Start_pos, xmax = End_pos),
            fill = "red", alpha = 0.4, col = "darkred")+
  geom_point(data = LMap, aes(x = av, y = as.numeric(gsub("LG", "", LG))), 
             col = "grey50", alpha = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = tmp, aes(x = map.pos, y = as.numeric(gsub("LG", "", LG)), colour = dataset),
             size = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = c("red", "blue", "green", "gold"))+
  scale_y_continuous(breaks = 1:17)+
  labs(x = "Position on crab map (cM)", y = "Linkage group", colour = "Dataset")+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour = "grey", size = 0.5),
        panel.grid.major.y = element_blank())
