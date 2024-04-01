########################################## Select SNPs for SNP-panel Round 2 ########################################
### This script identifies SNPs to include in a SNP-panel which will detect inversions in Swedish Littorina saxatilis.
### There are five different types of SNP marker which will be included in this panel:
###       1. Koster WGS inversion markers [Eva Koch & James Reeve 2021]
###       2. 10 Random SNPs from LGC14.3 [James Reeve 2022]
###       3. Inversion diagnostic markers from Parrallelogram approach [Westram et al. 2023, bioRxiv]
###       4. Outlier SNPs outside inversions from CZ-cline analysis [Westram et al. 2018, Evolution Letters]
###       5. Koster WGS genetic background markers [James Reeve 2021]
### This script refines the results of the first round of marker selection. An original set of ~1000 markers was 
### genotyped in 94 individuals. 856 markers passed genotyping. I will further refine these markers to select 364 
### markers for the second round of genotyping.

### James Reeve - University of Gothenburg
### 27/06/2022

### Preparation
rm(list = ls())
dev.off()
setwd("~")
PATH <- "/path/to/SNP/panel/results/"

### Packages
library("tidyverse")
library("ggpubr")
library("reshape")


#### 1: Access and format data ####

### A: Read SNP data results
# Summary data
KTsnp <- read.csv(paste0(PATH, "Genotyping-3366.015-03.csv"), skip = 14)[1:855,]

## Access information from SNP IDs
tmp <- strsplit(KTsnp$SNP, "_")
# Note: decimal points in SNP IDs have been replaced by "_"
# this loop finds "_" and replaces them with "." in the correct places
tmp2 <- lapply(tmp, function(X){
  # Fix "ikt14_3"
  if(length(X) == 7){
    X[4] <- paste(X[4], X[5], sep = ".")
    X[6] <- paste(X[6], X[7], sep = ".")
    X <- X[-c(5,7)]
  } else if(length(X) == 6 && X[6] == "3"){
    X[5] <- paste(X[5], X[6], sep = ".")
    X <- X[-c(6)]
    # Fix "18_9cM"
  } else if(length(X) == 6){
    X[4] <- paste(X[4], X[5], sep = ".")
    X <- X[-c(5)]
  } else {
    X
  }
})

# Save SNP info as new columns
KTsnp$Contig <- sapply(tmp, `[[`, 1) # contig
KTsnp$LG <- sapply(tmp, `[[`, 3) # linkage group
KTsnp$MapPos <- sapply(tmp2, `[[`, 4) # map position
KTsnp$MapPos <- as.numeric(gsub("cM", "", KTsnp$MapPos)) # convert to numeric
KTsnp$Source <- sapply(tmp2, `[[`, 5) # source of data

### Genotypes
KTsnpGT <- read.csv(paste0(PATH, "Genotyping-3366.015-03 Grid.csv"), skip = 7)[-c(1),] # Note 1st row is control marker

# Move sample IDs to rownames
rownames(KTsnpGT) <- KTsnpGT$DNA...Assay
KTsnpGT <- KTsnpGT %>% select(-DNA...Assay)

# Convert genotypes to 0/0, 0/1 & 1/1
for(i in 1:ncol(KTsnpGT)){
  snp <- KTsnpGT[,i]
  
  Ref <- KTsnp$Allele.X[i]
  Alt <- KTsnp$Allele.Y[i]
  
  KTsnpGT[snp == paste(Ref, Ref, sep = ":"), i] <- 0
  KTsnpGT[snp == paste(Alt, Alt, sep = ":"), i] <- 2
  KTsnpGT[snp == paste(Ref, Alt, sep = ":") | 
            snp == paste(Alt, Ref, sep = ":"), i] <- 1
}

# Change "?", "Uncallable" & "Bad" to NA
KTsnpGT[KTsnpGT == "?" | KTsnpGT == "Uncallable" | KTsnpGT == "Bad"] <- NA

# Convert genotypes to numeric
KTsnpGT <- as.data.frame(sapply(KTsnpGT, as.numeric))


### B: Inversion positions
Invs <- read.csv("/path/to/Inversion_positions_on_old_map.csv")
Invs <- Invs[Invs$Site == "consensus",]
Invs$Start_pos <- Invs[Invs$Start == TRUE, "Pos"]
Invs$End_pos <- Invs[Invs$Start == FALSE, "Pos"]
Invs <- Invs %>% select(-Pos, -Start) %>% distinct()


#### 2: General filters ####
### This are universal filters that can be applied to all types of data

### A: Filter out Unused samples < 5
KTsnp2 <- KTsnp[KTsnp$Unused <= 5,]

### B: Filter out fixed markers
# Make list of SNPs
snp_list <- unlist(KTsnp2$SNP)
# Search columns of genotype data to find SNPs with more than 1 genotype (excluding NAs)
logi <- apply(KTsnpGT[,snp_list], 2, function(x) {length(unique(na.omit(x))) > 1} )

KTsnp3 <- KTsnp2[logi, ]

### C: Find any cases where SNP are sharing a contig
# Randomly keep 1 of the SNPs
KTsnp4 <- KTsnp3 %>% group_by(LG, Contig) %>%
  slice_sample(n = 1)

#### 3: Inversion marker filters ####
### Inversion markers were filtered to retain only those most strongly correlated to the inversion
### The randomly selected markers on LGC14.3 were included

### A: Function to calcualte to make correlation matrix between SNPs in specific inversions
inversion.LD <- function(genotype.data, inversion, SNP.type){
  # Get inversion details
  LG <- Invs[Invs$INV == inversion, "LG"]
  Spos <- Invs[Invs$INV == inversion, "Start_pos"]
  Epos <- Invs[Invs$INV == inversion, "End_pos"]
  
  # Find SNPs in the inversion
  snp_list <- unlist(KTsnp4[KTsnp4$LG == LG &
                              KTsnp4$MapPos >= Spos &
                              KTsnp4$MapPos <= Epos, "SNP"])
  snp_list <- snp_list[grepl(SNP.type, snp_list)]
  
  # Filter genotypes
  invGT <- genotype.data[,snp_list]
  
  # Calculate correlation matrix
  invLD <- cor(invGT, use = "pairwise.complete.obs") # skips NAs in calculation
  
  # Remove impossible correlations from genotypes
  # These represent cases when a SNP is fixed for a genotype.
  logi <- apply(invLD, 2, function(snp){ any(complete.cases(snp)) })
  Res <- invLD[logi,logi]
  
  return(Res)
}

### B: Run function
KTinvLD <- sapply(Invs$INV, function(inv) { try(inversion.LD(KTsnpGT, inv, "ipr|ikt")) })

# Plot distribution of correlation scores
annotate_figure(ggarrange(plotlist = lapply(names(KTinvLD), function(x){
  ggplot(melt(KTinvLD[[x]]), aes(x = abs(value)))+
    geom_histogram(aes(y = ..count..), breaks = seq(0, 1, 0.05))+
    geom_density(aes(y = 0.05*..count..), 
                 colour = "purple", fill = "purple", alpha = 0.4)+
    labs(title = paste0(x,": (Nsnp =", nrow(KTinvLD[[x]]), ")"))+
    xlim(c(0,1))+
    theme_classic()+
    theme(axis.title = element_blank())
})), left = "SNP count", bottom = "Linkage disequilibirum (|r|)")

# Plot correlation matrixes
ggarrange(plotlist = lapply(names(KTinvLD), function(x){
  ggplot(reshape::melt(KTinvLD[[x]]))+
    geom_tile(aes(x = X1, y = X2, fill = value))+
    scale_fill_gradient2(low = "#8856a7", high = "#8856a7", 
                         mid = "#e0ecf4", na.value = "white")+
    ggtitle(x)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = 2),
          axis.text.x = element_blank())
}), common.legend = TRUE, legend = "right")

### C: Filter SNPs by an LD threshold
## Find N snps with tightest LD
filter.by.LD <- function(cor.matrix, max.num.snps){
  tmp <- cor.matrix
  if(nrow(tmp) < max.num.snps){stop("Too few SNPs for filtering")}
  
  # Set all autocorrelations to 1
  diag(tmp) <- NA
  
  # Calculate average LD for each SNP
  avgR <- apply(abs(tmp), 1, mean, na.rm = T)
  
  # Identify N snps with the highest average correlations
  logi <- avgR >= sort(avgR)[length(avgR)-(max.num.snps-1)]
  
  # Filter original data frame
  Res <- cor.matrix[logi,logi]
  
  return(Res)
}

## Set threshold based on the 5th best average correlation score
KTinvLD2 <- lapply(names(KTinvLD), function(inv){
  if(inv == "LGC6.1/2") return(filter.by.LD(KTinvLD[[inv]], max.num.snps = 7)) else
    if(inv == "LGC1.1" | inv == "LGC9.1") return(filter.by.LD(KTinvLD[[inv]], max.num.snps = 6)) else
    if(inv == "LGC14.1/2" | inv == "LGC14.3") return(KTinvLD[[inv]]) else
      return(filter.by.LD(KTinvLD[[inv]], max.num.snps = 5))
  })
names(KTinvLD2) <- names(KTinvLD)

# Plot correlations after filtering
# Note: when all values defult to 1 ggplot treats this as 0.
ggarrange(plotlist = lapply(names(KTinvLD2), function(x){
  ggplot(melt(KTinvLD2[[x]]))+
    geom_tile(aes(x = X1, y = X2, fill = value))+
    scale_fill_gradient2(low = "#8856a7", high = "#8856a7", 
                         mid = "#e0ecf4", na.value = "white")+
    ggtitle(x)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 2))
}), common.legend = TRUE, legend = "right")

## Make list of best inversion markers
snp_list_invMark <- unlist(lapply(KTinvLD2, function(x){rownames(x)}))


### D: Additional manual filters
# Manual edit: remove some SNPs which still had poor correlations based on correlation plots
snp_list_invMark  <- snp_list_invMark[!(snp_list_invMark  %in% c("Contig41294_34954_LG14_0_4cM_ikt", 
                                       "Contig39483_65696_LG12_48_7cM_ikt"))]

# Manual edit 2: randomly sample 6 SNPs from LGC1.1 & LGC11.1. 
# All remaining SNPs are tightly linked on the same map position.
snp_list_invMark  <- snp_list_invMark[!(snp_list_invMark %in% 
                                          sample(row.names(KTinvLD2[["LGC11.1"]]), 
                                                 nrow(KTinvLD2[["LGC11.1"]])- 6))]
snp_list_invMark  <- snp_list_invMark[!(snp_list_invMark %in% 
                                        sample(row.names(KTinvLD2[["LGC11.1"]]), 
                                               nrow(KTinvLD2[["LGC11.1"]])- 6))]

# Manual edit 3: Remove duplicates
snp_list_invMark[duplicated(snp_list_invMark)] # Identify the duplicates

# Assign Contig74870_55710_LG10_3_1cM_ikt to LGC10.1
mean(KTinvLD2[["LGC10.1"]][,"Contig74870_55710_LG10_3_1cM_ikt"])
mean(KTinvLD2[["LGC10.2"]][,"Contig74870_55710_LG10_3_1cM_ikt"])

# Assign Contig103432_8746_LG14_10_2cM_ikt to LGC14.1/2
mean(KTinvLD2[["LGC14.1/2"]][,"Contig103432_8746_LG14_10_2cM_ikt"])
mean(KTinvLD2[["LGC14.3"]][,"Contig103432_8746_LG14_10_2cM_ikt"])

# Assign Contig103677_3202_LG14_10_2cM_ikt to LGC14.1/2
mean(KTinvLD2[["LGC14.1/2"]][,"Contig103677_3202_LG14_10_2cM_ikt"])
mean(KTinvLD2[["LGC14.3"]][,"Contig103677_3202_LG14_10_2cM_ikt"])

# Assign Contig107495_6969_LG14_10_2cM_ikt to LGC14.1/2
mean(KTinvLD2[["LGC14.1/2"]][,"Contig107495_6969_LG14_10_2cM_ikt"])
mean(KTinvLD2[["LGC14.3"]][,"Contig107495_6969_LG14_10_2cM_ikt"])

snp_list_invMark <- snp_list_invMark[!(duplicated(snp_list_invMark))]

#### 4: Filter colinear outliers ####

### There is an oversaturation of colinear outliers on LG2
### 12 'cot' markers are randomly removed
cot_snp_list <- KTsnp4$SNP[grep("_cot", KTsnp4$SNP)]
cot_snp_list <- cot_snp_list[grep("_LG2_", cot_snp_list)]

LG2cot_snp <- KTsnp4[KTsnp4$SNP %in% cot_snp_list,]

tmp <- LG2cot_snp[LG2cot_snp$MapPos %>% duplicated(),] %>% 
  group_by(MapPos) %>% slice_sample(n = 1) %>% summarise(SNP)
tmp$SNP

LG2cot_snp <- LG2cot_snp[!(LG2cot_snp$MapPos %>% duplicated()) | LG2cot_snp$SNP %in% tmp$SNP,]

keeps <- sample(LG2cot_snp$SNP, size = 10)

cot_snp_list <- KTsnp4$SNP[grep("_cot", KTsnp4$SNP)]
cot_snp_list <- cot_snp_list[grep("_LG2_", cot_snp_list, invert = TRUE)]
cot_snp_list <- c(cot_snp_list, keeps)


#### 5: Filter genetic background markers ####
Ngbg <- 384 - length(snp_list_invMark) - length(cot_snp_list)

gbg_snp_list <- sample(KTsnp4$SNP[grep("_gbg", KTsnp4$SNP)], size = Ngbg)


#### 6: Merge togther filtered datasets ####

filtered_snp_list <- c(snp_list_invMark, cot_snp_list, gbg_snp_list)

RefinedSNPpanel <- KTsnp4[KTsnp4$SNP %in% filtered_snp_list,] ### Bug alert: I am loosing 4 SNPs with this approach!!!
### This bug is caused by SNP that occur in the overlap between inversions on LG10 & LG14
# Plot across linkage groups
ggplot()+
  geom_rect(data = Invs, aes(ymin = as.numeric(gsub("LG", "", LG))-0.3, 
                             ymax = as.numeric(gsub("LG", "", LG))+0.3, 
                             xmin = Start_pos, xmax = End_pos),
            fill = "red", alpha = 0.4, col = "darkred")+
  geom_point(data = KTsnp, 
             aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             alpha = 0.2, size = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = RefinedSNPpanel, aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = c("#E06B53", "#E0D24F", "#4AE0A5", "#4AE0A5", "#552DE0"))+
  scale_y_continuous(breaks = 1:17)+
  labs(x = "Position on crab map (cM)", y = "Linkage group", colour = "Dataset")+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour = "grey", size = 0.5),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

# Save SNP IDs for refined SNP panel
write.table(RefinedSNPpanel$SNP, paste0(PATH, "SNPs_for_final_round_genotyping.txt"),
            col.names = FALSE, quote = FALSE, row.names = FALSE)

### SNP counts
RefinedSNPpanel %>% group_by(LG) %>% summarise(n())
RefinedSNPpanel %>% group_by(Source) %>% summarise(n())
RefinedSNPpanel %>% group_by(LG, Source) %>% summarise(n())

RefinedSNPpanel$INV <- 'NA'
for(i in 1:nrow(Invs)){
  LG <- Invs[i,]$LG
  Inv <- Invs[i,]$INV
  Spos <- Invs[i,]$Start_pos
  Epos <- Invs[i,]$End_pos
  
  RefinedSNPpanel[RefinedSNPpanel$LG == LG &
                    RefinedSNPpanel$MapPos >= Spos &
                    RefinedSNPpanel$MapPos <= Epos, "INV"] <- Inv
};rm(i)

tmpN <- RefinedSNPpanel %>% group_by(Source, INV) %>% summarise(n())
