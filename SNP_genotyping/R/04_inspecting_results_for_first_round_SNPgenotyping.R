####################################### Investigating results of SNP panel ####################################
### This script investigates the results of an initial SNP panel of 856 sites. These were sequenced at LGC using 
### their KASP genotyping platform. Once filtered these SNPs will be used to select 364 inversion markers

### James Reeve - University of Gothenburg
### 17/06/2022

### Preparation
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library("tidyverse")
library("ggpubr")

### Read SNP data results
PATH <- "/path/to/SNP/panel/results/"
KTsnp <- read.csv(paste0(PATH, "Genotyping-3366.015-03.csv"), skip = 14)[1:855,]

### Access information from SNP IDs
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

### Summary stats
# Unused
KTsnp$Unused <- as.integer(KTsnp$Unused)
summary(KTsnp$Unused)
ggplot(KTsnp)+geom_density(aes(x = Unused), fill = "red", colour = "red", alpha = 0.4)+theme_classic()
nrow(KTsnp[KTsnp$Unused == 0,]) # 329 SNPs were 'used' in all individuals

# Missing
KTsnp$Missing <- as.integer(KTsnp$Missing)
summary(KTsnp$Missing)
unique(KTsnp$Missing) # No SNP is missing from any individual...

# Bad
KTsnp$Bad <- as.integer(KTsnp$Bad)
summary(KTsnp$Bad)
unique(KTsnp$Bad) # No SNP from any individual was 'bad'...

# Allele frequnecies
KTsnp$q <- KTsnp$Allele.Y./100
KTsnp$p <- KTsnp$Allele.X./100

summary(KTsnp$q > 0.05) # 315 SNPs had allele Q frequency below 0.05
summary(KTsnp$p > 0.05) # 30 SNPs had allele P frequency below 0.05
summary(KTsnp$p > 0.05 & KTsnp$q > 0.05) # 345 SNPs have a frequency < 0.05 for both alleles

ggplot(KTsnp)+
  geom_density(aes(p), colour = "cyan", fill = "cyan", alpha = 0.4)+
  geom_density(aes(q), colour = "magenta", fill = "magenta", alpha = 0.4)+
  geom_vline(xintercept = 0.05, lty = 3, colour = "magenta")+
  geom_vline(xintercept = 0.95, lty = 3, colour = "cyan")+
  theme_classic()

### Plot heterozygosity
ggplot(KTsnp)+
  geom_tile(aes(x = 0, y = SNP, fill = as.integer(Het)/nrow(KTsnp)))+
  labs(y = "SNPs", fill = "Het")+
  scale_fill_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  facet_grid(~factor(LG, levels = paste0("LG", 1:17)))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank())

### Set filter out Unused samples < 5
KTsnp2 <- KTsnp[KTsnp$Unused <= 5,]

# Find any cases where SNP are sharing a contig
# Randomly keep 1 of the SNPs
KTsnp3 <- KTsnp2 %>% group_by(LG, Contig) %>%
  slice_sample(n = 1)

### Summary counts of SNPs
# Count of SNPs in each LG
# LG15 has the fewest with just 15
KTsnp %>% group_by(LG) %>% summarise(n())

# Count of SNPs per source of data
KTsnp %>% group_by(Source) %>% summarise(n())

### Counts for filtered SNPs
KTsnp3 %>% group_by(LG) %>% summarise(n())
KTsnp3 %>% group_by(Source) %>% summarise(n())

p1 <- ggplot(KTsnp)+
  geom_bar(aes(x = LG, fill = Source), position = position_dodge(preserve = "single"))+
  scale_fill_manual(values = c("#2B5266", "#A8D4EA", "#61BAE6", "#495C66", "#4B90B3"))+
  labs(x = "Linkage group", y = "Number of SNPs", title = "Before filtering")+
  scale_x_discrete(limits = paste0("LG", 1:17))+
  ylim(c(0,65))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(KTsnp2)+
  geom_bar(aes(x = LG, fill = Source), position = position_dodge(preserve = "single"))+
  scale_fill_manual(values = c("#2B5266", "#A8D4EA", "#61BAE6", "#495C66", "#4B90B3"))+
  labs(x = "Linkage group", y = "Number of SNPs", title = "After filtering")+
  scale_x_discrete(limits = paste0("LG", 1:17))+
  ylim(c(0,65))+
  theme_bw()+
  theme(axis.title.y = element_blank())

annotate_figure(ggarrange(p1, p2, common.legend = TRUE, nrow = 2, ncol = 1), left = "Number of SNPs")


############## Plot position of SNPs on linkage map #############
# Inversion positions
Invs <- read.csv("/path/to/Inversion_positions_on_old_map.csv")
Invs <- Invs[Invs$Site == "consensus",]
Invs$Start_pos <- Invs[Invs$Start == TRUE, "Pos"]
Invs$End_pos <- Invs[Invs$Start == FALSE, "Pos"]
Invs <- Invs %>% select(-Pos, -Start) %>% distinct()

# Crab map markers
LMap <- read.table("/path/to/linkage/map_v11.txt", header = TRUE)
LMap$LG <-  paste0("LG", LMap$LG)

ggplot()+
  geom_rect(data = Invs, aes(ymin = as.numeric(gsub("LG", "", LG))-0.3, 
                             ymax = as.numeric(gsub("LG", "", LG))+0.3, 
                             xmin = Start_pos, xmax = End_pos),
            fill = "red", alpha = 0.4, col = "darkred")+
  #geom_point(data = LMap, aes(x = av, y = as.numeric(gsub("LG", "", LG))), 
  #           col = "grey50", alpha = 0.1, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = KTsnp[!(KTsnp$SNP %in% KTsnp3$SNP),], 
             aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             alpha = 0.2, size = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = KTsnp3, aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             size = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = c("red", "blue", "green", "green", "gold"))+
  scale_y_continuous(breaks = 1:17)+
  labs(x = "Position on crab map (cM)", y = "Linkage group", colour = "Dataset")+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour = "grey", size = 0.5),
        panel.grid.major.y = element_blank())


############## Filter based on correlation of SNPs to inversions #################
# Access SNP genotypes
KTsnpGT <- read.csv(paste0(PATH, "Genotyping-3366.015-03 Grid.csv"), skip = 7)[-c(1),] # Note 1st row is control marker
KTsnpGT <- KTsnpGT[,KTsnp3$SNP] # Filter based of summary stats

# Change exact genotypes to 0/0, 0/1 & 1/1
for(i in 1:ncol(KTsnpGT)){
  snp <- KTsnpGT[,i]
  
  Ref <- KTsnp3$Allele.X[i]
  Alt <- KTsnp3$Allele.Y[i]
  
  KTsnpGT[snp == paste(Ref, Ref, sep = ":"), i] <- 0
  KTsnpGT[snp == paste(Alt, Alt, sep = ":"), i] <- 2
  KTsnpGT[snp == paste(Ref, Alt, sep = ":") | 
            snp == paste(Alt, Ref, sep = ":"), i] <- 1
}

# Change "?", "Uncallable" & "Bad" to NA
KTsnpGT[KTsnpGT == "?" | KTsnpGT == "Uncallable" | KTsnpGT == "Bad"] <- NA

# Convert genotypes to numeric
KTsnpGT <- as.data.frame(sapply(KTsnpGT, as.numeric))

### Function to make LD matrix between SNPs in specific inversions
inversion.LD <- function(genotype.data, inversion, SNP.type){
  # Get inversion details
  LG <- Invs[Invs$INV == inversion, "LG"]
  Spos <- Invs[Invs$INV == inversion, "Start_pos"]
  Epos <- Invs[Invs$INV == inversion, "End_pos"]
  
  # Find SNPs in the inversion
  snp_list <- unlist(KTsnp3[KTsnp3$LG == LG &
                       KTsnp3$MapPos >= Spos &
                       KTsnp3$MapPos <= Epos, "SNP"])
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

### Loop through inversions for parallelogram SNPs
TEST <- sapply(Invs$INV, function(inv){
  try(inversion.LD(KTsnpGT, inv, "ipr|ikt"))
})

# Remove any correlation matrices that failed
TEST <- TEST[lapply(TEST, class) != "try-error"]

annotate_figure(ggarrange(plotlist = lapply(names(TEST), function(x){
  ggplot(reshape::melt(TEST[[x]]), aes(x = abs(value)))+
    geom_histogram(binwidth = 0.05)+
    geom_density(aes(y = 0.05*..count..), 
                 colour = "purple", fill = "purple", alpha = 0.4)+
    labs(title = paste0(x,": (Nsnp =", nrow(TEST[[x]]), ")"))+
    xlim(c(0,1))+
    theme_classic()+
    theme(axis.title = element_blank())
})), left = "SNP count", bottom = "Linkage disequilibirum (|r|)")

ggarrange(plotlist = lapply(names(TEST), function(x){
  ggplot(reshape::melt(TEST[[x]]))+
    geom_tile(aes(x = X1, y = X2, fill = value))+
    scale_fill_gradient2(low = "#8856a7", high = "#8856a7", 
                         mid = "#e0ecf4", na.value = "white")+
    ggtitle(x)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 2))
}), common.legend = TRUE, legend = "right")

###############################################
IKTinv <- sapply(Invs$INV, function(inv){
  try(inversion.LD(KTsnpGT, inv, "ikt"))
})

# Remove any correlation matrices that failed
IKTinv <- IKTinv[lapply(IKTinv, class) != "try-error"]

annotate_figure(ggarrange(plotlist = lapply(names(IKTinv), function(x){
  ggplot(reshape::melt(IKTinv[[x]]), aes(x = abs(value)))+
    geom_histogram(binwidth = 0.05)+
    geom_density(aes(y = 0.05*..count..), 
                 colour = "purple", fill = "purple", alpha = 0.4)+
    labs(title = paste0(x,": (Nsnp =", nrow(IKTinv[[x]]), ")"))+
    xlim(c(0,1))+
    theme_classic()+
    theme(axis.title = element_blank())
})), left = "SNP count", bottom = "Linkage disequilibirum (|r|)")

ggarrange(plotlist = lapply(names(IKTinv), function(x){
  ggplot(reshape::melt(IKTinv[[x]]))+
    geom_tile(aes(x = X1, y = X2, fill = value))+
    scale_fill_gradient2(low = "#8856a7", high = "#8856a7", 
                         mid = "#e0ecf4", na.value = "white")+
    ggtitle(x)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 2))
}), common.legend = TRUE, legend = "right")


#### Setting some kind of threshold
### Find N snps with tightest LD
### where N is the number to retain based on Roger's points
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

TEST2 <- lapply(TEST, filter.by.LD, max.num.snps = 5)

ggarrange(plotlist = lapply(names(TEST2), function(x){
  ggplot(reshape::melt(TEST2[[x]]))+
    geom_tile(aes(x = X1, y = X2, fill = value))+
    scale_fill_gradient2(low = "#8856a7", high = "#8856a7", 
                         mid = "#e0ecf4", na.value = "white")+
    ggtitle(x)+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 2))
}), common.legend = TRUE, legend = "right")

### Return best SNPs for each inversion
BestSNPs <- unlist(lapply(TEST2, function(x){rownames(x)}))
# Manual edit: remove some SNPs which still had poor correlations
BestSNPs <- BestSNPs[!(BestSNPs %in% c("Contig41294_55142_LG14_0_4cM_ikt", 
                                       "Contig_43489_81919_LG10_14_8cM_ikt", 
                                       "Contig39483_65696_LG12_48_7cM_ikt"))]
# Manual edit 2: randomly sample 5 SNPs from LGC11.1. All SNPs are tightly linked on the same map position.
# Randomly select 4 to remove
BestSNPs <- BestSNPs[!(BestSNPs %in% sample(BestSNPs[grep("LG11", BestSNPs)], 4))]

invMarks <- KTsnp3[KTsnp3$SNP %in% BestSNPs,]

ggplot()+
  geom_rect(data = Invs, aes(ymin = as.numeric(gsub("LG", "", LG))-0.3, 
                             ymax = as.numeric(gsub("LG", "", LG))+0.3, 
                             xmin = Start_pos, xmax = End_pos),
            fill = "red", alpha = 0.4, col = "darkred")+
  #geom_point(data = LMap, aes(x = av, y = as.numeric(gsub("LG", "", LG))), 
  #           col = "grey50", alpha = 0.1, position = position_jitter(width = 0.2, height = 0.2))+
  #geom_point(data = KTsnp[!(KTsnp$SNP %in% KTsnp3$SNP),], 
  #           aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
  #           alpha = 0.2, size = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = invMarks, aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = c("red", "blue", "green", "green", "gold"))+
  scale_y_continuous(breaks = 1:17)+
  labs(x = "Position on crab map (cM)", y = "Linkage group", colour = "Dataset")+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour = "grey", size = 0.5),
        panel.grid.major.y = element_blank())


####### Filter out 12 colinear outliers on LG2 #########
cot_snp_list <- KTsnp3$SNP[grep("_cot", KTsnp3$SNP)]
cot_snp_list <- cot_snp_list[grep("_LG2_", cot_snp_list)]

LG2cot_snp <- KTsnp3[KTsnp3$SNP %in% cot_snp_list,]

tmp <- LG2cot_snp[LG2cot_snp$MapPos %>% duplicated(),] %>% 
  group_by(MapPos) %>% slice_sample(n = 1) %>% summarise(SNP)
tmp$SNP

LG2cot_snp <- LG2cot_snp[!(LG2cot_snp$MapPos %>% duplicated()) | LG2cot_snp$SNP %in% tmp$SNP,]

keeps <- sample(LG2cot_snp$SNP, size = 10)

cot_snp_list <- KTsnp3$SNP[grep("_cot", KTsnp3$SNP)]
cot_snp_list <- cot_snp_list[grep("_LG2_", cot_snp_list, invert = TRUE)]
cot_snp_list <- c(cot_snp_list, keeps)

cotMarks <- KTsnp3[KTsnp3$SNP %in% cot_snp_list,]

####### Randomly sample genetic background SNPs
Ngbg <- 384 - nrow(cotMarks) - nrow(invMarks)

gbg_snp_list <- sample(KTsnp3$SNP[grep("_gbg", KTsnp3$SNP)], size = Ngbg)

gbgMarks <- KTsnp3[KTsnp3$SNP %in% gbg_snp_list,]

###### Merge filtered datasets for final SNP panel #####
RefinedSNPpanel <- rbind.data.frame(gbgMarks, cotMarks, invMarks)

ggplot()+
  geom_rect(data = Invs, aes(ymin = as.numeric(gsub("LG", "", LG))-0.3, 
                             ymax = as.numeric(gsub("LG", "", LG))+0.3, 
                             xmin = Start_pos, xmax = End_pos),
            fill = "red", alpha = 0.4, col = "darkred")+
  #geom_point(data = LMap, aes(x = av, y = as.numeric(gsub("LG", "", LG))), 
  #           col = "grey50", alpha = 0.1, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = KTsnp, 
             aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             alpha = 0.2, size = 0.4, position = position_jitter(width = 0.2, height = 0.2))+
  geom_point(data = RefinedSNPpanel, aes(x = MapPos, y = as.numeric(gsub("LG", "", LG)), colour = Source),
             position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = c("red", "blue", "green", "green", "gold"))+
  scale_y_continuous(breaks = 1:17)+
  labs(x = "Position on crab map (cM)", y = "Linkage group", colour = "Dataset")+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour = "grey", size = 0.5),
        panel.grid.major.y = element_blank())
