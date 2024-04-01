################# Plot PCAs coloured by heterozygosity ###############################
### This script combines PCA scores of genotypes within inversion's breakpoints with
### the average heterozygosity for these same genotypes. The original PCA script was
### developed by Eva Koch using an adgenet pipeline. 

### James Reeve - University of Gothenburg
### 08/12/2021

### Preparation
rm(list = ls())
dev.off()
setwd("~")
PATH <- "working/dir/path/"

### Packages
library(tidyverse)
library(ggpubr)

### A: Function to call heterozygosity across inversion ####
inversion.heterozygosity <- function(inversion){
  ### Access data
  # Genotypes:
  LG <- gsub("[.].*","",gsub("C","",inversion)) # Extract LG from inversion name
  GT <- read.table(paste0(PATH, "genotypes/KT_", LG, "_genotypes_BGI.txt")) # Access genotypes
  
  ### Data wrangling
  # Filter to diagnostic SNPs
  if(grepl("/", inversion)){
    dSNPs <- scan(paste0(PATH, "Evas_diagnostic_SNPs/", gsub("/", ".", inversion), "DiagSNPs.txt"), character())
  } else {
    dSNPs <- scan(paste0(PATH, "Evas_diagnostic_SNPs/", inversion, "DiagSNPs.txt"), character())
  }
  dGT <- GT[which(rownames(GT) %in% dSNPs),]
  
  # Create a list of sample IDs
  SnailIDs <- colnames(dGT)
  
  #Get SNP & contig names
  dGT$SNP <- rownames(dGT) # Add SNP names
  snps <- strsplit(dGT$SNP, "_") # Temporary vector storing SNP details
  dGT$contig <- sapply(snps, "[[", 1) # Add contigs names
  dGT$pos <- as.numeric(sapply(snps, "[[", 2)) # Add contig positions
  
  ### Apply loop to sum up the number of genotypes for each individual
  Het <- lapply(SnailIDs, function(i){
    dGT %>% select(contig, pos, i) %>%     # Select single individual
      rename("ind" = i) %>%               # Rename the name of the last column to "ind"
      na.omit() %>%                       # Subset data to specific contig and remove NAs
      summarise(Inv = inversion,               # Genotype counts
                ID = i,
                Nsnp = n(),
                AA = sum(ind == 2),
                Aa = sum(ind == 1),
                aa = sum(ind == 0)) %>%
      mutate(p = (2*aa + Aa) / (2*Nsnp))     # Get the frequency of the ref allele
  })
  # Transform to dataframe
  Het <- do.call(rbind.data.frame, Het)
  
  ### Return the genotype counts
  return(Het)
}



### B: Calculate heterozygosity of each inversion ####
# Access inversion breakpoints
InvPos <- read.csv(paste0(PATH, "Inversion_positions_on_old_map.csv"), header = T)

# Get vector of inversion names
Invs <- unique(InvPos$INV)

# Run the function for all inversions in an apply loop
He <- lapply(Invs, inversion.heterozygosity)
He <- do.call(rbind.data.frame, He)

### Save heterozygosity of each inversion
write.table(He, paste0(PATH, "heterozygosity/Inversion_Heterozygosity_DiagSNPs.txt"), 
            row.names = FALSE, quote = FALSE)



### C: Merge with PCA data ####
### Generate plot for each inversion in an apply loop
Plots <- lapply(Invs, function(INV){
  # Access PC1 & PC2 for given inversion
  if(grepl("/", INV)){
    PCA <- read.table(paste0(PATH, "Evas_diagnostic_SNPs/", gsub("/", ".", INV), "PCscores.txt"))[,1:2]
  } else {
    PCA <- read.table(paste0(PATH, "Evas_diagnostic_SNPs/", INV, "PCscores.txt"))[,1:2]
  }
  PCA$ID <- rownames(PCA) # Add column for snail IDs
  
  # Access heterozygosity of current inversion
  He.2 <- He[He$Inv == INV,]
  
  # Join heterozygosity with PCA
  He.PCA <- He.2 %>%
    group_by(ID) %>%
    summarise("Het" = sum(Aa) / sum(Nsnp)) %>%      # Calculate frequency of heterozygotes
    left_join(PCA, by = "ID") %>%                   # Join files
    drop_na()

  # Plot the results
  Res <- ggplot(He.PCA, aes(x = PC1, y = PC2))+
    geom_point(aes(colour = Het))+
    labs(colour = expression(H[E]), title = INV)+
    scale_colour_gradientn(colours = terrain.colors(40))+
    theme_bw()
  
  # Save plots
  ggsave(paste0(PATH, "heterozygosity/",
                if(grepl("/", INV)){gsub("/", ".", INV)}else{INV}, "_PCAplot.tiff"), 
                plot = Res, device = "tiff", units = "px", width = 2000, height = 1200)
  
  #### Alternative plot Het vs PC1 #####
  #Res <- ggplot(He.PCA, aes(x = Het, y = PC2))+
  #  geom_point()+
  #  labs(x = expression(H[E]), title = INV)+
  #  theme_bw()
  #####
  return(Res)
})

### Multipanel plot
ggarrange(plotlist = Plots, common.legend = TRUE, legend = "bottom")

