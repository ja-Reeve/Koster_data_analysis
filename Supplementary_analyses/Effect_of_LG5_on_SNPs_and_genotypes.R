##################### Analysing the effect of LG5 SNPs on results #######################
### SNPs on LG5 were included in the genetic background. However, there is a large region
### of tight linkage in the middle of LG5, sometimes referred to an inversion LGC5.1. This
### script checked the impact of SNPs on LG5 and those within LGC5.1, on the genetic 
### background PCA and checked if they are associated with phenotypic changes. Ultimately,
### there was an association, so we removed LG5 SNPs from the analysis. 

### James Reeve - DTU Aqua
### 2026-06-25

### Preparation
rm(list = ls())
try(dev.off())
options(stringsAsFactors = FALSE)

### Packages
library(ggpubr)
library(adegenet)
library(brms)

### Filepath
PATH <- "file/path/"


#### 1. Access data ####

### A: Genotypes
GT <- read.csv(paste0(PATH, "Genotypes/Koster_SNP_imputed_20221031.csv"), header = TRUE)

### B: SNP information
SNP_info <- read.csv(paste0(PATH, "Genotypes/Koster_SNP_info_FINAL.csv"), header = TRUE)

### C: Phenotype PCA scores
Pheno <- read.csv("Documents/snail-scape/PCAmix/principal_components/Pheno_PC1-8.V3.csv")
colnames(Pheno)[1] <- "snail"

### D: Dissection records - only needed from brms model
diss <- read.csv(paste0(PATH, "Phenotypes/KT_dissection_log.csv"))

### E: Get shellLenth from shellShaper parameters
ShS <- read.csv(paste0(PATH, "Phenotypes/KT_ShellShaper_parameters.v2.csv"))
ShS$snailID <- substr(ShS$snailID, 1, 6)
# Merge with Pheno
Pheno <- merge(Pheno, ShS[,c("snailID", "shellLength")], by.x = "snail", by.y = "snailID")



#### 2. Data wrangeling ####

### A: Filter genotypes
# Filter to just Koster-Tjärnö sampling
GT <- GT[grep("KT", GT$snail), ]
# Note: 3 snails (KT088b, KT128a & KT172f) were included in both the WGS and SNP genotyping
# I am filtering out the 2nd entry in GT, as this is the duplicate from the WGS data.
GT <- GT[!duplicated(GT$snail),]

# Move snail ID to be a rowname
rownames(GT) <- GT$snail
GT$snail <- NULL

# Creat new dataset with inversion genotypes
GTinv <- GT[,grep("LGC", colnames(GT))]

# Drop inversion genotypes
GT <- GT[,grep("Contig", colnames(GT))]


### B: Set "Contig5977_84235" as LGC14.1/2
### This inversion was tricky to genotype as most markers were tightly linked and located
### in the overlapping boundaries with LGC14.3. 
GTinv$LGC14.1 <- GT$Contig5977_84235
GTinv$LGC14.1[GTinv$LGC14.1 == 0] <- "RR"
GTinv$LGC14.1[GTinv$LGC14.1 == 1] <- "RA"
GTinv$LGC14.1[GTinv$LGC14.1 == 2] <- "AA"


### C: Extract LG5 genotypes
# Create SNP id in metadata
SNP_info$SNP <- paste(SNP_info$CHROM, SNP_info$POS, sep = "_")

# List of SNPs on LG5
LG5snp <- SNP_info[SNP_info$LG == 5, "SNP"]

# Genotypes on LG5
GT.LG5 <- GT[,LG5snp]


### D: Simplify the phenotypes
# Rescore juveniles to females
diss[diss$sex == "J", "sex"] <- "F"

# Rescore as c("adult", "juvenile")
diss$adult <- ifelse(diss$adult, "adult", "juvenile")



#### 3. Plot ####

### A: Function to plot Phenotype PCA coloured by SNP genotypes
GT.PC.plot <- function(SNP){
  # Merge SNPs with phenotype PCs
  dat <- Pheno
  dat$GT <- factor(GT.LG5[Pheno$snail, SNP])
  
  # Scatterplot
  ggplot(dat)+
    geom_point(aes(x = PC1, y = PC2, colour = GT), size = 0.4, alpha = 0.9)+
    labs(x = "PC1", y = "PC2", colour = "GT", title = SNP)+
    scale_colour_manual(values = c("#a6cee3", "orange", "#b2df8a"), na.value = "grey70")+
    scale_x_continuous(breaks = c(0), expand = FALSE)+
    scale_y_continuous(breaks = c(0), expand = FALSE)+
    theme_bw()+
    theme(plot.title = element_text(size = 8),
          panel.border = element_rect(colour = "grey80"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          margins = margin(0, 0.1, 1, 0.1))
}


### B: Loop through plots
annotate_figure(
  ggarrange(plotlist = lapply(LG5snp, GT.PC.plot), common.legend = TRUE, legend = "right"),
  left = "PC2", bottom = "PC1"
)



#### 4. Genotype linkage group ####
### Genotyping is a two step process. First all SNPs are put into a principal component
### analysis. Clusters will then be called using K-means. In some cases a derived arrangement
### (likely the inverted arrangement) will be rare. Both K = 2 and 3 will be run to check.

### Two different genotyping sets will be tried. First all SNPs will be combined. Second, only
### SNPs in the area containing LGC5.1 (15 - 47 cM) will be included. 

### A: PCA for all SNPs on LG5
# Convert genotypes to genind format
gen.LG5 <- df2genind(GT.LG5, ploidy = 2, type = "codom",
                     sep = "", pop = rownames(GT.LG5))
# Reformat as genpop
gp.LG5 <- genind2genpop(gen.LG5)

# Centre the scores
gp.LG5 <- scaleGen(gp.LG5, scale = FALSE)

# Run PCA
pca.LG5 <- dudi.pca(gp.LG5, scale = FALSE, center = FALSE, scannf = FALSE, nf = 2)
# Remove intermediate files
rm(gen.LG5, gp.LG5)


### B: Plot PCA results
LG5 <- as.data.frame(pca.LG5$li)
LG5$Snail <- rownames(pca.LG5$li)

# Variance explained (%)
Vexp.PC1 <- round(pca.LG5$eig[1] / sum(pca.LG5$eig) * 100, digits = 2)
Vexp.PC2 <- round(pca.LG5$eig[2] / sum(pca.LG5$eig) * 100, digits = 2)

# Scatter plot
ggplot(LG5)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = Axis1, y = Axis2),
             alpha = 0.4, colour = "dodgerblue")+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"))+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))


### C:K-means clustering
K2.LG5 <- kmeans(LG5$Axis1, centers = 2)
K3.LG5 <- kmeans(LG5$Axis1, centers = 3)

# Add clusters to PCA
LG5$K2 <- K2.LG5$cluster
LG5$K3 <- K3.LG5$cluster

# Plot with clusters
p.K2 <- ggplot(LG5)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = Axis1, y = Axis2, colour = factor(K2)), alpha = 0.4)+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"),
       colour = "Clusters")+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

p.K3 <- ggplot(LG5)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = Axis1, y = Axis2, colour = factor(K3)), alpha = 0.4)+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"),
       colour = "Clusters")+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

ggarrange(p.K2, p.K3)
# K3 has better fit, but it doesn't look like an inversion.


### D: Rename arrangements based on the most common genotype

# Sum up genotype counts
N.gt.LG5 <- data.frame("Snail" = LG5$Snail,
                       "Cluster" = LG5$K3,
                       "Nref" = apply(GT.LG5, 1, function(X){sum(X == 0)}),
                       "Nhet" = apply(GT.LG5, 1, function(X){sum(X == 1)}),
                       "Nalt" = apply(GT.LG5, 1, function(X){sum(X == 2)})) 


# Summary of genotype counts per cluster
N.gt.LG5 %>% group_by(Cluster) %>%
  summarise("P.ref" = sum(Nref) / (sum(Nref) + sum(Nhet) + sum(Nalt)),
            "P.het" = sum(Nhet) / (sum(Nref) + sum(Nhet) + sum(Nalt)),
            "P.alt" = sum(Nalt) / (sum(Nref) + sum(Nhet) + sum(Nalt)))



# Cenotype column based on the most cluster with the most of one type of genotype
LG5$GT <- "RR"
LG5$GT[LG5$cluster == 1] <- "RA"
LG5$GT[LG5$cluster == 2] <- "AA"


### E: Repeat for LGC5.1

# List of contigs between LGC5.1's boundaries
LGC5.1snps <- SNP_info[SNP_info$LG == 5 &
                         SNP_info$mp >= 15 &
                         SNP_info$mp < 47, "SNP"]

# Subset genotypes
GT.LGC5.1 <- GT.LG5[,LGC5.1snps]

# PCA
gen.LGC5.1 <- df2genind(GT.LGC5.1, ploidy = 2, type = "codom",
                     sep = "", pop = rownames(GT.LGC5.1))
gp.LGC5.1 <- genind2genpop(gen.LGC5.1)
gp.LGC5.1  <- scaleGen(gp.LGC5.1 , scale = FALSE)
pca.LGC5.1 <- dudi.pca(gp.LGC5.1 , scale = FALSE, center = FALSE, scannf = FALSE, nf = 4)
rm(gen.LGC5.1, gp.LGC5.1 )

# Create data frame
LGC5.1 <- as.data.frame(pca.LGC5.1$li)
LGC5.1$Snail <- rownames(pca.LGC5.1$li)

Vexp.PC1 <- round(pca.LGC5.1$eig[1] / sum(pca.LGC5.1$eig) * 100, digits = 2)
Vexp.PC2 <- round(pca.LGC5.1$eig[2] / sum(pca.LGC5.1$eig) * 100, digits = 2)
Vexp.PC3 <- round(pca.LGC5.1$eig[3] / sum(pca.LGC5.1$eig) * 100, digits = 2)
Vexp.PC4 <- round(pca.LGC5.1$eig[4] / sum(pca.LGC5.1$eig) * 100, digits = 2)

# Scatter plot
ggplot(LGC5.1)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = Axis1, y = Axis2),
             alpha = 0.4, colour = "dodgerblue")+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"))+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

# K-means clustering
K2.LGC5.1 <- kmeans(LGC5.1[,"Axis2"], centers = 2)
K3.LGC5.1 <- kmeans(LGC5.1[,"Axis2"], centers = 3)

# Add clusters to PCA
LGC5.1$K2 <- K2.LGC5.1$cluster
LGC5.1$K3 <- K3.LGC5.1$cluster

# Plot with clusters
p.K2 <- ggplot(LGC5.1)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = Axis1, y = Axis2, colour = factor(K2)), alpha = 0.4)+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"),
       colour = "Clusters")+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

p.K3 <- ggplot(LGC5.1)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = Axis1, y = Axis2, colour = factor(K3)), alpha = 0.4)+
  labs(x = paste0("PC1 (", Vexp.PC1, "%)"),
       y = paste0("PC2 (", Vexp.PC2, "%)"),
       colour = "Clusters")+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

ggarrange(p.K2, p.K3)
# Neither is particularly good at clustering points
# I've looked at PC3 & 4, and it's pretty similar


TEST <- prcomp(GT.LGC5.1, center = TRUE, scale. = FALSE)

LGC5.1 <- as.data.frame(TEST$x)
LGC5.1$Snail <- rownames(LGC5.1)

ggplot(LGC5.1)+
  geom_point(aes(x = PC1, y = PC2))


K2.LGC5.1 <- kmeans(LGC5.1$PC1, centers = 2)
K3.LGC5.1 <- kmeans(LGC5.1$PC1, centers = 3)

# Add clusters to PCA
LGC5.1$K2 <- K2.LGC5.1$cluster
LGC5.1$K3 <- K3.LGC5.1$cluster

# Plot with clusters
p.K2 <- ggplot(LGC5.1)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = PC1, y = PC2, colour = factor(K2)), alpha = 0.4)+
  labs(x = paste0("PC1"),
       y = paste0("PC2"),
       colour = "Clusters")+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

p.K3 <- ggplot(LGC5.1)+
  geom_vline(xintercept = 0, colour = "grey85")+
  geom_hline(yintercept = 0, colour = "grey85")+
  geom_point(aes(x = PC1, y = PC2, colour = factor(K3)), alpha = 0.4)+
  labs(x = paste0("PC1"),
       y = paste0("PC2"),
       colour = "Clusters")+
  theme_void()+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90))

ggarrange(p.K2, p.K3)

####### NOTE: adgenet and prcomp are giving noticably different results! ############

# Genotype
N.gt.LGC5.1 <- data.frame("Snail" = LGC5.1$Snail,
                       "Cluster" = LGC5.1$K3,
                       "Nref" = apply(GT.LGC5.1, 1, function(X){sum(X == 0)}),
                       "Nhet" = apply(GT.LGC5.1, 1, function(X){sum(X == 1)}),
                       "Nalt" = apply(GT.LGC5.1, 1, function(X){sum(X == 2)})) 
N.gt.LGC5.1 %>% group_by(Cluster) %>%
  summarise("P.ref" = sum(Nref) / (sum(Nref) + sum(Nhet) + sum(Nalt)),
            "P.het" = sum(Nhet) / (sum(Nref) + sum(Nhet) + sum(Nalt)),
            "P.alt" = sum(Nalt) / (sum(Nref) + sum(Nhet) + sum(Nalt)))
LGC5.1$GT <- "RR"
LGC5.1$GT[LGC5.1$K3 == 2] <- "RA"
LGC5.1$GT[LGC5.1$K3 == 3] <- "AA"



#### 5. Compare to phenotypes ####

### Merge LGC5.1 genotypes with phenotype PCs
Pheno.LGC5.1 <- merge(LGC5.1[,c("Snail", "GT")], Pheno,
                      by.x = "Snail", by.y = "snail")

### Means per phenotype
Pheno.means <- Pheno.LGC5.1 %>% 
  pivot_longer(cols = c(PC1, PC2, shellLength), values_to = "score", names_to = "PC") %>%
  group_by(GT, PC) %>%
  summarise("mean" = mean(score))

### Violin plots
p1 <- ggplot(Pheno.LGC5.1, aes(x = GT, y = PC1))+
  geom_violin()+
  geom_jitter(width = 0.3, height = 0, size = 0.8, alpha = 0.4)+
  geom_point(data = Pheno.means[Pheno.means$PC == "PC1",], aes(x = GT, y = mean),
             size = 3, pch = 15)+
  labs(x = "LGC5.1")+
  theme_classic()

p2 <- ggplot(Pheno.LGC5.1, aes(x = GT, y = PC2))+
  geom_violin()+
  geom_jitter(width = 0.3, height = 0, size = 0.8, alpha = 0.4)+
  geom_point(data = Pheno.means[Pheno.means$PC == "PC2",], aes(x = GT, y = mean),
             size = 3, pch = 15)+
  labs(x = "LGC5.1")+
  theme_classic()

p3 <- ggplot(Pheno.LGC5.1, aes(x = GT, y = shellLength))+
  geom_violin()+
  geom_jitter(width = 0.3, height = 0, size = 0.8, alpha = 0.4)+
  geom_point(data = Pheno.means[Pheno.means$PC == "shellLength",], aes(x = GT, y = mean),
             size = 3, pch = 15)+
  labs(x = "LGC5.1")+
  theme_classic()

ggarrange(p1, p2, p3, nrow = 1)



#### 6. Test of brms fit without LGC5.1 ####

### A: genetic covariance matrix from SNP genotypes
# Filter to out LGC5.1 from genotypes
GT_snp <- GT[, !(colnames(GT) %in% LGC5.1snps)]

# Filter to those with genotypes
GT_snp <- GT_snp[rownames(GT_snp) %in% Pheno$snail,]

# Transpose and divide by 2
GT_snp <- t(GT_snp / 2)

# Covariance matrix
GT_cov <- cov(GT_snp)

# Reset diagonals as 1
# Required to avoid a 'within-group covariance' error in brms
diag(GT_cov) <- diag(GT_cov) + min(GT_cov)/100


### B: Create joint dataset
GTinv$snail <- rownames(GTinv)
dat <- merge(Pheno, GTinv, by = "snail")


### C: Add on sex and 'adult' from dissection log
dat <- merge(dat, diss[,c("snail_ID", "sex", "adult")],
             by.x = "snail", by.y = "snail_ID")


### D: Run brms model
mod_sansLGC5.1 <- brm(PC1 ~ sex + adult +
                        LGC1.1 + LGC1.2 + LGC2.1 + LGC4.1 + LGC6.1.2 + 
                        LGC7.1 + LGC7.2 + LGC9.1 + LGC10.1 + LGC10.2 +
                        LGC11.1 + LGC12.1 + LGC12.2 + LGC12.3 + LGC12.4 +
                        LGC14.1 + LGC14.3 + LGC17.1 + (1|gr(snail, cov = GD)),
               data = dat, data2 = list(GD = GT_cov),
               iter = 1000, chains = 4)
mod_sansLGC5.1 <- add_criterion(mod_sansLGC5.1, criterion = "loo")

# Null model
mod_sansLGC5.1_null <- brm(PC1 ~ sex + adult + (1|gr(snail, cov = GD)),
                           data = dat, data2 = list(GD = GT_cov),
                           iter = 1000, chains = 4)
mod_sansLGC5.1_null <- add_criterion(mod_sansLGC5.1_null, criterion = "loo")



#### 7. Plot PCA with and without LGC5.1 SNPs ####

### Filter genotypes to genetic background
GTgbg <- GT[, SNP_info[SNP_info$inv == "background", "SNP"]]

### All SNPs
pca.ALL <- prcomp(GTgbg, center = TRUE, scale. = FALSE)

PC_all <- as.data.frame(pca.ALL$x)
PC_all$Snail <- rownames(PC_all)

PC_all <- merge(PC_all, LGC5.1[,c("Snail", "GT")], by = "Snail")

Vexp_all <- summary(pca.ALL)$importance

p.allSNPs <- ggplot(PC_all)+
  geom_point(aes(x = PC1, y = PC2, colour = GT))+
  scale_colour_manual(values = c("#a6cee3", "orange", "#b2df8a"))+
  labs(x = paste0("PC1 (", round(Vexp_all[2,1]*100, 2), "%)"), 
       y = paste0("PC1 (", round(Vexp_all[2,2]*100, 2), "%)"),
       title = "All SNPs")+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  theme_bw()+
  theme(axis.text = element_blank())

### Sans LGC5.1 SNPs
# Filter out LGC5.1
GT_sansLGC5.1 <- GTgbg[, !(colnames(GTgbg) %in% LGC5.1snps)]
# Run PCA
pca.sansLGC5.1 <- prcomp(GT_sansLGC5.1, center = TRUE, scale. = FALSE)

PC_sansLGC5.1 <- as.data.frame(pca.sansLGC5.1$x)
PC_sansLGC5.1$Snail <- rownames(PC_sansLGC5.1)

PC_sansLGC5.1 <- merge(PC_sansLGC5.1, LGC5.1[,c("Snail", "GT")], by = "Snail")

Vexp_sansLGC5.1 <- summary(pca.sansLGC5.1)$importance

p.sansLGC5.1 <- ggplot(PC_sansLGC5.1)+
  geom_point(aes(x = PC1, y = PC2, colour = GT))+
  scale_colour_manual(values = c("#a6cee3", "orange", "#b2df8a"))+
  labs(x = paste0("PC1 (", round(Vexp_sansLGC5.1[2,1]*100, 2), "%)"), 
       y = paste0("PC1 (", round(Vexp_sansLGC5.1[2,2]*100, 2), "%)"),
       title = "Excluding LGC5.1")+
  scale_x_continuous(breaks = 0)+
  scale_y_continuous(breaks = 0)+
  theme_bw()+
  theme(axis.text = element_blank())


ggarrange(p.allSNPs, p.sansLGC5.1, common.legend = TRUE)
