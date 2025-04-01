#=============================================================================#
# Research Project 1
# Section 1: Differentially expressed miRNAs analysis and visualization
#																	                                       		   															                                #
# Date: January 24, 2025											                                
# Author: Mai Le, ID: i6375777
# Maastricht University                                                       #
#=============================================================================#

setwd("D:/UM/Exp")
library(GEOquery)
library(affy)
library(tidyverse)

#-----------------------------------------------------------------------------#
# Section 1.a: Differentially expressed miRNAs analysis and visualization 
# between HCM and NF
#-----------------------------------------------------------------------------#
# Part 1: Import all necessary files from GEO database

#Access and import GSE36946 dataset for HCM
gse_accession <- "GSE36946"
eset <- getGEO(gse_accession, GSEMatrix = TRUE, AnnotGPL = TRUE);
View(eset)

#extract miRNA expression dataset
miR_data <- exprs(eset[[1]]);
View(miR_data)
miR_data <- as.data.frame(miR_data);

#extract metadata file
sample_metadata <- pData(eset[[1]]);
colnames(sample_metadata)[c(34, 36, 37)] <- c("age", "sample_type", "sex");
write.table(sample_metadata, file = "metadata_miR_HCM.txt", sep = "\t", row.names = F, quote = F)

#-----------------------------------------------------------------------------#
# Part 2: Quality control

# 2.1 Quality control - Box plots
if (!requireNamespace("tidyr", quietly = TRUE)) BiocManager::install("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library(tidyr)
library(ggplot2)

#Box plots of disease group
all(rownames(sample_metadata) == colnames(miR_data));
plot_HCM <- gather(miR_data[, sample_metadata$sample_type == "case"], 
                      key = "SampleID", value = "log2CPM")
HCM_distr <- ggplot(plot_HCM, aes(x = SampleID, y = log2CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(miR_data[, sample_metadata$sample_type == "case"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution in the HCM group") # add title
HCM_distr

#Box plots of control group
plot_NF <- gather(miR_data[, sample_metadata$sample_type == "control"], 
                   key = "SampleID", value = "log2CPM")
HCM_distr <- ggplot(plot_NF, aes(x = SampleID, y = log2CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(miR_data[, sample_metadata$sample_type == "case"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution of miRNA in the NF group") # add title
HCM_distr

# 2.2 Quality control - density plot

par(mar=c(1,1,1,1));
if (!requireNamespace("RColorBrewer", quietly = TRUE)) BiocManager::install("RColorBrewer")
library(RColorBrewer)

HCM_dens <- miR_data[, sample_metadata$sample_type == "case"];
nsamples <- ncol(HCM_dens);
# Create a color palette for all samples in this group
col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(nsamples);
# Density plots of HCM distribution
HCM_densityplot <- plot(density(HCM_dens[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.26), 
                las = 2, main = "", xlab = "")
title(main="Density plots of log2CPMs distribution in the HCM group ", xlab="Log2cpm")
for (i in 2:nsamples){
  den <- density(HCM_dens[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

#-----------------------------------------------------------------------------#
#Part 3: Annotation

# extract annotation information
sample_anno <- featureData(eset[[1]]);
View(sample_anno)
anno <- data(sample_anno[[2]]);
m <- sample_anno[[7]];  #sample_anno dataset = data dataset
View(m)
n<- sample_anno[[1]];
View(n)
anno <- data.frame(m,n);
View(anno)
colnames(anno) <- c ("miR_ID","PROBEID");
rownames(anno) <-anno$PROBEID;

#merge miRNA expression data with annotation
miR_HCM = cbind(miR_data,anno);
write.table(miR_HCM, file = "miR_HCM.txt", sep = "\t", row.names = F, quote = F)

#-----------------------------------------------------------------------------#
# Part 4: Differentially expressed miRNAa analysis 

#load required packages
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(edgeR)
library(limma)

# create a design matrix
design_2 <- model.matrix(~ 0 + sample_type , data = sample_metadata);
# Fit a linear model to the expression data for each gene 
fit_2 <- lmFit(miR_data, design_2);
# Create a contrast matrix for comparing expression levels between different conditions 
contrast_matrix_2<- makeContrasts(
  HCMvsNF = sample_typecase - sample_typecontrol, levels = design_2);
# Identify DEMiRNAs
fit_3 <- contrasts.fit(fit_2, contrast_matrix_2);
ebFit_2 <- eBayes(fit_3, trend = TRUE);
# Extract the top-ranked genes across conditions
dgeRes_miR <- topTable(ebFit_2, coef = 'HCMvsNF', number = nrow(miR_data));
#Annotate DEmiRNAs
miR_DEGM <- merge(dgeRes_miR, anno, by.x = "row.names", by.y = "PROBEID", all.x = TRUE);

#Extract list of up reg miRNAs and down reg miRNAs
miR_up <- miR_DEGM[miR_DEGM$P.Value < pval.cutoff & miR_DEGM$logFC > log2FC.cutoff,];
miR_down <- miR_DEGM[miR_DEGM$P.Value < pval.cutoff & miR_DEGM$logFC <  -log2FC.cutoff,];
DEmiR_HCM = rbind(miR_up,miR_down);
colnames(DEmiR_HCM)[1] <- "Probe_ID";
write.table(DEmiR_HCM, file = "DEmiR_HCM.txt", sep = "\t", row.names = F, quote = F)

#-----------------------------------------------------------------------------#
# Part 5: visualize the differentially expressed miRNAs in volcano plots

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!requireNamespace("ggrepel", quietly = TRUE)) BiocManager::install("ggrepel")
library(EnhancedVolcano)
library(ggrepel)


log2FC.cutoff <- 0.585;
pval.cutoff <- 0.05;
#Volcano plot of DEmiRNAs between HCM and NF
EnhancedVolcano(subset(miR_DEGM, select=c(2,3,5,6,8)), title = "miRNA HCM vs NF", 
                lab = NA, labSize = 3, x = 'logFC', 
                y = 'P.Value', pCutoff = pval.cutoff, FCcutoff = log2FC.cutoff,
                col = c("#999999", "#E69F00", "#56B4E9","#FF1493"))


#-----------------------------------------------------------------------------#
# Section 1.b: Differentially expressed miRNAs analysis and visualization 
# between DCM and NF
#-----------------------------------------------------------------------------#
# Part 1: Import all necessary files from GEO database

#Access and import GSE112556 dataset for HCM
gse_accession <- "GSE112556"
eset <- getGEO(gse_accession, GSEMatrix = TRUE, AnnotGPL = TRUE);
View(eset)

#extract miRNA expression dataset with available annotation
miR_data <- exprs(eset[[1]]);
View(miR_data)
miR_data <- as.data.frame(miR_data);
miR_data <- miR_data[-c(1:7),];
miR_DCM <- cbind(miR_data, rownames(miR_data));
write.table(miR_DCM, file = "DEmiR_DCM.txt", sep = "\t", row.names = F, quote = F)

#extract metadata file
sample_metadata <- pData(eset[[1]]);
sample_metadata$sample_type <- c("control", "control", "control","case", "case", "case");

#-----------------------------------------------------------------------------#
# Part 2: Quality control

# 2.1 Quality control - Box plots
#Box plots of disease group
all(rownames(sample_metadata) == colnames(miR_data));
plot_HCM <- gather(miR_data[, sample_metadata$sample_type == "case"], 
                   key = "SampleID", value = "log2CPM")
HCM_distr <- ggplot(plot_HCM, aes(x = SampleID, y = log2CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(miR_data[, sample_metadata$sample_type == "case"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution in the HCM group") # add title
HCM_distr

#Box plots of control group
plot_NF <- gather(miR_data[, sample_metadata$sample_type == "control"], 
                  key = "SampleID", value = "log2CPM")
NF_distr <- ggplot(plot_NF, aes(x = SampleID, y = log2CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(miR_data[, sample_metadata$sample_type == "case"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution of miRNA in the NF group") # add title
NF_distr

# 2.2 Quality control - density plot

#density plot
if (!requireNamespace("RColorBrewer", quietly = TRUE)) BiocManager::install("RColorBrewer")
library(RColorBrewer)

HCM_dens <- miR_data[, sample_metadata$sample_type == "case"] ;
nsamples <- ncol(HCM_dens) ;
# Create a color palette for all samples in this group
col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(nsamples);
# Density plots of HCM distribution
HCM_densityplot <- plot(density(HCM_dens[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.26), 
                        las = 2, main = "", xlab = "")
title(main="Density plots of log2CPMs distribution in the HCM group ", xlab="Log2cpm")
for (i in 2:nsamples){
  den <- density(HCM_dens[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

#-----------------------------------------------------------------------------#
# Part 4: Differentially expressed miRNAa analysis 

# Load required packages
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(edgeR)
library(limma)

# create a design matrix
design_2 <- model.matrix(~ 0 + sample_type , data = sample_metadata);
# Fit a linear model to the expression data for each gene 
fit_2 <- lmFit(miR_data, design_2);
# Create a contrast matrix for comparing expression levels between different conditions 
contrast_matrix_2<- makeContrasts(
  HCMvsNF = sample_typecase - sample_typecontrol, levels = design_2);
# Identify DEGs
fit_3 <- contrasts.fit(fit_2, contrast_matrix_2);
ebFit_2 <- eBayes(fit_3, trend = TRUE);
# Extract the top-ranked genes across conditions
dgeRes_miR <- topTable(ebFit_2, coef = 'HCMvsNF', number = nrow(miR_data));
miR_DEGM <- cbind(dgeRes_miR, rownames(dgeRes_miR));

#Extract list of up reg miRNAs and down reg miRNAs
miR_up <- miR_DEGM[miR_DEGM$P.Value < pval.cutoff & miR_DEGM$logFC > log2FC.cutoff,];
miR_down <- miR_DEGM[miR_DEGM$P.Value < pval.cutoff & miR_DEGM$logFC <  -log2FC.cutoff,];
miR_DCM =  rbind(miR_up,miR_down);
colnames(miR_DCM)[7] = "miR_ID";
write.table(miR_DCM, file = "DEmiR_DCM.txt", sep = "\t", row.names = F, quote = F)

#-----------------------------------------------------------------------------#
# Part 5: visualize the differentially expressed miRNAs in volcano plots

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!requireNamespace("ggrepel", quietly = TRUE)) BiocManager::install("ggrepel")
library(EnhancedVolcano)
library(ggrepel)


log2FC.cutoff <- 0.585;
pval.cutoff <- 0.05;
EnhancedVolcano(subset(miR_DEGM, select=c(7,1,4,5)), title = "miRNA DCM vs NF", 
                lab = NA, labSize = 3, x = 'logFC', 
                y = 'P.Value', pCutoff = pval.cutoff, FCcutoff = log2FC.cutoff,
                col = c("#999999", "#E69F00", "#56B4E9","#FF1493"))

#-----------------------------------End----------------------------------------#
