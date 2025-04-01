#=============================================================================#
# MSB1005_Final_Assignment.R                                                  #
#																	                                       		   															                                #
# Date: December 12, 2024											                                
# Author: Mai Le, ID: i6375777
# Maastricht University                                                       #
#=============================================================================#

#-----------------------------------------------------------------------------#
# Question 1: Import the data and export a publication-ready table
#-----------------------------------------------------------------------------#

# Question 1.a: Import all the data files
#First, set the active working directory to the folder containing the files
setwd("D:/UM/Exp")
#Import all necessary files
expression <-  read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", as.is = T, 
                          row.names = 1);
geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, 
                                 row.names = 1);
sampleInfo <-  read.csv("MAGNET_SampleData_18112022.csv", as.is = T, 
                        row.names = 1);
InfoDescription <- read.csv("MAGNET_SampleData_18112022_WithDescriptions.csv", as.is = T, 
                            row.names = 1);

# Question 1.b:Export a publication-ready table of participant characteristics, 
# including statistics comparing the four etiologies

# To create a publication-ready table, I use the 'gtsummary' package because
# it provides an elegant and flexible way to create publication-ready analytical 
# and summary tables using the R programming language. 

if (!requireNamespace("gtsummary", quietly = TRUE)) BiocManager::install("gtsummary")
if (!requireNamespace("flextable", quietly = TRUE)) BiocManager::install("flextable")
if (!requireNamespace("cardx", quietly = TRUE)) BiocManager::install("cardx")
if (!requireNamespace("tidyverse", quietly = TRUE)) BiocManager::install("tidyverse")
if (!requireNamespace("officer", quietly = TRUE)) BiocManager::install("officer")

# load required packages
library(gtsummary) 
library(flextable)
library(cardx)
library(tidyverse) 
library(officer)

# First, I set gtsummary theme to control many aspects of how a table is printed
# to meet specific publication standards and user preferences.
gtsummary::theme_gtsummary_compact(set_theme = TRUE, font_size = 9)

# Create participant characteristics table 
adv <- sampleInfo %>%
  tbl_summary(
    by = etiology,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous() ~ c(
        "{mean} ({sd})",
        "{median}",
        "{min}, {max}"
      ),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    missing = "no",
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels() %>%
  modify_header(
    label = "",
    all_stat_cols() ~ "**{level}**, N = {n} ({style_percent(p)}%)"
  ) %>%
  gtsummary::as_flex_table()

adv %>% 
  flextable::line_spacing(space = 1.1, part = "body") %>% 
  flextable::autofit()

# Typically, a participant characteristics table will need to be inserted into a
# Microsoft Word document where it can join the rest of a manuscript or report. 
# Thus, I use the 'gtsummary' package, with the 'flextable' and 'officer' packages,
# to export this tables to Microsoft Word for further publication purpose. 

# To export publication-ready table to Microsoft Word
flextable::save_as_docx(adv,
                        path = "Characteristics_Table.docx",
                        pr_section =
                          officer::prop_section(
                            page_size = officer::page_size(
                              orient = "landscape",
                              width = 8.3, height = 11.7
                            ),
                            type = "continuous",
                            page_margins = officer::page_mar()
                          )
)

#-----------------------------------------------------------------------------#
# Question 2: Diagnostic plots
#-----------------------------------------------------------------------------#

# Question 2.a: Plot data distribution figures (at least one figure) that enables
# comparing samples

# Check whether the order of the samples is the same in two objects:
all(rownames(sampleInfo) == colnames(expression)) # TRUE, so we can easily combine

# First, I use box plots to check the distribution of the CPMs on the log2 scale
# This plot is also used to check the quality control in each group by comparing
# the median of log2-transformed CPMs of each sample with the median of 
# that whole group.

# load required package
if (!requireNamespace("tidyr", quietly = TRUE)) BiocManager::install("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library(tidyr)
library(ggplot2)

# Boxplot_1: Check distribution of all samples in NF group using boxplots and 
# compare with median value of this group by adding a blue horizontal line that
# corresponds to the median log2CPM 
plotData_NF <- gather(expression[, sampleInfo$etiology == "NF"], 
                   key = "SampleID", value = "CPM")
NF_dist <- ggplot(plotData_NF, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(expression[, sampleInfo$etiology == "NF"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution in the NF group") # add title

# Boxplot_2: Check distribution of all samples for quality control in the DCM group 
plotData_DCM <- gather(expression[, sampleInfo$etiology == "DCM"], 
                      key = "SampleID", value = "CPM")
DCM_dist <- ggplot(plotData_DCM, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(expression[, sampleInfo$etiology == "DCM"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution in the DCM group") # add title

# Boxplot_3: Check distribution of all samples and quality control in the HCM group 
plotData_HCM <- gather(expression[, sampleInfo$etiology == "HCM"], 
                       key = "SampleID", value = "CPM")
HCM_dist <- ggplot(plotData_HCM, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(expression[, sampleInfo$etiology == "HCM"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution in the HCM group") # add title

# Boxplot_4: Check distribution of all samples and quality control in the PPCM group 
plotData_PPCM <- gather(expression[, sampleInfo$etiology == "PPCM"], 
                       key = "SampleID", value = "CPM")
PPCM_dist <- ggplot(plotData_PPCM, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic() + 
  geom_hline(
    yintercept = median(as.numeric(as.matrix(expression[, sampleInfo$etiology == "PPCM"]), 
                                   na.rm = TRUE)), color = "blue"
  ) + # add the horizontal line that corresponds to the median log2CPM
  ggtitle("Boxplots of log2CPMs distribution in the PPCM group") # add title

# Export figures as PDF files
ggsave("Boxplot_logCPMs_Distribution_of_NF.pdf", plot = NF_dist, 
       width = 10, height = 8, units = "in") #for NF distribution
ggsave("Boxplot_logCPMs_Distribution_of_DCM.pdf", plot = DCM_dist, 
       width = 10, height = 8, units = "in") #for DCM distribution
ggsave("Boxplot_logCPMs_Distribution_of_HCM.pdf", plot = HCM_dist, 
       width = 10, height = 8, units = "in") #for HCM distribution
ggsave("Box_plotlogCPMs_Distribution_of_PPCM.pdf", plot = PPCM_dist, 
       width = 10, height = 8, units = "in") #for HCM distribution

# From the boxplots, we see that the overall of the density distribution of 
# log2-transformed CPMs are not identical across samples but still not very 
# different, thus these group have good qualities 

# Next, I use the density plots to see the distribution of log2CPM values for all
# samples in the NF group
if (!requireNamespace("RColorBrewer", quietly = TRUE)) BiocManager::install("RColorBrewer")
library(RColorBrewer)

NF <- expression[, sampleInfo$etiology == "NF"] 
nsamples <- ncol(NF) 
# Create a color palette for all samples in this group
col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(nsamples)
# Density plots of NF distribution
NF_dens <- plot(density(NF[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.26), 
                las = 2, main = "", xlab = "")
title(main="Density plots of log2CPMs distribution in the NF group ", xlab="Log2cpm")
for (i in 2:nsamples){
  den <- density(NF[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

# Export figures as PDF files
ggsave("DensityPlot_logCPMs_Distribution_of_NF.pdf", plot = NF_dens, 
       device = "pdf", width = 11.69, height = 8.27, units = "in") 

# Do the same code for other groups to see their density plots

# Question 2.b: Plot at least one PCA figure showing the sample clustering 
# colored by relevant covariates.

# Install and load required package
if (!requireNamespace("pcaMethods", quietly = TRUE)) BiocManager::install("pcaMethods")
library(pcaMethods)
# Perform a Principal Component Analysis
pcaRes <- pca(t(expression), nPcs = 10)
plot(pcaRes)
# Visualizes the PCA scores based on the first two principal components
plotPcs(pcaRes, c(1,2))
# Again check that the order of samples is the same
all(rownames(pcaRes@scores) == rownames(sampleInfo)) # TRUE
# Create ggplot plotting data
PCA_Data <- cbind(data.frame(pcaRes@scores), sampleInfo)

# Why do I choose etiology, age, and gender, and race as covariates for CVD?
# Numerous studies about the epidemiology of CVD suggested that the etiology 
# of heart failure varies significantly among different populations and genders, 
# influencing both treatment approaches and prognoses.  
# Moreover, age is another critical covariate since the risk of developing heart 
# failure increases significantly with advancing age. Race has been shown to 
# affect both the prevalence of heart failure and its outcomes. 
# For gender, women typically present with CVD at an older age than men and often
# have different underlying causes and symptom profiles. Therefore, incorporating
# etiology, age, gender, and race as covariates in CVD studies is vital for
# understanding the multifaceted nature of this condition and reducing the 
# confounding risk in outcome representation. 

# Create PCA plots colored by covariates
etio <- ggplot(PCA_Data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(size = age, col = etiology)) + # colored by covariate etiology
  labs(title = "PCA of Gene Expression Data")
gend <- ggplot(PCA_Data, aes(x = PC1, y = PC3)) + 
  geom_point(aes(size = age, col = gender)) + # colored by covariate gender
  labs(title = "PCA of Gene Expression Data")
race <- ggplot(PCA_Data, aes(x = PC2, y = PC3)) + 
  geom_point(aes(size = age, col = race)) + # colored by covariate race
  labs(title = "PCA of Gene Expression Data")

# Export the figures to PDF files
ggsave("PCA_etiology.pdf", plot = etio, width = 10, 
       height = 8, units = "in") 
ggsave("PCA_gender.pdf", plot = gend, width = 10, 
       height = 8, units = "in") 
ggsave("PCA_race.pdf", plot = race, width = 10, 
       height = 8, units = "in") 

#-----------------------------------------------------------------------------#
# Question 3: Statistical analysis: Perform a differential gene expression 
# analysis and correct for relevant co-variates. 
#-----------------------------------------------------------------------------#

# Load required packages
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(edgeR)
library(limma)

# Create a design matrix that is corrected for covariates: etiology, gender and
# age (reason is explained above)
design <- model.matrix(~ 0 + etiology + gender + age, data = sampleInfo)
# Fit a linear model to the expression data for each gene 
fit <- lmFit(expression, design)
# Create a contrast matrix for comparing expression levels between different conditions 
contrast_matrix <- makeContrasts(
  DCMvsNF = etiologyDCM - etiologyNF, #compare DCM and NF
  HCMvsNF = etiologyHCM - etiologyNF, # compare HCM and NF
  PPCMvsNF = etiologyPPCM - etiologyNF, # compare PPCM and NF
  levels = design)
# Identify DEGs
fit2 <- contrasts.fit(fit, contrast_matrix)
ebFit <- eBayes(fit2, trend = TRUE)
# Extract the top-ranked genes across conditions
dgeRes <- topTable(ebFit, coef = c('DCMvsNF','HCMvsNF','PPCMvsNF'), number = nrow(expression))

# Retrieve the top 200 DEGs and all DEGs for further analysis
write.table(cbind(rownames(dgeRes)[1:200]), "top200.txt", sep = "\t", 
            row.names = F, quote = F)
write.table(cbind(rownames(dgeRes)), "background.txt", sep = "\t", 
            row.names = F, quote = F)

#-----------------------------------------------------------------------------#
# Question 4: Gene annotation: Retrieve gene symbols and gene names based on 
# Ensembl gene IDs
#-----------------------------------------------------------------------------#

# Question 4a: Retrieve gene symbols and gene names based on Ensembl gene ID

# Load required packages
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
library(org.Hs.eg.db)
library(biomaRt)

gene_IDs <- rownames(expression)
# List available database in biomaRt
avail_dataset <- listEnsembl()
# Select the available database to use in biomaRt
ensembl <- useEnsembl(biomart = "ensembl", mirror = "www")
# List the available data and using it from biomaRt
datasets <- listDatasets(ensembl)
# Connect the dataset by using useMart function
ensembl_conn <- useEnsembl(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl",mirror = "www")
# To see attributes and filters in biomaRt database
attr <- listAttributes(ensembl_conn)
filters <- listFilters(ensembl_conn)
# Build the query
annotated_ids2 <- getBM(attributes = c("uniprot_gn_symbol","ensembl_gene_id","description"), 
                        filters = "ensembl_gene_id", 
                        values = gene_IDs, 
                        mart = ensembl_conn)

# Ensure unique rows in annotation data - avoid duplicate entries in the Annotation files
annotated_ids2 <- annotated_ids2[!duplicated(annotated_ids2$ensembl_gene_id), ]

# Question 4b: Merge this annotation with the gene expression data object

# Merge this annotation with expression data based on the same gene_id
expression_2 <- merge(expression,annotated_ids2, by.x = "row.names", 
                      by.y = "ensembl_gene_id", all.x = TRUE)
# Remove the first column - gene IDs 
annotated_expression <- expression_2[,-1]

#-----------------------------------------------------------------------------#
# Question 5: Relative expression levels
#-----------------------------------------------------------------------------#

# Question 5a: Transform the data to FPKM values

all(rownames(geneTotExonLengths) == rownames(expression)) # TRUE (just a check)
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
  #	.t <- 2^(x) * 1E3 / geneTotExonLengths[, 1] # this does the same, but shorter
}
expression_fpkm <- cpm2fpkm(expression)

# Question 5b+c: Assess for each gene in the dataset whether it is expressed above
# background (noise) level - which is the average expression of Y chromosome genes 
# in female subjects.

# First, I identify the noise level, which is the average expression of all genes 
# located in Y chromosome and expressed in females. 

# Retrieve list of genes on chromosome Y from Ensembl database
genes_Y <- getBM(attributes = c("chromosome_name", "ensembl_gene_id", "hgnc_symbol"),
                 filters = "chromosome_name", values = "Y", mart = ensembl_conn)
rownames(genes_Y) <- genes_Y$ensembl_gene_id

# Filter these Y chromosome genes among all sequenced genes 
exp_genes_Y <- expression_fpkm[rownames(expression_fpkm)  %in% rownames(genes_Y),]

# Because exp_genes_Y data includes the fpkm values of all participants - males and
# females, I will filter females participants from this data

# Filter by gender to find the Y chromosome genes expressed in females 
exp_genes_Y_females <- exp_genes_Y[, colnames(exp_genes_Y) %in% 
                      rownames(sampleInfo[sampleInfo$gender == "Female",])]
# Calculate average expression of Y chromosome genes in females as noise 
noise <- mean(rowMeans(exp_genes_Y_females))
# Compare the average expression of each gene to the noise value
above_noise <- rowMeans(expression_fpkm) > noise
# Merge comparison result with the expression_fpkm data
expression_above_noise <- data.frame(expression_fpkm,above_noise)

#-----------------------------------------------------------------------------#
# Question 6: Export the results
#-----------------------------------------------------------------------------#

# Question 6b: The file should contain all the additionally generated data 
# (e.g. statistics,annotation, expressed above background) with clear column names.

# Statistics data after DEGA is dgeRes data
# Annotation data is annotated_ids2 
# Expression above background data is in the expression_above_noise data (both
# fpkm values and above_noise status)

# First, I will create a dataframe which include generated statistics data, 
# annotation data, and expressed above background data,
# Because the orders of gene_id in each data are different, I also need to 
# rearrange them in the same order as the dgeRes data (as this data is in order
# of the top-ranked differentially expressed genes)

rownames(annotated_ids2) <- annotated_ids2$ensembl_gene_id
# Rearrange the annotation data with the same order of gene_id as dgeRes data
annotated_ids2 <- annotated_ids2[rownames(dgeRes), ]
# Rearrange the expression_above_noise data with the same order of gene_id as dgeRes data
expression_above_noise <- expression_above_noise[rownames(dgeRes), ]
# Extract above_noise results
Expression_above_Background <- expression_above_noise$above_noise

# Question 6c: Include the average expression value (fpkm) for DCM, HCM, PPCM, 
# and controls of each gene into export file for further visualization purposes

# Average expression values of each gene for each group
average_DCM <- rowMeans(expression_above_noise[, colnames(expression_above_noise) %in% 
                                     rownames(sampleInfo[sampleInfo$etiology == "DCM",])])
average_HCM <- rowMeans(expression_above_noise[, colnames(expression_above_noise) %in% 
                                                 rownames(sampleInfo[sampleInfo$etiology == "HCM",])])
average_PPCM <- rowMeans(expression_above_noise[, colnames(expression_above_noise) %in% 
                                                 rownames(sampleInfo[sampleInfo$etiology == "PPCM",])])
average_NF <- rowMeans(expression_above_noise[, colnames(expression_above_noise) %in% 
                                                  rownames(sampleInfo[sampleInfo$etiology == "NF",])])

# Question 6c: Export file with results in a tab-delimited text file.

# Combine all above data to the export file of results
export_file <- cbind(dgeRes,average_DCM,average_HCM, average_PPCM,
                          average_NF,Expression_above_Background, annotated_ids2)

# Specify the file path where you want to export the file
output_file_path <- "export_file.txt"

# Write the data frame to a tab-delimited text file
write.table(export_file, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

#-----------------------------------End----------------------------------------#

