#=============================================================================#
# Research Project 1
# Section 2: Differentially expressed gene analysis and visualization
#																	                                       		   															                                #
# Date: January 24, 2025											                                
# Author: Mai Le, ID: i6375777
# Maastricht University                                                       #
#=============================================================================#

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
#-----------------------------------------------------------------------------#
#Part 1:Retrieve gene symbols and gene names based on Ensembl gene ID

# Load required packages
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
library(org.Hs.eg.db)
library(biomaRt)

gene_IDs <- rownames(expression);
# List available database in biomaRt
avail_dataset <- listEnsembl();
# Select the available database to use in biomaRt
ensembl <- useEnsembl(biomart = "ensembl", mirror = "www");
# List the available data and using it from biomaRt
dataset <- listDatasets(ensembl);
# Connect the dataset by using useMart function
ensembl_conn <- useEnsembl(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl",mirror = "www");
# To see attributes and filters in biomaRt database
attr <- listAttributes(ensembl_conn);
filters <- listFilters(ensembl_conn);
# Build the query
annotated_ids2 <- getBM(attributes = c("uniprot_gn_symbol","ensembl_gene_id","description"), 
                        filters = "ensembl_gene_id", 
                        values = gene_IDs, 
                        mart = ensembl_conn)
# Ensure unique rows in annotation data - avoid duplicate entries in the Annotation files
annotated_ids2 <- annotated_ids2[!duplicated(annotated_ids2$ensembl_gene_id), ];

# Merge this annotation with expression data based on the same gene_id
expression_2 <- merge(expression,annotated_ids2, by.x = "row.names", 
                      by.y = "ensembl_gene_id", all.x = TRUE);
# Remove the first column - gene IDs 
annotated_expression <- expression_2
write.table(annotated_expression, file = "annotated_expression.txt", sep = "\t", row.names = F, quote = F)

#-----------------------------------------------------------------------------#
# Part 2: Differentially expressed gene analysis

# Load required packages
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(edgeR)
library(limma)

# Create a design matrix that is corrected for covariates: etiology, gender and
# age (reason is explained above)
design <- model.matrix(~ 0  + etiology+ gender + age, data = sampleInfo);
# Fit a linear model to the expression data for each gene 
fit <- lmFit(expression, design);
# Create a contrast matrix for comparing expression levels between different conditions 
contrast_matrix <- makeContrasts(
  DCMvsNF = etiologyDCM - etiologyNF, #compare DCM and NF
  HCMvsNF = etiologyHCM - etiologyNF, # compare HCM and NF
  PPCMvsNF = etiologyPPCM - etiologyNF, # compare PPCM and NF
  levels = design);
# Identify DEGs
fit2 <- contrasts.fit(fit, contrast_matrix);
ebFit <- eBayes(fit2, trend = TRUE);
# Extract the top-ranked genes across conditions
dgeRes_1 <- topTable(ebFit, coef = 'DCMvsNF', number = nrow(expression));
dgeRes_2 <- topTable(ebFit, coef = 'HCMvsNF', number = nrow(expression));
dgeRes_3 <- topTable(ebFit, coef = 'PPCMvsNF', number = nrow(expression));
dgeRes <- cbind(dgeRes_1, dgeRes_2, dgeRes_3)

rownames(annotated_ids2) <- annotated_ids2$ensembl_gene_id;
# Rearrange the annotation data with the same order of gene_id as dgeRes data
annotated_ids2 <- annotated_ids2[rownames(dgeRes), ];
network_dataset <- cbind(dgeRes, annotated_ids2)
#export data file for network analysis 
write.table(network_dataset, file = "network_dataset.txt", sep = "\t", row.names = F, quote = F)
network_dataset <- read.delim("network_dataset.txt", as.is = T)
                                 
#-----------------------------------------------------------------------------#
# Part 3: visualize the differentially expressed genes in volcano plots

log2FC.cutoff <- 0.585;
pval.cutoff <- 0.05;

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!requireNamespace("VennDiagram", quietly = TRUE)) BiocManager::install("VennDiagram")
if (!requireNamespace("grid", quietly = TRUE)) BiocManager::install("grid")
if (!requireNamespace("futile.logger", quietly = TRUE)) BiocManager::install("futile.logger")

library(EnhancedVolcano)
library(VennDiagram)
library(grid)
library(futile.logger)

colnames(network_dataset) <- c("logFC", "AveExpr" , "t" , "P.Value", "adj.P.Val",        
                       "B","logFC_2","AveExpr_2","t_2", "P.Value_2",          
                       "adj.P.Val_2","B_2","logFC_3","AveExpr_3","t_3",                
                       "P.Value_3", "adj.P.Val_3","B_3","uniprot_gn_symbol","ensembl_gene_id",  
                       "description");
#Volcano plot of DEGs between DCM and NF
EnhancedVolcano(subset(network_dataset, select=c(19,20,1,4,5)), title = "DCM vs NF", 
                lab = NA, labSize = 3, x = 'logFC', 
                y = 'P.Value', pCutoff = pval.cutoff, FCcutoff = log2FC.cutoff,
                col = c("#999999", "#E69F00", "#56B4E9","#FF1493")) # Color-blind-friendly scheme
#Volcano plot of DEGs between HCM and NF
EnhancedVolcano(subset(network_dataset, select=c(19,20,7,10,11)), title = "HCM vs NF", 
                lab = NA, labSize = 3, x = 'logFC_2', 
                y = 'P.Value_2', pCutoff = pval.cutoff, FCcutoff = log2FC.cutoff,
                col = c("#999999", "#E69F00", "#56B4E9","#FF1493"))

#-----------------------------------------------------------------------------#
# Part 4: visualize the overlaps between differentially expressed genes  and 
# targets of differentially expressed miRNAs in venn diagrams

# retrieve the differentially expressed genes for HF and NF
deg.DCM <- network_dataset[network_dataset$P.Value < 
                     pval.cutoff & abs(network_dataset$logFC) 
                   > log2FC.cutoff,c(1,2,3,4,5,6,19,20,21)];
deg.HCM <- network_dataset[network_dataset$P.Value_2 < 
                     pval.cutoff & abs(network_dataset$logFC_2) 
                   > log2FC.cutoff,c(7,8,9,10,11,12,19,20,21)];

#remove NA and blank values
deg.HCM_1 <-deg.HCM[!is.na(deg.HCM$uniprot_gn_symbol), ];
deg.HCM_clean <- deg.HCM_1[deg.HCM_1$uniprot_gn_symbol  != "", ];
write.table(deg.HCM_clean, file = "DEGs_HCM.txt", sep = "\t", row.names = F, quote = F)

deg.DCM_1 <-deg.DCM[!is.na(deg.DCM$hgnc_symbol), ];
deg.DCM_clean <- deg.DCM_1[deg.DCM_1$hgnc_symbol  != "", ];
write.table(deg.DCM_clean, file = "DEGs_DCM.txt", sep = "\t", row.names = F, quote = F)

#extract significant DEGs after removing NA and blank values
deg.DCM.down <- deg.DCM_clean[deg.DCM_clean$P.Value < pval.cutoff & 
                                deg.DCM_clean$logFC <  -log2FC.cutoff,];
deg.DCM.up <- deg.DCM_clean[deg.DCM_clean$P.Value < pval.cutoff & 
                              deg.DCM_clean$logFC > log2FC.cutoff,];
deg.HCM.up <- deg.HCM_clean[deg.HCM_clean$P.Value_2 < pval.cutoff & 
                              deg.HCM_clean$logFC_2 > log2FC.cutoff,];
deg.HCM.down <- deg.HCM_clean[deg.HCM_clean$P.Value_2 < pval.cutoff 
                                & deg.HCM_clean$logFC_2 < -log2FC.cutoff,];

# Part 4.1: HCM 
#load DEM targets files
setwd("D:/UM/Exp/Project Period 3/Target");
up_DEM_HCM <-  read.csv("up_DEMs_HCM.csv", as.is = T);
down_DEM_HCM <-  read.csv("down_DEMs_HCM.csv", as.is = T);
target_HCM = rbind(up_DEM_HCM,down_DEM_HCM);

#find dual-role genes which are targets of both up and down-regulated miRNAs
dual.role_HCM = Reduce(intersect, list(down_DEM_HCM$DEMs, up_DEM_HCM$DEMs));

#filter DEGs which overlap with target_HCM
overlap.DEGs <-  deg.HCM_clean %>%
filter(deg.HCM_clean$uniprot_gn_symbol %in% target_HCM$DEMs);
write.table(overlap.DEGs, file = "overlap.DEGs_HCM.txt", sep = "\t", row.names = F, quote = F)

#filter overlapping upregulated DEGs with target_HCM
overlap.up.DEGs <-  deg.HCM.up %>%
  filter(deg.HCM.up$uniprot_gn_symbol %in% up_DEM_HCM$DEMs);

#filter overlapping downregulated DEGs with target_HCM
overlap.down.DEGs <-  deg.HCM.down %>%
  filter(deg.HCM.down$uniprot_gn_symbol %in% down_DEM_HCM$DEMs);

#merge overlapping up and downregulated DEGs into overlapping DEGs
overlap_DEGs = rbind(overlap.up.DEGs,overlap.down.DEGs);
write.table(overlap_DEGs, file = "overlap_DEGs_HCM.txt", sep = "\t", row.names = F, quote = F);

# check output file - venn_genes.png
library(VennDiagram)
library(scales)

#Venn diagram to find overlap between up reg DEGs and DEM targets
venn.diagram(
  x = list(deg.HCM.up$uniprot_gn_symbol, up_DEM_HCM$DEMs),
  category.names = c(" DEGs", " miR targets"),
  filename = 'up_overlap_venn.png',
  output=FALSE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 1.5)

#Venn diagram to find overlap between down reg DEGs and DEM targets
venn.diagram(
  x = list(deg.HCM.down$uniprot_gn_symbol, down_DEM_HCM$DEMs),
  category.names = c(" DEGs", " miR targets"),
  filename = 'down_overlap_vene.png',
  output=FALSE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 1.5)

#Overlap of both up and down regulated genes with DEM targets
venn.diagram(
  x = list(deg.HCM.up$uniprot_gn_symbol, deg.HCM.down$uniprot_gn_symbol, 
           up_DEM_HCM$DEMs, down_DEM_HCM$DEMs),
  category.names = c("Up DEGs"," Down DEGs","Up miR Targets","Down miR Targets"),
  filename = 'venn_overlapping_genes.png',
  output=FALSE,
  col=c("#440154ff","#440154ff", '#21908dff','#21908dff' ),
  fill = c(alpha("#440154ff",0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3),alpha('#21908dff',0.3)),
  cex = 1.5)

#Find overlap genes
#overlap_up_HCM <- Reduce(intersect, list(deg.HCM_clean$uniprot_gn_symbol, up_DEM_HCM$DEMs))
#write.table(overlap_up_HCM, file = "overlap_upregulation_HCM.txt", sep = "\t", row.names = F, quote = F)
#overlap_down_HCM <- Reduce(intersect, list(deg.HCM_clean$uniprot_gn_symbol, down_DEM_HCM$DEMs))
#write.table(overlap_down_HCM, file = "overlap_downregulation_HCM.txt", sep = "\t", row.names = F, quote = F)

# Find overlap miRs
#overlap_miR_HCM <- Reduce(intersect, list(up_DEM_HCM$DEMs, down_DEM_HCM$DEMs))
#write.table(overlap_miR_HCM, file = "overlap_miRNA_HCM.txt", sep = "\t", row.names = F, quote = F)


#Part 4.2 DCM

#load DEM targets files of DCM
setwd("D:/UM/Exp/Project Period 3/Target");
up_DEM_DCM <-  read.csv("up_DEM_DCM.csv", as.is = T);
down_DEM_DCM <-  read.csv("down_DEM_DCM.csv", as.is = T);
target_DCM = rbind(up_DEM_DCM,down_DEM_DCM);

#find dual-role genes which are targets of both up and down-regulated miRNAs
dual.role_DCM = Reduce(intersect, list(up_DEM_DCM$DEMs, down_DEM_DCM$DEMs))

#filter overlapping upregulated DEGs with target_HCM
overlap.up.DEGs_DCM <-  deg.DCM.up %>%
  filter(deg.DCM.up$uniprot_gn_symbol %in% up_DEM_DCM$DEMs);

#filter overlapping downregulated DEGs with target_HCM
overlap.down.DEGs_DCM <-  deg.DCM.down %>%
  filter(deg.DCM.down$uniprot_gn_symbol %in% down_DEM_DCM$DEMs);

#merge overlapping up and downregulated DEGs into overlapping DEGs
overlap_DEGs_DCM = rbind(overlap.up.DEGs_DCM,overlap.down.DEGs_DCM);
write.table(overlap_DEGs_DCM, file = "overlap_DEGs_DCM.txt", sep = "\t", row.names = F, quote = F)

# check output file - venn_genes.png
library(VennDiagram)
library(scales)

#Venn diagram to find overlap between up reg DEGs and DEM targets
venn.diagram(
  x = list(deg.DCM.up$uniprot_gn_symbol, up_DEM_DCM$DEMs),
  category.names = c(" DEGs", " miR targets"),
  filename = 'up_overlap_venn_DCM.png',
  output=FALSE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 1.5)

#Venn diagram to find overlap between down reg DEGs and DEM targets
venn.diagram(
  x = list(deg.DCM.down$uniprot_gn_symbol, down_DEM_DCM$DEMs),
  category.names = c(" DEGs", " miR targets"),
  filename = 'down_overlap_vene_DCM.png',
  output=FALSE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 1.5)

#Overlap of both up and down regulated genes
venn.diagram(
  x = list(deg.DCM.up$uniprot_gn_symbol, deg.DCM.down$uniprot_gn_symbol, 
           up_DEM_DCM$DEMs, down_DEM_DCM$DEMs),
  category.names = c("Up DEGs"," Down DEGs","Up miR Targets","Down miR Targets"),
  filename = 'venn_overlapping_genes_DCM.png',
  output=FALSE,
  col=c("#440154ff","#440154ff", '#21908dff','#21908dff' ),
  fill = c(alpha("#440154ff",0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3),alpha('#21908dff',0.3)),
  cex = 1.5)

#-----------------------------------End----------------------------------------#

