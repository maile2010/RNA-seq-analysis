#=============================================================================#
# Research Project 1
# Section 3:Pathway enrichment analysis
#																	                                       		   															                                #
# Date: January 24, 2025											                                
# Author: Mai Le, ID: i6375777
# Maastricht University                                                       #
#=============================================================================#

# Install required libraries (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rstudioapi", quietly = TRUE)) BiocManager::install("rstudioapi")
if (!requireNamespace("readxl", quietly = TRUE)) BiocManager::install("readxl")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!requireNamespace("VennDiagram", quietly = TRUE)) BiocManager::install("VennDiagram")
if (!requireNamespace("rWikiPathways", quietly = TRUE)) BiocManager::install("rWikiPathways")
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("RCy3", quietly = TRUE)) BiocManager::install("RCy3")

library(rstudioapi)
library(readxl)
library(EnhancedVolcano)
library(VennDiagram)
library(rWikiPathways)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(RCy3)

# Step 1 Retrieve the differentially expressed genes for DCM and HCM
setwd("D:/UM/Exp")
dataset_2 <- read.delim("DEGA_output_nw.txt", as.is = T);
dataset <- read_excel("magnet_de_dataset.xlsx");

setwd("D:/UM/Exp/Project Period 3/DEGs")
DEGs_DCM <- read.delim("DEGs_DCM.txt", as.is = T);
DEGs_HCM <- read.delim("DEGs_HCM.txt", as.is = T);

log2FC.cutoff <- 0.585;
pval.cutoff <- 0.05;

#Retrieve the up or down-regulated genes for both DCM and HCM
deg.DCM.up <- DEGs_DCM[DEGs_DCM$P.Value < pval.cutoff & DEGs_DCM$logFC > log2FC.cutoff,c(8,7)];
deg.DCM.down <- DEGs_DCM[DEGs_DCM$P.Value < pval.cutoff & DEGs_DCM$logFC < -log2FC.cutoff,c(8,7)];
deg.HCM.up <- DEGs_HCM[DEGs_HCM$P.Value_2 < pval.cutoff & DEGs_HCM$logFC_2 > log2FC.cutoff,c(8,7)];
deg.HCM.down <- DEGs_HCM[DEGs_HCM$P.Value_2 < pval.cutoff & DEGs_HCM$logFC_2 < -log2FC.cutoff,c(8,7)];

#-----------------------------------------------------------------------------#
# Step 2: Pathway collection and identifier mapping
# Option 1: Wikipathway
#We will perform pathway enrichment with the gene sets of all pathway models in WikiPathways (human only).

gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", date = "20241110", format = "gmt")
wp2gene <- readPathwayGMT(gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- unique(wp2gene %>% dplyr::select(wpid,name)) #TERM2NAME
bkgd.genes <- unique(dataset_2[,c(20,19)])

#The dataset has Ensembl identifiers but the pathway annotations are using Entrez Gene identifiers (NCBI Gene). 
#ClusterProfiler provides a function to map the identifiers.

bkgd.genes.entrez <- bitr(bkgd.genes$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
deg.HCM.entrez <- bitr(DEGs_HCM$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
deg.DCM.entrez <- bitr(DEGs_DCM$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- data.frame(ENTREZID = unique(wpid2gene$gene))

#-----------------------------------------------------------------------------#
# Option 2: KEGG Pathway enrichment analysis

library(clusterProfiler)

# For HCM 
ewp.HCM_2 <- enrichKEGG(
  gene = deg.HCM.entrez$ENTREZID,
  pAdjustMethod = "fdr",
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02
)
ewp.HCM.res_2 <- as.data.frame(ewp.HCM_2)

# For DCM
ewp.DCM_2 <- enrichKEGG(
  gene = deg.DCM.entrez$ENTREZID,
  pAdjustMethod = "fdr",
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02
)
ewp.DCM.res_2 <- as.data.frame(ewp.DCM_2)

# number of genes measured that are present in pathways
length(ewp.HCM_2@universe)
length(ewp.DCM_2@universe)

# number of DEG in pathways
length(deg.HCM.entrez$ENTREZID[deg.HCM.entrez$ENTREZID %in% unique(wp2gene$gene)])
length(deg.DCM.entrez$ENTREZID[deg.DCM.entrez$ENTREZID %in% unique(wp2gene$gene)])
num.pathways.HCM <- dim(ewp.HCM.res_2)[1]
num.pathways.DCM <- dim(ewp.DCM.res_2)[1]

# barchart for the top 10 overrepresented pathways
#ewp.HCM_2 <- ewp.HCM_2[order(ewp.HCM_2$pvalue), ]
#top_10_HCMpathways <- ewp.HCM_2[1:14, ]
#ggplot(top_10_HCMpathways, aes(x = reorder(Description, -pvalue), y = Count)) +
#  geom_bar(stat = "identity", fill = "#BA8CD7") +
#  coord_flip() +
#  labs(x = "", y = "HCM DEG gene count", fill = "") +
#  theme_minimal()

#For HCM 
ggplot(top_10_HCMpathways, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue)) +
  coord_flip() +
  labs(
    x = "",
    y = "-log10(p-value)",
    size = "Gene Count",
    color = "p-value"
  ) +
  theme_minimal() +
  scale_color_gradient(low = "#FF1493", high = "#0000FF")+ 
  theme(
    axis.text.x = element_text(size = 12),      # Increase x-axis text size
    axis.text.y = element_text(size = 12),      # Increase y-axis text size
    axis.title = element_text(size = 14),       # Increase axis title size
    legend.text = element_text(size = 12),      # Increase legend text size
    legend.title = element_text(size = 14),     # Increase legend title size
    plot.title = element_text(size = 16, face = "bold")  # Adjust title size if needed
  )

#For DCM
ewp.DCM_2 <- ewp.DCM_2[order(ewp.DCM_2$pvalue), ]
top_14_pathways <- ewp.DCM_2[1:14, ]
#ggplot(top_20_pathways, aes(x = reorder(Description, -pvalue), y = Count)) +
#  geom_bar(stat = "identity", fill = "#BA8CD7") +
#  coord_flip() +
#  labs(x = "", y = "DCM DEG gene count", fill = "") +
#  theme_minimal()

ggplot(top_14_pathways, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue)) +
  coord_flip() +
  labs(
    x = "",
    y = "-log10(p-value)",
    size = "Gene Count",
    color = "p-value"
  ) +
  theme_minimal() +
  scale_color_gradient(low = "#FF1493", high = "#0000FF") +
  theme(
    axis.text.x = element_text(size = 12),      # Increase x-axis text size
    axis.text.y = element_text(size = 12),      # Increase y-axis text size
    axis.title = element_text(size = 14),       # Increase axis title size
    legend.text = element_text(size = 12),      # Increase legend text size
    legend.title = element_text(size = 14),     # Increase legend title size
    plot.title = element_text(size = 16, face = "bold")  # Adjust title size if needed
  )

# another visualization is a treeplot that also shows similarity between pathways (gene overlap)
library(enrichplot)
ewp.DCM.sim <- enrichplot::pairwise_termsim(ewp.DCM)
treeplot(ewp.DCM.sim, label_format = 0.4)

ewp.HCM.sim <- enrichplot::pairwise_termsim(ewp.HCM_2)
treeplot(ewp.HCM.sim, label_format = 0.4)

#-----------------------------------------------------------------------------#
# Part 3: pathway enrichment analysis for DEGs in hub module

#Extract DEGs in huhb module
deg.hub.module.HCM.entrez <- bitr(node_HCM_FC$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
deg.hub.module.DCM.entrez <- bitr(node_DCM_FC$ensembl_gene_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Pathway enrichment analysis for hub module of HCM
ewp.DEGs_HCM <- clusterProfiler::enricher(
  deg.hub.module.HCM.entrez$ENTREZID,
  universe = bkgd.genes.entrez$ENTREZID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.DEGs_HCM.res <- as.data.frame(ewp.DEGs_HCM)

#Pathway enrichment analysis for hub module in DCM
ewp.DEGs_DCM <- clusterProfiler::enricher(
  deg.hub.module.DCM.entrez$ENTREZID,
  universe = bkgd.genes.entrez$ENTREZID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.DEGs_DCM.res <- as.data.frame(ewp.DEGs_DCM)

# number of genes measured that are present in pathways
length(ewp.DEGs_HCM@universe)
length(ewp.DCM@universe)

# number of DEG in pathways
length(deg.hub.module.HCM.entrez$ENTREZID[deg.hub.module.HCM.entrez$ENTREZID %in% unique(wp2gene$gene)])
length(deg.DCM.entrez$ENTREZID[deg.DCM.entrez$ENTREZID %in% unique(wp2gene$gene)])
num.pathways.HCM_hubmodule <- dim(ewp.DEGs_HCM.res)[1]
num.pathways.DCM <- dim(ewp.DCM.res)[1]

# barchart for the overrepresented pathways
ggplot(ewp.DEGs_HCM[1:num.pathways.HCM_hubmodule,], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="HCM DEG gene count in hub module", fill="")
write.table(ewp.male.res, file="DCM_Male_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)

#-----------------------------------End----------------------------------------#

