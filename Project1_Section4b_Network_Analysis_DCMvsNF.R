#=============================================================================#
# Research Project 1
# Section 4b: Network analysis: WGCNA, PPI network, mRNA-miRNA interaction network - DCM vs NF
#																	                                       		   															                                #
# Date: January 24, 2025											                                
# Author: Mai Le, ID: i6375777
# Maastricht University                                                       #
#=============================================================================#

#-----------------------------------------------------------------------------#
# Part 1: WGCNA and PPI network for DEGs between DCM and NF
#-----------------------------------------------------------------------------#
setwd("D:/UM/Exp")
expression <-  read.delim("annotated_expression.txt", as.is = T, 
                          row.names = 1);
sampleInfo <-  read.csv("MAGNET_SampleData_18112022.csv", as.is = T, 
                        row.names = 1);

#load list of DEGs
setwd("D:/UM/Exp/Project Period 3/DEGs")
DEGs_DCM <- read.delim("DEGs_DCM.txt", as.is = T);

expression_DCM <- expression[DEGs_DCM$ensembl_gene_id,]
expression_DCM$uniprot_gn_symbol = DEGs_DCM$uniprot_gn_symbol

# Load required package
options(stringsAsFactors = F)
if (!requireNamespace("writexl", quietly = TRUE)) BiocManager::install("writexl")
if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")
if (!requireNamespace("rstudioapi", quietly = TRUE)) BiocManager::install("rstudioapi")
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr")
if (!requireNamespace("RCy3", quietly = TRUE)) BiocManager::install("RCy3")
if (!requireNamespace("readr", quietly = TRUE)) BiocManager::install("readr")
if (!requireNamespace("readxl", quietly = TRUE)) BiocManager::install("readxl")

library(writexl)
library(WGCNA)
library(rstudioapi)
library(dplyr)
library(RCy3)
library(readr)
library(readxl)

if (!requireNamespace("tidyr", quietly = TRUE)) BiocManager::install("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library(tidyr)
library(ggplot2)

#-----------------------------------------------------------------------------#
# Step 1: extract expression data and sample information of HCM and NF groups

#extract expression data
data_DCM <- expression_DCM[, sampleInfo$etiology == "DCM"]
data_NF <- expression_DCM[, sampleInfo$etiology == "NF"]
data_NF = data_NF[,-c(167,168)]
data.log_DCM <- rbind(t(data_DCM), t(data_NF))

#extract metadata 
traitdata_DCM <- sampleInfo[rownames(data.log_DCM),c(2,3,5)]
traitdata_DCM$sampleID <- rownames(traitdata_DCM)

#-----------------------------------------------------------------------------#
# Step 2: Form a data frame analogous to expression data that will hold the clinical traits.
sample_DCM <- rownames(data.log_DCM)
traitrows_DCM = match(sample_DCM, traitdata_DCM$sampleID);
dat_traits_DCM = traitdata_DCM[traitrows_DCM, -4];
rownames(dat_traits_DCM) = traitrows_DCM

dat_traits_DCM$age <- as.numeric(dat_traits_DCM$age)
dat_traits_DCM[dat_traits_DCM =="Male"]<-0
dat_traits_DCM[dat_traits_DCM =="Female"]<-1
dat_traits_DCM[dat_traits_DCM =="DCM"]<-1
dat_traits_DCM[dat_traits_DCM =="NF"]<-0

dat_traits_DCM <- mutate_all(dat_traits_DCM, function(x) as.numeric(as.character(x)))

collectGarbage()

#-----------------------------------------------------------------------------#
## Step 3: Check dataset
  
# Cluster samples
sample_tree_DCM = hclust(dist(data.log_DCM), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
trait_colors_DCM = numbers2colors(dat_traits_DCM, signed = FALSE);
sizeGrWindow(12,12)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sample_tree_DCM, trait_colors_DCM,
                    groupLabels = names(dat_traits_DCM), cex.dendroLabels = 0.5, 
                    main = "Sample dendrogram and trait heatmap")

# ------------------------------------------------------------------------
## Step 4: Find soft threshold

# Choose a set of soft-thresholding powers
powers = seq(1,15, by=1)

# Call the network topology analysis function
sft_DCM = pickSoftThreshold(data.log_DCM, powerVector = powers, verbose = 5)
save(sft_DCM, file = "WGCNA-sft_Project.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

par(mar=c(1,1,1,1))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_DCM$fitIndices[,1], -sign(sft_DCM$fitIndices[,3])*sft_DCM$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft_DCM$fitIndices[,1], -sign(sft_DCM$fitIndices[,3])*sft_DCM$fitIndices[,2], labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft_DCM$fitIndices[,1], sft_DCM$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft_DCM$fitIndices[,1], sft_DCM$fitIndices[,5], labels=powers, cex=cex1,col="red")

# looking at both - soft threshold and mean connectivity 
# I decided to go with power 5  for this small example dataset
net_DCM = blockwiseModules(data.log_DCM, power = 5,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "expTOM", 
                       verbose = 3)

save(net_DCM, file = "WGCNA-net_project.RData")

#Let's visualize the modules dendrogram next.

# open a graphics window
sizeGrWindow(15, 9)

# Convert labels to colors for plotting
mergedColors_DCM = labels2colors(net_DCM$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_DCM$dendrograms[[1]], mergedColors_DCM[net_DCM$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels_DCM = net_DCM$colors
moduleColors_DCM = labels2colors(net_DCM$colors)
table(moduleColors_DCM)
MEs_DCM = net_DCM$MEs;
gene_tree_DCM = net_DCM$dendrograms[[1]]

#-----------------------------------------------------------------------------#
# Step 5: Find signficiant modules for trait data

# Define numbers of genes and samples
nGenes_DCM = ncol(data.log_DCM);
nSamples_DCM = nrow(data.log_DCM);

# Recalculate MEs with color labels
MEs0_DCM = moduleEigengenes(data.log_DCM, moduleColors_DCM)$eigengenes
MEs_DCM = orderMEs(MEs0_DCM)
moduleTraitCor_DCM = cor(MEs_DCM, dat_traits_DCM, use = "p");
moduleTraitPvalue_DCM = corPvalueStudent(moduleTraitCor_DCM, nSamples_DCM)

sizeGrWindow(20,20)

# Will display correlations and their p-values
textMatrix_DCM =  paste(signif(moduleTraitCor_DCM, 2), "\n(",
                    signif(moduleTraitPvalue_DCM, 1), ")", sep ="");
dim(textMatrix_DCM) = dim(moduleTraitCor_DCM)
par(mar = c(8, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
par(mar=c(1,1,1,1))
labeledHeatmap(Matrix = moduleTraitCor_DCM,
               xLabels = names(dat_traits_DCM),
               yLabels = names(MEs_DCM),
               ySymbols = names(MEs_DCM),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_DCM,
               setStdMargins = FALSE,
               cex.text = 0.9,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#-----------------------------------------------------------------------------#
#  Step 6: Get module membership and significance

# Define variable time containing the time column of datTrait
disease_DCM = as.data.frame(dat_traits_DCM$etiology);
names(disease_DCM) = "disease"

# names (colors) of the modules
modNames_DCM = substring(names(MEs_DCM), 3)

geneModuleMembership_DCM = as.data.frame(cor(data.log_DCM, MEs_DCM, use = "p"));
MMPvalue_DCM = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_DCM), nSamples_DCM));

names(geneModuleMembership_DCM) = paste("MM", modNames_DCM, sep="");
names(MMPvalue_DCM) = paste("p.MM", modNames_DCM, sep="");

geneTraitSignificance_DCM = as.data.frame(cor(data.log_DCM, disease_DCM, use = "p"));
GSPvalue_DCM = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_DCM), nSamples_DCM));

names(geneTraitSignificance_DCM) = paste("GS.", names(disease_DCM), sep="");

modules_DCM = c("black","blue","brown")
#modules_DCM = c("black","blue","brown","green","grey","magenta","pink","purple","red","turquoise","yellow")
#modules_DCM = c("black","blue","brown","green","grey","pink","red","turquoise","yellow")

sizeGrWindow(12, 4);
par(mfrow = c(1,3));
par(mar=c(1,1,1,1))
for(module in modules_DCM) {
  column_DCM = match(module, modNames_DCM)
  moduleGenes_DCM = moduleColors_DCM == module
  verboseScatterplot(abs(as.numeric(geneModuleMembership_DCM[moduleGenes_DCM, column_DCM])),
                     abs(as.numeric(geneTraitSignificance_DCM[moduleGenes_DCM, 1])),
                     xlab = paste("Module membership (MM,", module, ")", sep=""),
                     ylab = "Gene significance (GS) for disease",
                     main = paste("MM vs. GS\n", module, "Module"),
                     cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.2, col = module)
}

# Create the starting data frame
geneInfo0_DCM = data.frame(Gene.ID_DCM = colnames(data.log_DCM),
                       moduleColor_DCM = moduleColors_DCM,
                       geneTraitSignificance_DCM,
                       GSPvalue_DCM)
colnames(geneInfo0_DCM)[4] <- "GS.disease"

# Order modules by their significance for time
modOrder_DCM = order(-abs(cor(MEs_DCM, disease_DCM, use = "p")));

# Add module membership information in the chosen order
for (mod_DCM in 1 : ncol(geneModuleMembership_DCM))
{
  oldNames_DCM = names(geneInfo0_DCM)
  geneInfo0_DCM = data.frame(geneInfo0_DCM, geneModuleMembership_DCM[, modOrder_DCM[mod_DCM]], 
                         MMPvalue_DCM[, modOrder_DCM[mod_DCM]]);
  names(geneInfo0_DCM) = c(oldNames_DCM, paste("MM.", modNames_DCM[modOrder_DCM[mod_DCM]], sep=""),
                       paste("p.MM.", modNames_DCM[modOrder_DCM[mod_DCM]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_DCM = order(geneInfo0_DCM$moduleColor_DCM,-abs(geneInfo0_DCM$GS.disease));
geneInfo_DCM = geneInfo0_DCM[geneOrder_DCM, ]
write.csv(geneInfo_DCM, file = "geneInfo_project.csv", row.names = FALSE)

RCy3::cytoscapePing()
if (!grepl("status: Installed", RCy3::getAppStatus("stringApp"))) RCy3::installApp("stringApp")

#-----------------------------------------------------------------------------#
# Step 7: Recalculate topological overlap if needed - this takes quite a while, 
# I save the Rdata object and you can directly load it for this example

TOM_DCM = TOMsimilarityFromExpr(data.log_DCM, power = 5);
save(TOM_DCM, file = "WGCNA-TOM.RData")

# Select modules
modules = c("black");
#modules = c("brown");
#modules = c("grey");
#modules = c("red","black","brown");

# Select module probes
probes_DCM = colnames(data.log_DCM)
inModule_DCM = is.finite(match(moduleColors_DCM, modules));
modProbes_DCM = probes_DCM[inModule_DCM];

# Select the corresponding Topological Overlap
modTOM_DCM = TOM_DCM[inModule_DCM, inModule_DCM];
dimnames(modTOM_DCM) = list(modProbes_DCM, modProbes_DCM)

# Export the network into edge and node list files Cytoscape can read
cyt_DCM = exportNetworkToCytoscape(modTOM_DCM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.08,
                               nodeNames = modProbes_DCM,
                               nodeAttr = moduleColors_DCM[inModule_DCM])

edges_DCM <- read.table(paste0("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt"), sep="")
edges_DCM <- edges_DCM[-c(1),c(1:3)]
colnames(edges_DCM) <- c("source", "target", "weight")
write.table(edges_DCM, file = "edges_blue_DCM_power5.txt", sep = "\t", row.names = F, quote = F)

#import list of differentially expression genes between DCM and NF
setwd("D:/UM/Exp")
network_data <-  read.delim("DEGA_output_nw.txt", as.is = T);

#annotate sources and targets of network
gene_name_id = cbind(expression_DCM$uniprot_gn_symbol, rownames(expression_DCM))
colnames(gene_name_id) <- c("Gene_Name", "Gene_ID")

edges_DCM_source <- edges_DCM %>%
  left_join(as.data.frame(gene_name_id), by = c("source" = "Gene_ID"))

edges_DCM_target <- edges_DCM_source %>%
  left_join(as.data.frame(gene_name_id), by = c("target" = "Gene_ID")) 

edges_DCM_final <- edges_DCM_target[c(4,5,3)]
colnames(edges_DCM_final) <- c("source", "target", "weight")

#Connect to cytoscape
RCy3::createNetworkFromDataFrames(edges = edges_DCM_final, title=paste0("module", " ", paste(modules_DCM, collapse="-")))

#label node color for up and down genes using FC values
setwd("D:/UM/Exp/Project Period 3")
node_DCM <-  read.csv("node_DCM.csv", as.is = T)
node_DCM_FC <- node_DCM %>%
  left_join(DEGs_DCM, by = c("Gene" = "uniprot_gn_symbol")) 

write_xlsx(node_DCM_FC[,c(1,2)], path = "FC_DCM.xlsx")
write_xlsx(node_DCM_FC, path = "hub_module_DCM.xlsx")

#filter hub genes of both WGCNA and PPI 
black_MM = as.data.frame(cbind(geneModuleMembership_DCM[,8],rownames(geneModuleMembership_DCM)))

#filter hub genes with MM over than 0.8
# MM > 0.8
black_MM_2 <- black_MM[abs(as.numeric(black_MM[,1])) >= 0.8, ]
colnames(black_MM_2) = c("MM", "gene_id")
black_MM_2 = as.data.frame(black_MM_2)

GS_2 = cbind(geneTraitSignificance_DCM, rownames(geneTraitSignificance_DCM))
GS_3 = GS_2[abs(GS_2$GS.disease) >= 0.2,]
colnames(GS_3) = c("GS", "gene_id")
GS_3 = as.data.frame(GS_3)
hubgene_black<- as.data.frame(Reduce(intersect, list(black_MM_2$gene_id, GS_3$gene_id)))
colnames(hubgene_black) = "gene_id"

hubgene_black_name <- hubgene_black %>%
  left_join(as.data.frame(gene_name_id), by = c("gene_id" = "Gene_ID"))

# find overlaping gene between genes in hub module and target of DEmiRs
setwd("D:/UM/Exp/Project Period 3/Target")
up_target_DCM <-  read.csv("up_DEM_DCM.csv", as.is = T)
down_target_DCM <-  read.csv("down_DEM_DCM.csv", as.is = T)
targetofDCM = rbind(up_target_DCM, down_target_DCM)  
overlap.hub_DCM <- as.data.frame(Reduce(intersect, list(targetofDCM$DEMs, node_DCM$Gene)))

#-----------------------------------------------------------------------------#
# Part 2:miRNA-mRNA interaction network
#-----------------------------------------------------------------------------#

# merge overlapping targets miR for DCM
setwd("D:/UM/Exp/Project Period 3/Target")
up_target_DCM <-  read.csv("up_DEM_DCM.csv", as.is = T)
down_target_DCM <-  read.csv("down_DEM_DCM.csv", as.is = T)
targetofDCM = rbind(up_target_DCM, down_target_DCM)    

#load miR-mRNA network file
setwd("D:/UM/Exp/Project Period 3")
miRNA_DEGs_nw_DCM <- read.csv("miRNA_mRNA_nw_DCM.csv", as.is = T)
colnames(miRNA_DEGs_nw_DCM) = c("source", "target")

# create interaction network in cytoscape
RCy3::createNetworkFromDataFrames(edges = miRNA_DEGs_nw_DCM, title=paste0("module", " ", paste(modules, collapse="-")))

#-----------------------------------End----------------------------------------#

