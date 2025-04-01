#=============================================================================#
# Research Project 1
# Section 4a: Network analysis: WGCNA, PPI network, mRNA-miRNA interaction network - HCM vs NF
#																	                                       		   															                                #
# Date: January 24, 2025											                                
# Author: Mai Le, ID: i6375777
# Maastricht University                                                       #
#=============================================================================#

#-----------------------------------------------------------------------------#
# Part 1: WGCNA and PPI network for DEGs between HCM and NF
#-----------------------------------------------------------------------------#
setwd("D:/UM/Exp")
expression <-  read.delim("annotated_expression.txt", as.is = T, 
                          row.names = 1);
sampleInfo <-  read.csv("MAGNET_SampleData_18112022.csv", as.is = T, 
                        row.names = 1);

#load list of DEGs
setwd("D:/UM/Exp/Project Period 3/DEGs")
DEGs <- read.delim("DEGs_HCM.txt", as.is = T);

expression_2 <- expression[DEGs$ensembl_gene_id,]
expression_2$uniprot_gn_symbol = DEGs$uniprot_gn_symbol

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
data_HCM <- expression_2[, sampleInfo$etiology == "HCM"]
data_NF <- expression_2[, sampleInfo$etiology == "NF"]
data_NF = data_NF[,-c(167,168)]
data.log <- rbind(t(data_HCM), t(data_NF))

#extract metadata 
traitdata <- sampleInfo[rownames(data.log),c(2,3,5)]
traitdata$sampleID <- rownames(traitdata)

#-----------------------------------------------------------------------------#
# Step 2: Form a data frame analogous to expression data that will hold the clinical traits.

sample <- rownames(data.log)
traitrows = match(sample, traitdata$sampleID);
dat_traits = traitdata[traitrows, -4];
rownames(dat_traits) = traitrows

dat_traits$age <- as.numeric(dat_traits$age)
dat_traits[dat_traits=="Male"]<-0
dat_traits[dat_traits=="Female"]<-1
dat_traits[dat_traits=="HCM"]<-1
dat_traits[dat_traits=="NF"]<-0

dat_traits <- mutate_all(dat_traits, function(x) as.numeric(as.character(x)))

collectGarbage()

#-----------------------------------------------------------------------------#
## Step 3: Check dataset

# Cluster samples
  
# Cluster samples
sample_tree = hclust(dist(data.log), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
trait_colors = numbers2colors(dat_traits, signed = FALSE);
sizeGrWindow(12,12)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sample_tree, trait_colors,
                    groupLabels = names(dat_traits), cex.dendroLabels = 0.5, 
                    main = "Sample dendrogram and trait heatmap")

# ------------------------------------------------------------------------
## Step 4: Find soft threshold

# Choose a set of soft-thresholding powers
powers = seq(1,15, by=1)

# Call the network topology analysis function
sft = pickSoftThreshold(data.log, powerVector = powers, verbose = 5)

save(sft, file = "WGCNA-sft_Project.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

par(mar=c(1,1,1,1))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# looking at both - soft threshold and mean connectivity 
# I decided to go with power 6 for this small example dataset
net = blockwiseModules(data.log, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "expTOM", 
                       verbose = 3)

save(net, file = "WGCNA-net_project.RData")

#Let's visualize the modules dendrogram next.

# open a graphics window
sizeGrWindow(15, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;
gene_tree = net$dendrograms[[1]]

#-----------------------------------------------------------------------------#
# Step 5: Find signficiant modules for trait data

# Define numbers of genes and samples
nGenes = ncol(data.log);
nSamples = nrow(data.log);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(data.log, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, dat_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(20,20)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep ="");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
par(mar=c(1,1,1,1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dat_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#-----------------------------------------------------------------------------#
#  Step 6: Get module membership and significance

# Define variable time containing the time column of datTrait
disease = as.data.frame(dat_traits$etiology);
names(disease) = "disease"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data.log, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(data.log, disease, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(disease), sep="");

modules = c("grey","yellow","blue")
#modules = c("blue","brown","green","grey","red","turquoise","yellow")

sizeGrWindow(9, 3);
par(mfrow = c(1,3));
par(mar=c(1,1,1,1))
for(module in modules) {
  column = match(module, modNames)
  moduleGenes = moduleColors == module
  
  verboseScatterplot(abs(as.numeric(geneModuleMembership[moduleGenes, column])),
                     abs(as.numeric(geneTraitSignificance[moduleGenes, 1])),
                     xlab = paste("Module membership (MM,", module, ")", sep=""),
                     ylab = "Gene significance (GS) for disease",
                     main = paste("MM vs. GS\n", module, "Module"),
                     cex.main = 1.2, cex.lab = 1.5, cex.axis = 1.2, col = module)
}

# Create the starting data frame
geneInfo0 = data.frame(Gene.ID = colnames(data.log),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for time
modOrder = order(-abs(cor(MEs, disease, use = "p")));

# Add module membership information in the chosen order
for (mod in 1 : ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.disease));
geneInfo = geneInfo0[geneOrder, ]

RCy3::cytoscapePing()
if (!grepl("status: Installed", RCy3::getAppStatus("stringApp"))) RCy3::installApp("stringApp")

#-----------------------------------------------------------------------------#
# Step 7: Recalculate topological overlap if needed - this takes quite a while, 
# I save the Rdata object and you can directly load it for this example

TOM = TOMsimilarityFromExpr(data.log, power = 6);
save(TOM, file = "WGCNA-TOM.RData")

# Select modules
modules = c("blue");
#modules = c("grey");
#modules = c("blue","yellow");

# Select module probes
probes = colnames(data.log)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.08,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])

edges <- read.table(paste0("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt"), sep="")
edges <- edges[-c(1),c(1:3)]
colnames(edges) <- c("source", "target", "weight")
write.table(edges, file = "edges_yellow_blue_grey.txt", sep = "\t", row.names = F, quote = F)

#import list of differentially expression genes between DCM and NF
setwd("D:/UM/Exp")
network_HCM <-  read.delim("DEGA_output_nw.txt", as.is = T);

#annotate source and target
gene_name_id_HCM = cbind(expression_2$uniprot_gn_symbol, rownames(expression_2))
colnames(gene_name_id_HCM) <- c("Gene_Name", "Gene_ID")

edges_HCM_source <- edges %>%
  left_join(as.data.frame(gene_name_id_HCM), by = c("source" = "Gene_ID"))

edges_HCM_target <- edges_HCM_source %>%
  left_join(as.data.frame(gene_name_id_HCM), by = c("target" = "Gene_ID")) 

edges_HCM_final <- edges_HCM_target[c(4,5,3)]
colnames(edges_HCM_final) <- c("source", "target", "weight")

#Connect to cytoscape
RCy3::createNetworkFromDataFrames(edges = edges_HCM_final, title=paste0("module", " ", paste(modules, collapse="-")))

#label node color for up and down genes using FC values
setwd("D:/UM/Exp/Project Period 3")
node_HCM <-  read.csv("node_HCM.csv", as.is = T)
node_HCM_FC <- node_HCM %>%
  left_join(DEGs, by = c("Gene" = "uniprot_gn_symbol")) 

write_xlsx(node_HCM_FC[,c(1,2)], path = "hub_module_HCM.xlsx")
write_xlsx(node_HCM_FC[,c(1,2)], path = "hub_module_HCM.xlsx")
write_xlsx(node_HCM_FC, path = "hub_module_HCM.xlsx")

#filter hub genes of both WGCNA and PPI 
blue_MM = as.data.frame(cbind(geneModuleMembership[,5],rownames(geneModuleMembership)))

#filter hub genes with MM over than 0.8
# MM > 0.8
blue_MM_2 <- blue_MM[abs(as.numeric(blue_MM[,1])) >= 0.8, ]
colnames(blue_MM_2) = c("MM", "gene_id")
blue_MM_2 = as.data.frame(blue_MM_2)

#filter hub genes with geneTraitSignificance > 0.2
GS_2 = cbind(geneTraitSignificance, rownames(geneTraitSignificance))
GS_3 = GS_2[abs(GS_2$GS.disease) >= 0.2,]
colnames(GS_3) = c("GS", "gene_id")
GS_3 = as.data.frame(GS_3)

hubgene_blue<- as.data.frame(Reduce(intersect, list(blue_MM_2$gene_id, GS_3$gene_id)))
colnames(hubgene_blue) = "gene_id"
hubgene_blue_name <- hubgene_blue %>%
  left_join(as.data.frame(gene_name_id_HCM), by = c("gene_id" = "Gene_ID"))

# find overlaping gene between genes in hub module and target of DEmiRs
setwd("D:/UM/Exp/Project Period 3/Target")
up_target <-  read.csv("up_DEMs_HCM.csv", as.is = T)
down_target <-  read.csv("down_DEMs_HCM.csv", as.is = T)
target = rbind(up_target, down_target)                     
overlap.hub_HCM <- as.data.frame(Reduce(intersect, list(target$DEMs, node_HCM$Gene)))

#-----------------------------------------------------------------------------#
# Part 2:miRNA-mRNA interaction network
#-----------------------------------------------------------------------------#

# merge overlapping targets miR for DCM
setwd("D:/UM/Exp/Project Period 3/Target")
up_target_HCM <-  read.csv("up_DEM_HCM.csv", as.is = T)
down_target_HCM <-  read.csv("down_DEM_HCM.csv", as.is = T)
targetofHCM = rbind(up_target_HCM, down_target_HCM)    

#load miR-mRNA network file
setwd("D:/UM/Exp/Project Period 3")
miRNA_DEGs_nw_HCM <- read.csv("miRNA_mRNA_nw_HCM.csv", as.is = T)
colnames(miRNA_DEGs_nw_HCM) = c("source", "target")

# create interaction network in cytoscape
RCy3::createNetworkFromDataFrames(edges = miRNA_DEGs_nw_HCM, title=paste0("module", " ", paste(modules, collapse="-")))

#-----------------------------------End----------------------------------------#
