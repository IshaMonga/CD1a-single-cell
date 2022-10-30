#########################Staphylococcal phosphatidylglycerol antigens activate human T cells via CD1a: Monnot et.al.,###############
##########Authors_script: Monnot and Monga et.al.,##########
##############Stepwise reproducible Code for the scRNAseq analysis 


## set working directory
setwd("~/Documents/CD1a-work/Processed_files")
projDir<-("~/Documents/CD1a-work/Processed_files/scRNAseq_outputs_harmony_batch_2022/")
##cleans global environment
rm(list=ls())
##installs necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repo = "https://mac.R-project.org")
BiocManager::install()
BiocManager::install(c("SingleCellExperiment"), site_repository = "https://mac.R-project.org")
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("svglite")
install.packages("spatstat.sparse")

if (!require("here")) {
  install.packages("here", dependencies = TRUE)
  library(here)
}

if (!require("useful")) {
  install.packages("useful", dependencies = TRUE)
  library(useful)
}
if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
if (!require("Seurat")) {
  install_version('Seurat', version = '2.3.2')
  library(Seurat)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("Matrix")) {
  install.packages("Matrix", dependencies = TRUE)
  library(Matrix)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("openxlsx")) {
  install.packages("openxlsx", dependencies = TRUE)
  library(openxlsx)
}
if (!require("tidyr")) {
  install.packages("tidyr", dependencies = TRUE)
  library(tidyr)
}
if (!require("plotly")) {
  install.packages("plotly", dependencies = TRUE)
  library(plotly)
}
if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("htmlwidgets")) {
  install.packages("htmlwidgets", dependencies = TRUE)
  library(htmlwidgets)
}
if (!require("htmltools")) {
  install.packages("htmltools", dependencies = TRUE)
  library(htmltools)
}
if (!require("ggrepel")) {
  install.packages("ggrepel", dependencies = TRUE)
  library(ggrepel)
}
if (!require("gridExtra")) {
  install.packages("gridExtra", dependencies = TRUE)
  library(gridExtra)
}
if (!require("varSelRF")) {
  install.packages("varSelRF", dependencies = TRUE)
  library(varSelRF)
}
if (!require("biomaRt")) {
  BiocManager::install('biomaRt')
  library(biomaRt)
}
if (!require("matrixStats")) {
  install.packages("matrixStats", dependencies = TRUE)
  library(matrixStats)
}
if (!require("quantmod")) {
  install.packages("quantmod", dependencies = TRUE)
  library(quantmod)
}
if (!require("GGally")) {
  install.packages("GGally", dependencies = TRUE)
  library(GGally)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}



#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages('Seurat')


##loads necessary packages
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(cowplot)
library(RColorBrewer)
library(scales)

#READ DATA ######################################################

##reads the count matrix of 921 and creates a data.frame "count_df"
counts_df <- read.csv(file = "./umi_counts_921_2.csv")
head(counts_df)
##adds gene names as row names
row.names(counts_df)<- counts_df$X
##remove first row which is non numeric and has row names
counts_df<- counts_df[,-(1)]
counts <- as.matrix(counts_df)

##reads the meta data matrix and creates a data.frame "cell_info"
cell_info_921 <- read.csv(file="metadata_921_2_new.csv", header = TRUE)
head(cell_info_921,5)
row.names(cell_info_921)<- cell_info_921$X
cell_info_921<- cell_info_921[,-(1)]
mat_cell_info_921 <- as.matrix(cell_info_921)

##Create seurat object with 921
seurat_921 <- CreateSeuratObject(counts = counts, meta.data = cell_info_921, min.cells = 3, min.features = 200)
seurat_921


##reads the count matrix of 325 and creates a data.frame "count_df"
counts_df <- read.csv(file = "all_plates_umi_counts_325.csv")
##adds gene names as row names
row.names(counts_df)<- counts_df$X
##remove first row which is non numeric and has row names
counts_df<- counts_df[,-(1)]
counts <- as.matrix(counts_df)

##reads the meta data matrix and creates a data.frame "cell_info"
cell_info_325 <- read.csv(file="metadata_325_new.csv", header = TRUE)
head(cell_info_325,10)
row.names(cell_info_325)<- cell_info_325$X
cell_info_325<- cell_info_325[,-(1)]
mat_cell_info_325 <- as.matrix(cell_info_325)

##Create seurat object with 325
seurat_325 <- CreateSeuratObject(counts = counts, meta.data = cell_info_325, min.cells = 3, min.features = 200)
seurat_325


###############################NORMALIZationa nd mito genes addition###########################################
#1 seurat
##Add percentage of mitochondrial genes to meta data of 1st Seurat object
seurat_921[["percent.mt"]] <- PercentageFeatureSet(seurat_921, pattern = "^MT-")

##To check the metadata of the first 5 cells
head(seurat_921@meta.data, 5)

##Visualize the features of the "outliers" on a plot
plot1 <- FeatureScatter(seurat_921, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'cell_type')
plot1

##remove the cells that are outside of the "norm"
seurat_921 <- subset(seurat_921, subset = nCount_RNA > 200 & nCount_RNA < 3200 & percent.mt < 12)

##check again
plot1 <- FeatureScatter(seurat_921, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'cell_type')
plot1

##Normalize data
seurat_921 <- NormalizeData(seurat_921)

##Identify most variable genes
seurat_921 <- FindVariableFeatures(seurat_921, selection.method = "vst", nfeatures = 2000)

#2nd seurat
##Add percentage of mitochondrial genes to meta data of second Seurat object
seurat_325[["percent.mt"]] <- PercentageFeatureSet(seurat_325, pattern = "^MT-")

##To check the metadata of the first 5 cells
head(seurat_325@meta.data, 5)

##Visualize the features of the "outliers" on a plot
plot2 <- FeatureScatter(seurat_325, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'cell_type')
plot2

##remove the cells that are outside of the "norm"
seurat_325 <- subset(seurat_325, subset = nCount_RNA > 200 & nCount_RNA < 5000 & percent.mt < 11)

##check
plot2 <- FeatureScatter(seurat_325, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'cell_type')
plot2

##Normalize data
seurat_325 <- NormalizeData(seurat_325)

##Identify most variable genes
seurat_325 <- FindVariableFeatures(seurat_325, selection.method = "vst", nfeatures = 2000)

#MERGE ###################################

##create list of seurat objects
seurat.list <- list("donor_325" = seurat_325, "donor_921" = seurat_921)

##create anchors and merge datasets into one seurat object
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)
seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
DefaultAssay(seurat.integrated) <- "integrated"

seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)
table(seurat.integrated$sample)
#library(harmony)
#seurat.integrated <- RunHarmony(seurat.integrated, "sample",assay.use = "integrated")
#seurat.integrated <- RunUMAP(seurat.integrated, reduction = "harmony", dims = 1:30)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30)

p1 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "donor")
p1
p2 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "cell_type")
p2
plot(p1+p2)
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.9)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30)

p3 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters")
p3

##calculate most differentially expressed markers by clusters and write csv file
markers_integrated <- FindAllMarkers(seurat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(markers_integrated, paste0(projDir, 'merged_all_clusters_diff_expr_all_2022.csv'))

##Finds the top10 genes and 
#top_genes = 20
#top_genes = 10
top10<-sub_markers_integrated <- markers_integrated %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#top_deg = degClusters %>% group_by(cluster) %>% top_n (top_genes, avg_log2FC)
#top10 <- markers_integrated %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(seurat.integrated, features = top10$gene) 

##Dotplot##
features<-c("IL2", "IL3", "IL4", "IL5", "IL6", "IL10", "IL13", "IL21", "IL26",  "IFNG", "TNF", "CSF1", "CSF2", "OSM", "GZMB","GZMH", "PRF1")
DotPlot(seurat.integrated, features = features) + RotatedAxis()
#features<-seurat.integrated@assays$RNA


feat_p = FeaturePlot (seurat.integrated, features = unique(top10$gene), keep.scale = 'all', order=F,
                      combine = FALSE, pt.size = .01, col=viridis::turbo(100))
for(i in 1:length(feat_p)) {
  print(length(feat_p))
  feat_p[[i]] = feat_p[[i]] + theme_void() + NoLegend() + NoAxes()
  # print(feat_p[[i]])
  
}
getwd()
pdf("./scRNAseq_outputs_harmony_batch_2022/cluster_deg_fplots.pdf", useDingbats = F, width = 22, height = 10)
print (wrap_plots (feat_p))
dev.off()

######STEPS TO REMOVE THE CLUSTER 7 AND CD1A MOCK#####
##To remove cluster 7
seurat.subset <- subset(seurat.integrated, subset = seurat_clusters != 7)
abc<-table(seurat.subset$seurat_clusters)
#0   1   2   3   4   5   6   7 
#350 226 129 111  80  42  25   0 

write.csv(abc, paste0(projDir,'Cells_numbers_without_cluster7.csv'))
##To remove CD1a mock
seurat.subset <- subset(seurat.subset, subset = cell_type != "CD1a mock")
seurat.subset <- FindNeighbors(seurat.subset, dims = 1:30)
seurat.subset <- FindClusters(seurat.subset, resolution = 0.9)
abc<-table(seurat.subset$seurat_clusters)

write.csv(abc, paste0(projDir, 'Cells_numbers_without_cluster7_and_CD1a_mock.csv'))
abc1<-table(seurat.subset$seurat_clusters,seurat.subset$cell_type)

write.csv(abc1,paste0(projDir,'Cells_numbers_per_cluster_by_celltype_cd1aLPG_Tetramer.csv'))


seurat.subset <- RunUMAP(seurat.subset, dims = 1:30)
DimPlot(seurat.subset, reduction = "umap", pt.size = 0.5)
p1<-DimPlot(seurat.subset, reduction = "umap", group.by = "cell_type", pt.size = 0.5)
p2<-DimPlot(seurat.subset, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5)
p3<-DimPlot(seurat.subset, reduction = "umap", group.by = "sample", pt.size = 0.5)
pdf("./scRNAseq_outputs_harmony_batch_2022/UMAPs_subset_noCD1a.mock.cluster7.pdf", useDingbats = F, width = 22, height = 10)
plot(p1+p2+p3)
dev.off()

new.markers_integrated <- FindAllMarkers(seurat.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

##Finds the top10 genes and 
top10 <- new.markers_integrated %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- new.markers_integrated %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("./scRNAseq_outputs_harmony_batch_2022/doheatmap_seurat_subset_nocd1amock_clusetr7.pdf",useDingbats = F, width = 22, height = 19)
DoHeatmap(seurat.subset, features = top10$gene)
dev.off()

pdf("./scRNAseq_outputs_harmony_batch_2022/dotplot_seurat_subset_nocd1amock_clusetr7.pdf",useDingbats = F, width = 22, height = 19)
DotPlot(seurat.subset, features = features, group.by= "seurat_clusters") + RotatedAxis()
dev.off()
#VlnPlot(seurat.integrated, features = c( "TRAC", "TRBC1", "CD3D"), group.by = "seurat_clusters")

#FeaturePlot(seurat.integrated, features = c("IL13", "IFNG"))

###write CSV new markers
write.csv(new.markers_integrated, file="./scRNAseq_outputs_harmony_batch_2022/seurat_subset_clusters_diff_expr_all.csv")

# ##Name Clusters
# clusters.ids <- c("Resting cells", "Resting?", "Th2", "?", "Tregs", "Cycling cells", 
#                   "TEMRA", "contaminant")
# names(clusters.ids) <- levels(seurat.integrated)
# seurat.integrated <- RenameIdents(seurat.integrated, clusters.ids)
# DimPlot(seurat.integrated, reduction = "umap", label = TRUE, pt.size = 0.5)

##clusters.ids <- c("Resting cells", "Resting?", "Th2", "?", "Tregs", "Cycling cells", 
##                                     "TEMRA")
##names(clusters.ids) <- levels(seurat.subset)
##seurat.subset <- RenameIdents(seurat.subset, clusters.ids)
#p1 <- DimPlot(seurat.subset, reduction = "umap", pt.size = 0.5)
p1<-DimPlot(seurat.subset, reduction = "umap", group.by = "cell_type", pt.size = 0.5)
p2<-DimPlot(seurat.subset, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5)
p3<-DimPlot(seurat.subset, reduction = "umap", group.by = "sample", pt.size = 0.5)
pdf("./scRNAseq_outputs_harmony_batch_2022/New_UMAPs_subset_nocd1a.mock.cluster7.pdf", useDingbats = F, width = 22, height = 10)
plot(p1+p2+p3)
dev.off()

table.clusters <- table(seurat.subset@active.ident, seurat.subset@meta.data$cell_type)
write.csv(table.clusters, file="./scRNAseq_outputs_harmony_batch_2022/table_clusters.csv")

#Hierarchical clustering ##########

##finds markers for each cluster and takes top 10
DP.markers <- FindMarkers(seurat.subset, ident.1="CD1a double pos", group.by = "cell_type")
genes <- rownames(DP.markers)
DP.markers <- cbind(genes, DP.markers)
DP.top10 <- DP.markers %>% top_n(n = 10, wt = avg_logFC)

LPG.markers <- FindMarkers(seurat.subset, ident.1="CD1a LPG+", group.by = "cell_type")
genes <- rownames(LPG.markers)
LPG.markers <- cbind(genes, LPG.markers)
LPG.top10 <- LPG.markers %>% top_n(n = 10, wt = avg_logFC)

mock.markers <- FindMarkers(seurat.subset, ident.1="CD1a mock", group.by = "cell_type")
genes <- rownames(mock.markers)
mock.markers <- cbind(genes, mock.markers)
mock.top10 <- mock.markers %>% top_n(n = 10, wt = avg_logFC)

neg.markers <- FindMarkers(seurat.subset, ident.1="Tetramer neg", group.by = "cell_type")
genes <- rownames(neg.markers)
neg.markers <- cbind(genes, neg.markers)
neg.top10 <- neg.markers %>% top_n(n = 10, wt = avg_logFC)

write.csv(DP.markers, file="DPmarkers.csv")
write.csv(DP.top10, file="DPmarkers_top10.csv")
write.csv(LPG.markers, file="LPGmarkers.csv")
write.csv(LPG.top10, file="LPGmarkers_top10.csv")
write.csv(mock.markers, file="mockmarkers.csv")
write.csv(mock.top10, file="mockmarkers_top10.csv")
write.csv(neg.markers, file="negmarkers.csv")
ewrite.csv(neg.top10, file="negmarkers_top10.csv")

#FIGURES FOR PAPER###########
pdf("./scRNAseq_outputs_harmony_batch_2022/Dimplot_by_Donor_nomock_clus7.pdf", useDingbats = F, width = 22, height = 10)
p3 <- DimPlot(seurat.subset, reduction = "umap", group.by = "donor", pt.size = 0.7)
p3
dev.off()
pdf("./scRNAseq_outputs_harmony_batch_2022/Dimplot_by_seuClusters_nomock_clus7.pdf", useDingbats = F, width = 22, height = 10)
p4 <- DimPlot(seurat.subset, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.7)
p4
dev.off()
##p4bis <- DimPlot(seurat.subset, reduction = "umap", group.by = "cell_type", pt.size = 0.7, 
             ## cols= c("CD1a double pos"= "Salmon" , "CD1a LPG+"="BlueViolet", "CD1a mock"="DeepSkyBlue", "Tetramer neg" = "grey"))
##p4bis
##clusters.ids <- c("Resting cells", "Resting?", "Th2", "?", "Tregs", "Cycling cells", 
##                  "TEMRA")
##names(clusters.ids) <- levels(seurat.subset)
##seurat.subset <- RenameIdents(seurat.subset, clusters.ids)
pdf("./scRNAseq_outputs_harmony_batch_2022/Dimplot_by_celltype.pdf", useDingbats = F, width = 22, height = 10)
p5 <- DimPlot(seurat.subset, reduction = "umap", pt.size = 0.7, group.by = "cell_type", cols= c("CD1a LPG+"="BlueViolet", "Tetramer neg" = "sandybrown"))
p5
dev.off()
cluster.averages <- AverageExpression(seurat.subset, return.seurat = TRUE, add.ident = "donor")
pdf("./scRNAseq_outputs_harmony_batch_2022/doHeatmap_avgexpression_clusters.pdf", useDingbats = F, width = 22, height = 10)
p6 <-DoHeatmap(cluster.averages, features = top5$gene, size = 3, draw.lines = FALSE)
p6
dev.off()
p7 <- DoHeatmap(seurat.subset, features = top10$gene)
p7
#p7 <-DoHeatmap(cluster.averages, features = top10$gene, size = 3, draw.lines = FALSE)
#p7

##VIOLIN PLOTS#######
##TH2
p8 <-VlnPlot(seurat.subset, features = c( "IL13", "IL2", "IL4", "CCR4", "FASLG", "CD40LG"), assay = "RNA", group.by = "seurat_clusters", pt.size = 0)
##Treg
p9 <- VlnPlot(seurat.subset, features = c( "CTLA4", "PDCD1", "TNFRSF4", "GATA3", "TIGIT", "FOXP3"), assay = "RNA", group.by = "seurat_clusters", pt.size = 0)
##CTX
p10 <- VlnPlot(seurat.subset, features = c( "CCL3","CCL4", "XCL1", "IFNG", "GZMB", "PRF1"), assay = "RNA", group.by = "seurat_clusters", pt.size = 0)

p8
p9
p10

##extraviolinplots
p11 <- VlnPlot(seurat.subset, features = c( "GNLY", "KLRD1", "LRRC32", "HLA-DRB1", "IL2RA", "TNFRSF18"), assay = "RNA", group.by = "seurat_clusters", pt.size = 0)
p12 <- VlnPlot(seurat.subset, features = c( "IL22", "LTA", "LIF", "CSF2", "STAT6", "IL5"), assay ="RNA", group.by = "seurat_clusters", pt.size = 0)
p13 <- VlnPlot(seurat.subset, features = c( "NKG7", "CRTAM", "IKZF2", "SLAMF7", "CISH", "IL13"), assay = "RNA", group.by = "seurat_clusters", pt.size = 0)


# Open a pdf file
pdf("rplot20210426.pdf") 
# 2. Create a plot
p3; p5; p4; p7; p8; p9; p10; p11; p12; p13
# Close the pdf file
dev.off()


saveRDS(seurat.subset, paste0(projDir, 'Main_seurat_analyses_cd1a_tet_positive.rds'))