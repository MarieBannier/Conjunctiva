## ======================================================================= ##
##                  integration of human conjunctiva datasets              ##
## ======================================================================= ##

# Load libraries
library(crayon,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(labeling,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(farver,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(withr,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(dplyr,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(SeuratObject,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(ggplot2,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(cowplot,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(patchwork, lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(usethis,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
install.packages("Seurat", lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(Seurat,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")

#library(devtools,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")

directory <- "/hpc/hub_clevers/MB/conjunctiva/Conjunctiva_10X/"
setwd(directory)


#########################################################################
# Load all datasets
#########################################################################
dir.create("QC_single_dataset")

LX237_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX237/an_392/outs/filtered_feature_bc_matrix/")
LX237 <- CreateSeuratObject(counts = LX237_before, project = "tissue", min.cells = 3, min.features = 200)
rm(LX237_before)
LX237[["percent.mt"]] <- PercentageFeatureSet(LX237, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset1.csv")

LX243_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX243/an_394/outs/filtered_feature_bc_matrix/")
LX243 <- CreateSeuratObject(counts = LX243_before, project = "ALId17_hCjM04", min.cells = 3, min.features = 200)
rm(LX243_before)
LX243[["percent.mt"]] <- PercentageFeatureSet(LX243, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset2.csv")

LX244_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX244/an_395/outs/filtered_feature_bc_matrix/")
LX244 <- CreateSeuratObject(counts = LX244_before, project = "ALId17_hCjM16", min.cells = 3, min.features = 200)
rm(LX244_before)
LX244[["percent.mt"]] <- PercentageFeatureSet(LX244, pattern = "^MT-")

TX294_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX294/count/sample_filtered_feature_bc_matrix/")
TX294 <- CreateSeuratObject(counts = TX294_before$`Gene Expression`, project = "Organoids_EM", min.cells = 3, min.features = 200)
rm(TX294_before)
TX294[["percent.mt"]] <- PercentageFeatureSet(TX294, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset4.csv")

TX295_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX295/count/sample_filtered_feature_bc_matrix/")
TX295 <- CreateSeuratObject(counts = TX295_before$`Gene Expression`, project = "Organoids_DM", min.cells = 3, min.features = 200)
rm(TX295_before)
TX295[["percent.mt"]] <- PercentageFeatureSet(TX295, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset5.csv")

TX296_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX296/count/sample_filtered_feature_bc_matrix/")
TX296 <- CreateSeuratObject(counts = TX296_before$`Gene Expression`, project = "ALId0", min.cells = 3, min.features = 200)
rm(TX296_before)
TX296[["percent.mt"]] <- PercentageFeatureSet(TX296, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset6.csv")

TX297_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX297/count/sample_filtered_feature_bc_matrix/")
TX297 <- CreateSeuratObject(counts = TX297_before$`Gene Expression`, project = "ALId3", min.cells = 3, min.features = 200)
rm(TX297_before)
TX297[["percent.mt"]] <- PercentageFeatureSet(TX297, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset7.csv")

TX298_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX298/count/sample_filtered_feature_bc_matrix/")
TX298 <- CreateSeuratObject(counts = TX298_before$`Gene Expression`, project = "ALId17_IL4+13", min.cells = 3, min.features = 200)
rm(TX298_before)
TX298[["percent.mt"]] <- PercentageFeatureSet(TX298, pattern = "^MT-")
write.csv(c(0,0,1,2), "dataset8.csv")

#save.image("envir_v1.Rdata")

#########################################################################
# QC on individual datasets
#########################################################################
dir.create("QC_single_dataset")

pdf("QC_single_dataset/QC.pdf")
VlnPlot(LX237, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("LX237")
VlnPlot(LX243, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("LX243")  
VlnPlot(LX244, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("LX244")
VlnPlot(TX294, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("TX294")
VlnPlot(TX295, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("TX295")
VlnPlot(TX296, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("TX296")
VlnPlot(TX297, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("TX297")
VlnPlot(TX298, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("TX298")
dev.off()

#########################################################################
# Integrate datasets
#########################################################################
   
dataset <- merge(LX237, y = c(LX243, LX244, TX294, TX295, TX296, TX297, TX298),
                 add.cell.ids = c("Tissue", "ALId17_hCjM04", "ALId17_hCjM16",
                                "Organoids_EM", "Organoids_DM", "ALId0", "ALId3",
                                "ALId17_IL4+13"),
                 project = "conjunctiva")
write.csv(c(0,0,1,2), "datasetmerged.csv")
                 
dataset.list <- SplitObject(dataset, split.by = "orig.ident")
write.csv(c(0,0,1,2), "datasetsplitted.csv")

# normalize and identify variable features for each dataset independently
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
    x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
write.csv(c(0,0,1,2), "datasetnormalized.csv")

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(dataset.list, nfeatures = 2000)

dataset.list <- lapply(X = dataset.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
write.csv(c(0,0,1,2), "datasetscaled.csv")

anchors <- FindIntegrationAnchors(object.list = dataset.list, anchor.features = features,
                                  reduction = "rpca")
write.csv(c(0,0,1,2), "datasetanchored.csv")

dataset_combined <- IntegrateData(anchorset = anchors)
write.csv(c(0,0,1,2), "datasetcombined.csv")
rm(dataset.list)
rm(anchors)
rm(features)
rm(LX243)
rm(LX244)
rm(TX294)
rm(TX295)
rm(TX296)
rm(TX297)
rm(TX298)
rm(dataset)

# Analyse the integrated data
DefaultAssay(dataset_combined) <- "integrated"

#########################################################################
# Perform QC on integrated dataset
#########################################################################
dir.create("QC_merged_dataset")

pdf("QC_merged_dataset/QC.pdf")
VlnPlot(dataset_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts')
VlnPlot(dataset_combined, group.by = "orig.ident", 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0.1, ncol = 3) +
    NoLegend()
dev.off()

plot1 <- FeatureScatter(dataset_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dataset_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("QC_merged_dataset/Scatter-counts.pdf", width = 14)
plot1 + plot2
dev.off()

#########################################################################
# Save merged environment
#########################################################################
save.image("envir_with_variables_v1.RData")

#########################################################################
# Visualization and Clustering
#########################################################################
dir.create("Analysis")

## Dimension 10
dir.create("Analysis/Resolution_0,5_Dims_10")

# Analysis
dataset_combined <- ScaleData(dataset_combined, verbose = FALSE)
dataset_combined <- RunPCA(dataset_combined, npcs = 30, verbose = FALSE)
dataset_combined <- RunUMAP(dataset_combined, reduction = "pca", dims = 1:10)
dataset_combined <- FindNeighbors(dataset_combined, reduction = "pca", dims = 1:10)
dataset_combined <- FindClusters(dataset_combined, resolution = 0.5)

# Elbow plot
pdf("Analysis/Resolution_0,5_Dims_10/elbow_plot.pdf")
ElbowPlot(dataset_combined)
dev.off()

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_10/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_10/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()

## Dimension 9
dir.create("Analysis/Resolution_0,5_Dims_9")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:9)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_9/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_9/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 8
dir.create("Analysis/Resolution_0,5_Dims_8")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:8)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_8/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_8/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 7
dir.create("Analysis/Resolution_0,5_Dims_7")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:7)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_7/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_7/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 6
dir.create("Analysis/Resolution_0,5_Dims_6")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:6)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_6/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_6/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 6
dir.create("Analysis/Resolution_0,5_Dims_6")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:6)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_6/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_6/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 5
dir.create("Analysis/Resolution_0,5_Dims_5")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:5)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_5/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_5/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 4
dir.create("Analysis/Resolution_0,5_Dims_4")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:4)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_4/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_4/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 11
dir.create("Analysis/Resolution_0,5_Dims_11")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:11)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_11/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_11/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 12
dir.create("Analysis/Resolution_0,5_Dims_12")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:12)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_12/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_12/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 13
dir.create("Analysis/Resolution_0,5_Dims_13")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:13)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_13/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_13/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 14
dir.create("Analysis/Resolution_0,5_Dims_14")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:14)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_14/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_14/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


## Dimension 15
dir.create("Analysis/Resolution_0,5_Dims_15")

# Analysis
dataset_combined <- RunUMAP(dataset_combined, dims = 1:15)

# Visualization
p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_15/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_15/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()


#########################################################################
# Save merged environment
#########################################################################
save.image("envir_with_Resolutions_v1.RData")



