## ======================================================================= ##
##                  integration of human conjunctiva datasets              ##
## ======================================================================= ##

# Load libraries
library(crayon)
library(labeling)
library(farver)
library(withr)
library(dplyr)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(patchwork)
library(usethis)
library(Seurat)
library(hdf5r)


#########################################################################
# Load gut dataset (normalized from https://www.gutcellatlas.org/)
#########################################################################

# ---- Load converted gut data from .h5ad to several .csv ---- #
message("loading data")
genenames <- read.csv("for_Marie/var.csv", header = T, row.names = 1)
metadata <- read.csv("for_Marie/obs.csv", header = T, row.names = 1)
data <- read.csv("for_Marie/X.csv", header = F)
rownames(data) <- rownames(metadata)
colnames(data) <- rownames(genenames)

data <- t(data)

# ---- Create Seurat object---- #
dataset_combined <- CreateSeuratObject(counts = data, project = "gut", min.cells = 3, min.features = 200, meta.data = metadata)
print("data loaded...")

# ---- Retrieve cell type identity ---- #
Idents(dataset_combined) <- dataset_combined$annotation

# ---- Keep only Healthy adult gut cells ---- #
dataset_combined <- subset(dataset_combined, subset = Diagnosis == "Healthy adult")

# ---- Save downsampled object ---- #
downsampled_gut <- subset(dataset_combined, downsample = 500)
saveRDS(downsampled_gut, "gutDownsampled.rds")


#########################################################################
# Load lung dataset (from https://www.synapse.org/#!Synapse:syn21041850)
#########################################################################

# ---- Load data ---- #
metadata <- read.csv("krasnow_hlca_10x_metadata.csv", row.names = 1, header = T)
data <- read.csv("krasnow_hlca_10x_UMIs.csv", row.names = 1, header =T)
print("data loaded...")

# ---- Create Seurat object---- #
dataset_combined <- CreateSeuratObject(counts = data, project = "lung", min.cells = 3, min.features = 200, meta.data = metadata)
dataset_combined[["percent.mt"]] <- PercentageFeatureSet(dataset_combined, pattern = "^MT-")
print("Seurat object loaded...")

# ---- Retrieve cell type identity ---- #
Idents(dataset_combined) <- dataset_combined$free_annotation

# ---- Save downsampled object ---- #
downsampled_lung <- subset(dataset_combined, downsample = 200)
saveRDS(downsampled_lung, "lungDownsampled.rds")


#########################################################################
# Load stomach dataset (from GSE183904)
#########################################################################

# ---- Create directory ---- #
dir.create("QC_single_dataset")

# ---- Load individual normal samples and create Seurat objects ---- #
df <- read.csv("GSM5573466_sample1.csv.gz", row.names = 1, header =T)
sample1 <- CreateSeuratObject(counts = df, project = "sample1", min.cells = 3, min.features = 200)
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")
print("sample1 loaded...")

df <- read.csv("GSM5573469_sample4.csv.gz", row.names = 1, header =T)
sample4 <- CreateSeuratObject(counts = df, project = "sample4", min.cells = 3, min.features = 200)
sample4[["percent.mt"]] <- PercentageFeatureSet(sample4, pattern = "^MT-")
print("sample4 loaded...")

df <- read.csv("GSM5573471_sample6.csv.gz", row.names = 1, header =T)
sample6 <- CreateSeuratObject(counts = df, project = "sample6", min.cells = 3, min.features = 200)
sample6[["percent.mt"]] <- PercentageFeatureSet(sample6, pattern = "^MT-")
print("sample6 loaded...")

df <- read.csv("GSM5573474_sample9.csv.gz", row.names = 1, header =T)
sample9 <- CreateSeuratObject(counts = df, project = "sample9", min.cells = 3, min.features = 200)
sample9[["percent.mt"]] <- PercentageFeatureSet(sample9, pattern = "^MT-")
print("sample9 loaded...")

df <- read.csv("GSM5573486_sample21.csv.gz", row.names = 1, header =T)
sample21 <- CreateSeuratObject(counts = df, project = "sample21", min.cells = 3, min.features = 200)
sample21[["percent.mt"]] <- PercentageFeatureSet(sample21, pattern = "^MT-")
print("sample21 loaded...")

df <- read.csv("GSM5573488_sample23.csv.gz", row.names = 1, header =T)
sample23 <- CreateSeuratObject(counts = df, project = "sample23", min.cells = 3, min.features = 200)
sample23[["percent.mt"]] <- PercentageFeatureSet(sample23, pattern = "^MT-")
print("sample23 loaded...")

df <- read.csv("GSM5573490_sample25.csv.gz", row.names = 1, header =T)
sample25 <- CreateSeuratObject(counts = df, project = "sample25", min.cells = 3, min.features = 200)
sample25[["percent.mt"]] <- PercentageFeatureSet(sample25, pattern = "^MT-")
print("sample25 loaded...")

df <- read.csv("GSM5573496_sample31.csv.gz", row.names = 1, header =T)
sample31 <- CreateSeuratObject(counts = df, project = "sample31", min.cells = 3, min.features = 200)
sample31[["percent.mt"]] <- PercentageFeatureSet(sample31, pattern = "^MT-")
print("sample31 loaded...")

df <- read.csv("GSM5573500_sample35.csv.gz", row.names = 1, header =T)
sample35 <- CreateSeuratObject(counts = df, project = "sample35", min.cells = 3, min.features = 200)
sample35[["percent.mt"]] <- PercentageFeatureSet(sample35, pattern = "^MT-")
print("sample35 loaded...")

# ---- QC on individual datasets ---- #
dir.create("QC_single_dataset")

pdf("QC_single_dataset/QC.pdf")
VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample1")
VlnPlot(sample4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample4")  
VlnPlot(sample6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample6")
VlnPlot(sample9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample9")
VlnPlot(sample21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample21")
VlnPlot(sample23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample23")
VlnPlot(sample25, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample25")
VlnPlot(sample31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample31")
VlnPlot(sample35, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = 'counts') +
  xlab("sample35")
dev.off()

# ---- Integrate datasets ---- #
dataset <- merge(sample1, y = c(sample4, sample6, sample9, sample21, sample23, sample25, sample31, sample35),
                 add.cell.ids = c("sample1", "sample4", "sample6", "sample9", "sample21", "sample23", "sample25", "sample31", "sample35"),
                 project = "stomach")
print("dataset merged...")
                 
dataset.list <- SplitObject(dataset, split.by = "orig.ident")
print("dataset splitted...")

dataset.list <- lapply(X = dataset.list, FUN = function(x) {
    x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 10)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
print("dataset normalized...")

features <- SelectIntegrationFeatures(dataset.list, nfeatures = 2000)

dataset.list <- lapply(X = dataset.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
print("dataset scaled...")

anchors <- FindIntegrationAnchors(object.list = dataset.list, anchor.features = features,
                                  reduction = "rpca")
print("dataset anchored...")

dataset_combined <- IntegrateData(anchorset = anchors)
print("dataset combined...")

rm(dataset.list)
rm(anchors)
rm(features)
rm(df)
rm(dataset)

DefaultAssay(dataset_combined) <- "integrated"

# ---- QC on integrated dataset ---- #
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

# ---- Visualization and Clustering ---- #
dir.create("Analysis")
dir.create("Analysis/Resolution_0,5_Dims_10")

dataset_combined <- ScaleData(dataset_combined, verbose = FALSE)
dataset_combined <- RunPCA(dataset_combined, npcs = 30, verbose = FALSE)
dataset_combined <- RunUMAP(dataset_combined, reduction = "pca", dims = 1:10)
dataset_combined <- FindNeighbors(dataset_combined, reduction = "pca", dims = 1:10)
dataset_combined <- FindClusters(dataset_combined, resolution = 0.5)

pdf("Analysis/Resolution_0,5_Dims_10/elbow_plot.pdf")
ElbowPlot(dataset_combined)
dev.off()

p1 <- DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dataset_combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("Analysis/Resolution_0,5_Dims_10/UMAPclusters.pdf")
p1 + p2
dev.off()

pdf("Analysis/Resolution_0,5_Dims_10/UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, features = c("MUC5AC"), pt.size = 1, order = T)
dev.off()

# NB: identity of Goblet cells is cluster 3

# ---- Save downsampled object ---- #
downsampled_stomach <- subset(dataset_combined, downsample = 500)
saveRDS(downsampled_stomach, "stomachDownsampled.rds")
