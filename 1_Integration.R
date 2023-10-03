## ======================================================================= ##
##       integration of human conjunctiva single-cell RNAseq datasets      ##
## ======================================================================= ##

#########################################################################
# Load libraries
#########################################################################

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
library(remotes,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/")
library(Seurat,lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/") 

directory <- "/hpc/hub_clevers/MB/conjunctiva/Conjunctiva_10X/"
setwd(directory)

dir.create("Final_analysis/")
setwd("Final_analysis/")

directory <- getwd()


#########################################################################
# Load all datasets
#########################################################################

dir.create("QC_single_dataset/")
setwd("QC_single_dataset/")

# ---- Load tissue dataset ---- #
LX237_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX237/an_392/outs/raw_feature_bc_matrix/")
colnames(LX237_before) <- paste0(colnames(LX237_before), "-1")

LX237 <- CreateSeuratObject(counts = LX237_before, 
                            project = "tissue", 
                            min.cells = 3, 
                            min.features = 200)
rm(LX237_before)
LX237[["percent.mt"]] <- PercentageFeatureSet(LX237, pattern = "^MT-")

# ---- Add metadata to tissue: genotype and doublet information ---- #
# Load file with genotype information (obtained from Souporcell)
# https://github.com/wheaton5/souporcell
genotype_df <- read.csv("/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX237/an_393/k4/clusters.tsv", 
                         sep = "\t", 
                         row.names = 1)
print("genotype loaded")

# Subset the dataset to only contain the data for which we have genotype information
LX237$CellName <- colnames(LX237)
LX237 <- subset(LX237, 
                subset = CellName %in% rownames(genotype_df))
print("LX237 subsetted")

# Subset the genotype dataframe to only contain data for cells in the filtered object
genotype_df_sub <- genotype_df[rownames(genotype_df) %in% LX237$CellName,]
print("genotype subsetted")

# Add metadata column to object
LX237$genotype <- genotype_df_sub$assignment # genotype ("0", "1", "2", "3", or mixed)
LX237$doubletstatus <- genotype_df_sub$status # doublet, singlet or unassigned
print("metadata added")

# Plot QC
pdf("QC_tissue_beforedoubletfiltering.pdf")
VlnPlot(LX237, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        slot = 'counts', 
        group.by = "doubletstatus",
        pt.size = 0.1)
VlnPlot(LX237, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        slot = 'counts',
        pt.size = 0.1)
dev.off()

# Keep only the singlets
LX237 <- subset(LX237, 
                subset = doubletstatus == "singlet")
print("Tissue dataset: doublets removed...")

# ---- Load ALI day17 M04 dataset ---- #
LX243_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX243/an_394/outs/filtered_feature_bc_matrix/")
LX243 <- CreateSeuratObject(counts = LX243_before, 
                            project = "ALId17_hCjM04", 
                            min.cells = 3, 
                            min.features = 200)
rm(LX243_before)
LX243[["percent.mt"]] <- PercentageFeatureSet(LX243, pattern = "^MT-")


# ---- Load ALI day17 M16 dataset ---- #
LX244_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX244/an_395/outs/filtered_feature_bc_matrix/")
LX244 <- CreateSeuratObject(counts = LX244_before, 
                            project = "ALId17_hCjM16", 
                            min.cells = 3, 
                            min.features = 200)
rm(LX244_before)
LX244[["percent.mt"]] <- PercentageFeatureSet(LX244, pattern = "^MT-")


# ---- Load Organoid EM dataset ---- #
TX294_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX294/count/sample_filtered_feature_bc_matrix/")
TX294 <- CreateSeuratObject(counts = TX294_before$`Gene Expression`, 
                            project = "Organoids_EM", 
                            min.cells = 3, 
                            min.features = 200)
rm(TX294_before)
TX294[["percent.mt"]] <- PercentageFeatureSet(TX294, pattern = "^MT-")


# ---- Load Organoid DM dataset ---- #
TX295_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX295/count/sample_filtered_feature_bc_matrix/")
TX295 <- CreateSeuratObject(counts = TX295_before$`Gene Expression`, 
                            project = "Organoids_DM",
                            min.cells = 3, 
                            min.features = 200)
rm(TX295_before)
TX295[["percent.mt"]] <- PercentageFeatureSet(TX295, pattern = "^MT-")


# ---- Load ALI day 0 dataset ---- #
TX296_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX296/count/sample_filtered_feature_bc_matrix/")
TX296 <- CreateSeuratObject(counts = TX296_before$`Gene Expression`, 
                            project = "ALId0", 
                            min.cells = 3, 
                            min.features = 200)
rm(TX296_before)
TX296[["percent.mt"]] <- PercentageFeatureSet(TX296, pattern = "^MT-")


# ---- Load ALI day 3 dataset ---- #
TX297_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX297/count/sample_filtered_feature_bc_matrix/")
TX297 <- CreateSeuratObject(counts = TX297_before$`Gene Expression`, 
                            project = "ALId3", 
                            min.cells = 3, 
                            min.features = 200)
rm(TX297_before)
TX297[["percent.mt"]] <- PercentageFeatureSet(TX297, pattern = "^MT-")


# ---- Load ALI day 17 with IL4/13 dataset ---- #
TX298_before <- Read10X(data.dir = "/hpc/hub_clevers/sequencing_data/MB_single_cell/analysis/LX255_LX256/an_403/outs/per_sample_outs/TX298/count/sample_filtered_feature_bc_matrix/")
TX298 <- CreateSeuratObject(counts = TX298_before$`Gene Expression`, 
                            project = "ALId17_IL4+13", 
                            min.cells = 3, 
                            min.features = 200)
rm(TX298_before)
TX298[["percent.mt"]] <- PercentageFeatureSet(TX298, pattern = "^MT-")


print("datasets loaded...")



#########################################################################
# QC on individual datasets
#########################################################################

dir.create(paste0(directory, "/QC_single_dataset/"))
setwd(paste0(directory, "/QC_single_dataset/"))

pdf("QC_list.pdf")
lapply(list(LX237, LX243, LX244, TX294, TX295, TX296, TX297, TX298), 
       VlnPlot, 
       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
       ncol = 3, 
       slot = "counts",
       pt.size = 0.1)
dev.off()



#########################################################################
# Integrate datasets
#########################################################################
   
# ---- Merge datasets ---- #
dataset <- merge(LX237, y = c(LX243, LX244, TX294, TX295, TX296, TX297, TX298),
                 add.cell.ids = c("Tissue", "ALId17_hCjM04", "ALId17_hCjM16",
                                "Organoids_EM", "Organoids_DM", "ALId0", "ALId3",
                                "ALId17_IL4+13"),
                 project = "conjunctiva")

# ---- Create dataset list ---- #                 
dataset.list <- SplitObject(dataset, split.by = "orig.ident")
print("datasets merged...")

# ---- Normalize and identify variable features for each dataset independently---- #
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
    x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 10)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
print("datasets normalized...")

# ---- Select variable features for integration ---- # 
features <- SelectIntegrationFeatures(dataset.list, nfeatures = 3000)

dataset.list <- lapply(X = dataset.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
print("datasets scaled...")

# ---- Find integration anchors ---- # 
anchors <- FindIntegrationAnchors(object.list = dataset.list, 
                                  anchor.features = features,
                                  reduction = "rpca")
print("dataset anchored...")

# ---- Integrate datasets ---- #
dataset_combined <- IntegrateData(anchorset = anchors)
print("dataset integrated")

# ---- Remove heavy objects ---- #
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

# ---- Analyse the integrated data ---- # 
DefaultAssay(dataset_combined) <- "integrated"



#########################################################################
# Perform QC on integrated dataset
#########################################################################

dir.create("QC_merged_dataset")

# ---- Violin plot - QC ---- # 
pdf("QC_merged_dataset/QC.pdf", width = 14)
VlnPlot(dataset_combined,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        slot = 'counts',
        pt.size = 0.1)
VlnPlot(dataset_combined, 
        group.by = "orig.ident", 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0.1, 
        ncol = 3) +
    NoLegend()
dev.off()


# ----  Scatter plot - QC ---- # 
pdf("QC_merged_dataset/Scatter-counts.pdf", width = 14)
FeatureScatter(dataset_combined, 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mt",
               pt.size = 0.1)
FeatureScatter(dataset_combined, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               pt.size = 0.1)
dev.off()



#########################################################################
# Visualization and Clustering - dimensions = 40 
#########################################################################

# ----  Prepare working directory ---- # 
dir.create(paste0(directory, "/Analysis/"))
dir.create(paste0(directory,"/Analysis/Dim40/"))
setwd(paste0(directory,"/Analysis/Dim40/"))

# ----  Scale, UMAP ---- # 
dataset_combined <- ScaleData(dataset_combined, verbose = FALSE)
dataset_combined <- RunPCA(dataset_combined, verbose = FALSE)
print("Dataset scaled...")

# ----  Elbow plot ---- #  
pdf("elbow_plot.pdf")
ElbowPlot(dataset_combined, ndims = 50)
dev.off()

# ----  UMAP, neighbors ---- # 
dataset_combined <- RunUMAP(dataset_combined, reduction = "pca", dims = 1:40)
dataset_combined <- FindNeighbors(dataset_combined, reduction = "pca", dims = 1:40)
print("Neighbors found...")

# ----  Clustering with variable resolutions ---- # 
dataset_combined <- FindClusters(dataset_combined, resolution = c(0.7, 0.9, 1.2, 1.5, 2, 3))
print("Dataset clustered")


#########################################################################
# Save merged environment
#########################################################################
save.image("envir_with_Resolutions.RData")


#########################################################################
# Plot UMAP with each resolution showing clusters and MUC5AC expression
#########################################################################

Idents(dataset_combined) <- dataset_combined$integrated_snn_res.0.7
p1 <- DimPlot(dataset_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE) +
        ggtitle('Res_0.7')
        
Idents(dataset_combined) <- dataset_combined$integrated_snn_res.0.9
p2 <- DimPlot(dataset_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE) +
        ggtitle('Res_0.9')

Idents(dataset_combined) <- dataset_combined$integrated_snn_res.1.2
p3 <- DimPlot(dataset_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE) +
        ggtitle('Res_1.2')

Idents(dataset_combined) <- dataset_combined$integrated_snn_res.1.5
p4 <- DimPlot(dataset_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE) +
        ggtitle('Res_1.5')

Idents(dataset_combined) <- dataset_combined$integrated_snn_res.2
p5 <- DimPlot(dataset_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE) +
        ggtitle('Res_2')
                  
Idents(dataset_combined) <- dataset_combined$integrated_snn_res.3
p6 <- DimPlot(dataset_combined, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE) +
        ggtitle('Res_3')       
                  
pdf("Clusters_Res.pdf", width = 15, height = 15)
p1 + p2 + p3 + p4 + p5 + p6 + NoLegend()
dev.off()

pdf("UMAP_MUC5AC.pdf")
FeaturePlot(dataset_combined, 
            features = c("MUC5AC"), 
            pt.size = 0.1, 
            order = T)
dev.off()

