## ======================================================================= ##
##                  integration of human conjunctiva datasets              ##
## ======================================================================= ##

#########################################################################
# Load libraries
#########################################################################

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
library(remotes)
library(Seurat) 
library(tidyr) 


#########################################################################
# Load and split STOMACH dataset (Goblet cluster = "3")
#########################################################################

# ---- Read in data ---- #
stomach <- readRDS("stomachDownsampled.rds")

# ---- Add tissue type as orig.ident (to facilitate splitting) ---- #
stomach$condition <- stomach$orig.ident
stomach$orig.ident <- "stomach"

# ---- Add cell type as additional meta.data ---- #
stomach$CellType <- stomach$seurat_clusters

# ---- Revert to RNA assay & Remove previous integration data---- #
DefaultAssay(stomach) <- "RNA"
stomach[["integrated"]] <- NULL

print(dim(stomach))

message("Stomach loaded...")


#########################################################################
# Load and split LUNG dataset (Goblet cluster called "Goblet")
#########################################################################

# ---- Read in data ---- #
lung <- readRDS("lungDownsampled.rds")

# ---- Add tissue type as orig.ident (to facilitate splitting and integration) ---- #
lung$condition <- lung$orig.ident
lung$orig.ident <- "lung"

# ---- Add cell type as additional meta.data ---- #
lung$CellType <- lung$free_annotation

# ---- Set to RNA assay ---- #
DefaultAssay(lung) <- "RNA"
print(dim(lung))

message("Lung loaded...")


#########################################################################
# Load and split GUT dataset (Goblet clusters called "Goblet cell" and "BEST2+ Goblet cell")
#########################################################################

# ---- Read in data ---- #
gut <- readRDS("gutDownsampled.rds")

# ---- Remove mt genes ---- #
gut <- subset(gut, features = grep("^MT-", rownames(gut), invert = TRUE, value = TRUE))

# ---- Keep cells found in the small & large intestine (i.e. exclude lymph nodes etc.) ---- #
gut <- subset(gut, subset = Region %in% c("LargeInt", "SmallInt"))

# ---- Add tissue type as orig.ident (to facilitate splitting) ---- #
gut$condition <- gut$Sample.name
gut$orig.ident <- gut$Region

# ---- Add cell type as additional meta.data ---- #
gut$CellType <- gut$annotation

# ---- Set to RNA assay ---- #
DefaultAssay(gut) <- "RNA"
print(dim(gut))

message("Gut loaded...")


#########################################################################
# Load and split CONJUNCTIVA dataset (Goblet cluster called "goblet cells")
#########################################################################

# ---- Read in data ---- #
conjunctiva <- readRDS("conjunctivaDownsampled.rds")

# ---- Set to RNA assay ---- #
DefaultAssay(conjunctiva) <- "RNA"
conjunctiva[["integrated"]] <- NULL

# ---- Modify conjunctiva conditions to be able to integrate few cells ---- #
conjunctiva$CellType <- conjunctiva@active.ident

conjunctiva$condition <- conjunctiva$orig.ident
conjunctiva$orig.ident <- "conjunctiva"

print(dim(conjunctiva))

message("Conjunctiva loaded...")


#########################################################################
# QC on individual datasets
#########################################################################

dir.create("Merged_Goblets/")
setwd("Merged_Goblets/")

directory <- getwd()

dir.create("QC_single_dataset/")
setwd("QC_single_dataset/")

pdf("QC_list.pdf")
lapply(list(stomach, lung, gut, conjunctiva), 
       VlnPlot, 
       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
       ncol = 3, 
       slot = "counts",
       pt.size = 0.3)
dev.off()


#########################################################################
# Merge datasets
#########################################################################

# ---- Merge datasets ---- #
dataset_combined <- merge(stomach, y = c(lung, gut, conjunctiva),
                          add.cell.ids = c("Stomach",
                                           "Lung", 
                                           "Gut",
                                           "Conjunctiva"),
                          project = "Goblets")

message("Datasets merged...")

saveRDS(dataset_combined, "mergedDatasets.rds")

# ---- Normalize and identify variable features for each dataset independently---- #
dataset_combined <- subset(dataset_combined, subset = nFeature_RNA > 500 & nFeature_RNA < 7500)
dataset_combined <- NormalizeData(dataset_combined)
dataset_combined <- FindVariableFeatures(dataset_combined, selection.method = "vst", nfeatures = 3000)
dataset_combined <- ScaleData(dataset_combined, verbose = FALSE)
dataset_combined <- RunPCA(dataset_combined, verbose = FALSE, npcs = 50)

message("datasets scaled & normalized...")

# ---- Remove heavy objects ---- #
rm(dataset.list)
rm(anchors)
rm(features)
rm(stomach)
rm(lung)
rm(gut)
rm(conjunctiva)
rm(dataset)


#########################################################################
# Perform QC on integrated dataset
#########################################################################

dir.create(paste0(directory, "/QC_merged_dataset"))
setwd(paste0(directory, "/QC_merged_dataset"))

# ---- Violin plot - QC ---- # 
pdf("QC.pdf", width = 14)
VlnPlot(dataset_combined,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        slot = 'counts',
        pt.size = 1)
VlnPlot(dataset_combined, 
        group.by = "orig.ident", 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 1, 
        ncol = 3) +
  NoLegend()
dev.off()


# ----  Scatter plot - QC ---- # 
pdf("Scatter-counts.pdf", width = 14)
FeatureScatter(dataset_combined, 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mt",
               pt.size = 1)
FeatureScatter(dataset_combined, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               pt.size = 1)
dev.off()


#########################################################################
# Visualization and Clustering - dimensions = 20
#########################################################################

# ----  Prepare working directory ---- # 
dir.create(paste0(directory, "/Analysis/"))
setwd(paste0(directory,"/Analysis/"))

# ----  UMAP, neighbors ---- # 
dataset_combined <- RunUMAP(dataset_combined, reduction = "pca", dims = 1:20)
dataset_combined <- FindNeighbors(dataset_combined, reduction = "pca", dims = 1:20)
print("Neighbors found...")


#########################################################################
# Plot UMAP of tissue of origin in whole dataset
#########################################################################
# ---- Reorder identities ---- #              
dataset_combined$orig.ident <- factor(dataset_combined$orig.ident, 
                                      levels=rev(c("conjunctiva", 
                                                   "SmallInt",
                                                   "LargeInt", 
                                                   "lung", 
                                                   "stomach")))

# ---- UMAP of tissue of origin ---- #  
pdf("UMAP_tissueOforigin.pdf")
print(DimPlot(dataset_combined, 
              reduction = "umap",
              pt.size = 2,
              group.by = "orig.ident",
              cols = c("coral2", "cornflowerblue", "darkred", "darkolivegreen3", "aquamarine")))
dev.off()


#########################################################################
# Subset Goblet cells only
#########################################################################

goblet_cells <- subset(dataset_combined, subset = CellType %in% c("Goblet", "3", "Goblet cell", "BEST2+ Goblet cell", "goblet cells"))

# ---- UMAP of tissue of origin ---- #  
pdf("UMAP_tissueOforigin_goblet.pdf")
print(DimPlot(goblet_cells, 
              reduction = "umap",
              pt.size = 2,
              group.by = "orig.ident",
              cols = c("coral2", "cornflowerblue", "darkred", "darkolivegreen3", "aquamarine")))
dev.off()

markers <- c("TFF1", "TFF3", "SPDEF", "MUC5AC", "BEST2", "WFDC2", "LCN2", 
             "MUC12", "F3", "C3", "S100A8", "S100A9", "S100A11", 
             "LYPD2", "LYPD", "PSCA", "SLPI", "FCGBP", "ZG16", "AQP5", "AQP3", 
             "SCGB3A1", "SFTPC", "BPIFB1", "LIPF", "TFF2", "PGC", "GKN1", "GKN2", 
             "MUC2", "MUC5B", "MUC1", "MKI67", "SCGB1A1",
             "SELENBP1", "SET", "CCL24", "LYZ",
             "ITLN1", "SOX4", "MSMB", "S100P", 
             "CCSER1", "AZGP1", "HES6", "POU2F3",
             "NREP", "FOXI1", "STMN1", "RASSF6")


pdf("UMAP_Gobletmarkers_goblet.pdf")
print(FeaturePlot(goblet_cells, features = markers, 
                  order = T, pt.size = 1, combine = F, min.cutoff = 0))
dev.off()


pdf("Dot_Gobletmarkers.pdf", height = 3.5, width = 15)
print(DotPlot(goblet_cells, features = markers, 
              dot.scale=6, group.by = "orig.ident", assay = "RNA")+ 
        scale_color_gradient(low = "grey98", high = "black") +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
        RotatedAxis())
dev.off()

pdf("Dot_Gobletmarkers_limited.pdf", height = 3, width = 6)
print(DotPlot(goblet_cells, features = c("TFF1", "TFF2", "TFF3", "SPDEF", "SOX4", 
                                         "MKI67", "MUC1", "MUC2", "MUC5AC", "MUC5B"), 
              dot.scale=6, group.by = "orig.ident", assay = "RNA")+ 
        scale_color_gradient(low = "grey98", high = "black") +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
        RotatedAxis())
dev.off()

# ---- Find tissue-specific markers ---- #
Idents(goblet_cells) <- "orig.ident"
DEG <- FindAllMarkers(goblet_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(DEG, "DEG.csv")

# ---- Identify top markers per cluster and remove mitochondrial genes ---- # 
DEG <- DEG[grep("^MT-", rownames(DEG), invert = T),]
print(head(DEG))

DEG %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
print("top10 markers identified")


# ---- Plot top markers per cluster ---- # 
pdf("Dot_top10_DEG_per_cluster_hi.pdf", width = 13, height = 2.6)
DotPlot(goblet_cells, 
        features = unique(top10$gene), 
        dot.scale=6, 
        group.by = "orig.ident")+ 
  scale_color_gradient(low = "grey98", high = "black") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
  RotatedAxis()
DotPlot(goblet_cells, 
        features = unique(top10$gene), 
        dot.scale=6, 
        group.by = "orig.ident") + 
  scale_color_gradient(low = "grey98", high = "black") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
  RotatedAxis() + NoLegend()
dev.off()

print("dotplot of top10 markers: done")


# ---- Plot goblet subclusters markers ---- # 
goblet_subcl_genes <- read.csv("Goblet_subcl_genes.csv", header =T)

pdf("Dot_Goblet_subcl.pdf", width = 10, height = 3)
DotPlot(dataset_combined, 
        features = goblet_subcl_genes, 
        dot.scale=6, 
        group.by = "orig.ident")+ 
  scale_color_gradient(low = "grey98", high = "black") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
  RotatedAxis()
DotPlot(dataset_combined, 
        features = goblet_subcl_genes, 
        dot.scale=6, 
        group.by = "orig.ident") + 
  scale_color_gradient(low = "grey98", high = "black") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
  RotatedAxis() + NoLegend()
dev.off()
print("dotplot of subcluster markers: done")

# ---- Obtain correlation between clusters ---- # 
av.exp <- AverageExpression(goblet_cells)$RNA
head(av.exp)
message("Avg Exp calculated...")

cor.exp <- as.data.frame(cor(av.exp))
head(cor.exp)
cor.exp$x <- rownames(cor.exp)
head(cor.exp)

cor.df <- tidyr::gather(data = cor.exp, y, correlation, c('stomach', 'lung', 'LargeInt', 'SmallInt', 'conjunctiva'))
head(cor.df)

cor.df$x <- as.character(cor.df$x)
cor.df$y <- as.character(cor.df$y)
cor.df$x <- factor(cor.df$x, levels = c('lung', 'stomach', 'conjunctiva', 'LargeInt', 'SmallInt'))
cor.df$y <- factor(cor.df$y, levels = c( 'lung', 'stomach','conjunctiva', 'LargeInt', 'SmallInt'))

pdf("Heatmap_cluster_cor.pdf")
print(ggplot(cor.df, aes(x, y, fill = correlation)) +
        geom_tile() +
        scale_fill_gradientn(colors = c("white", "#d53e4f"))+
        geom_raster()+
        #scale_fill_continuous(expand = c(0,0)) +
        coord_fixed()+
        theme_classic()
)
dev.off()



