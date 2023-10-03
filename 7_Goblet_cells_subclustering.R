## ======================================================================= ##
##         Seurat analysis of human conjunctival goblet 10X scSeq          ##
## ======================================================================= ##

# Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(umap)
library(Matrix)
library(RColorBrewer)
library(goseq)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(clusterProfiler)


#########################################################################
# Open goblet cells object
#########################################################################
obj <- readRDS("goblet.rds")

dir.create("Goblet/")
setwd("Goblet/")
directory <- getwd()


#########################################################################
# Clustering analysis
#########################################################################
obj <- FindVariableFeatures(obj, nfeatures = 3000, assay = "RNA")
obj <- ScaleData(obj, assay = "integrated", features = rownames(obj))
obj <- RunPCA(obj, assay = "integrated", features = obj[['RNA']]@var.features)
obj <- RunUMAP(obj, dims = 1:15)
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5)


#########################################################################
# Plot UMAPs
#########################################################################

# --- create directory to save data --- #
dir.create(paste0(directory,"/Res0.5/noScaling/renormalized/"))
setwd(paste0(directory,"/Res0.5/noScaling/renormalized/"))

# --- Plot UMAPs with clusters and conditions --- #
pdf("UMAP_clusters.pdf")
DimPlot(obj, reduction = "umap", 
        cols = c("orange", "pink", "red", "yellow"),
        pt.size = 6)
DimPlot(obj, reduction = "umap", cols = c("orange", "pink", "red", "yellow"), pt.size = 6) +
  NoLegend()
dev.off()

colors_conditions <- c("#B5DFD0", #ALI day17 M04
                       "#95BCAD", #ALI day17 M16
                       "#1D92A9", #ALI day17 IL4+13
                       "#EE8B64", #DM
                       "#808181") #tissue

pdf("CellOrigin.pdf")
DimPlot(obj, group.by = "orig.ident", cols = colors_conditions, pt.size = 6)
DimPlot(obj, group.by = "orig.ident", cols = colors_conditions, pt.size = 6)+
  NoLegend()
dev.off()

#########################################################################
# Calculate DEGs and plot some
#########################################################################

# --- Calculate DEGs --- #
allmarkers <- FindAllMarkers(obj)
write.csv(allmarkers, "DEG_res0.5.tsv")

# --- UMAP of gene expression --- #
pdf("Gobletmarkers.pdf")
FeaturePlot(obj, c("MUC5AC", "TFF1", "TFF3", "SPDEF"), combine = F, min.cutoff = 0, coord.fixed = T, pt.size = 5)
dev.off()

significant_markers <- unique(filter(allmarkers, p_val_adj < 0.01)$gene)
pdf("Clustermarkers.pdf")
FeaturePlot(obj, unique(significant_markers), combine = F, min.cutoff = 0, coord.fixed = T, pt.size = 5)
dev.off()


#########################################################################
# Rename clusters
#########################################################################

# ----  Save data in a new folder ---- #
dir.create("NewIdentities/")
setwd("NewIdentities/")

# --- Rename clusters --- #
new.cluster.ids <- c("S100A8+ goblet", #0
                     "TFF3+ goblet") #1
Idents(obj) <- "seurat_clusters"
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# --- UMAP of new clusters --- #
pdf("UMAP_clusters_newID.pdf")
DimPlot(obj, reduction = "umap", cols = c("orange", "pink"), pt.size = 5)
DimPlot(obj, reduction = "umap", cols = c("orange", "pink"), pt.size = 5) +
  NoLegend()
dev.off()


#########################################################################
# Recalculate DE genes per cluster and save it
#########################################################################
      
allmarkers <- FindAllMarkers(obj, logfc.threshold = 0)
write.table(allmarkers, file = paste0("DEgenes_res0.5_allClusters_newIDs.tsv"))


#########################################################################
# Plot top20 DEGs
#########################################################################

# ---- Dot plot of DEGs ---- #  
top20 <- allmarkers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(-20, avg_log2FC) 

pdf("Dot_DEG_per_cluster.pdf", height = 2.1, width = 12)
DotPlot(obj,
        features = top20$gene, assay = "RNA",
        dot.scale=6, scale = T, scale.min = 0) +
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +
  scale_color_gradient(low = "white", high = "black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()

# ---- Violin plot of selected markers ---- #  
pdf("Vln_selected.pdf")
VlnPlot(obj, c("TFF3", "LCN2", "S100A8", "SLPI",  "MUC5AC", "SPDEF", "TFF1", "LYPD2"), 
        cols = c("orange","pink"), 
        combine = F,
        assay = "RNA", 
        pt.size = 0)
dev.off()


#########################################################################
# GO term analysis
#########################################################################

# --- Function to plot GO terms --- #
plot_go_subcl <- function(DEGdataframe, celltype) {
  
  genes_to_test <<- filter(DEGdataframe, p_val_adj < 0.01, avg_log2FC > 0.25, cluster == celltype)$gene
  
  enrich_go <<- enrichGO(gene = genes_to_test,
                         OrgDb = "org.Hs.eg.db",
                         keyType = "SYMBOL",
                         ont = "BP")
  selected_pathways <- enrich_go[c(1, 2, 5, 6, 9, 10, 11, 13, 14, 15),]$Description
  
  pdf(paste0("GOterms_", celltype, ".pdf"), width = 6, height = 4)
  print(dotplot(enrich_go, font.size = 8, showCategory = selected_pathways))
  #, showCategory = selected_pathways))
  dev.off()
}

# --- Plot GO terms --- #
plot_go_subcl(allmarkers, "S100A8+ goblet")
plot_go_subcl(allmarkers, "TFF3+ goblet")


