## ======================================================================= ##
##         Seurat analysis of human conjunctival tuft 10X scSeq          ##
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
obj <- readRDS("tuft.rds")

dir.create("Tuft/")
setwd("Tuft/")
directory <- getwd()


#########################################################################
# Clustering analysis
#########################################################################

obj <- FindVariableFeatures(obj, nfeatures = 3000, assay = "RNA")
obj <- ScaleData(obj, assay = "integrated", features = rownames(obj))
obj <- RunPCA(obj, assay = "integrated", features = obj[['RNA']]@var.features)
obj <- RunUMAP(obj, dims = 1:20)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 1)


#########################################################################
# Plot UMAPs
#########################################################################

# --- create directory to save data --- #
dir.create(paste0(directory,"/Res1/"))
dir.create(paste0(directory,"/Res1/renormalized/"))
setwd(paste0(directory,"/Res1/renormalized/"))

# --- Plot UMAPs of clusters --- #
pdf("UMAP_clusters.pdf")
DimPlot(obj, reduction = "umap", 
        cols = c("orange", "pink", "red", "yellow", "purple"),
        pt.size = 6)
DimPlot(obj, reduction = "umap", cols = c("orange", "pink", "red", "yellow", "purple"), pt.size = 6) +
  NoLegend()
dev.off()

# --- Plot UMAPs of conditions --- #
colors_conditions <- c("#C1234E", #ALI day0
                       "#B5DFD0", #ALI day17 M04
                       "#95BCAD", #ALI day17 M16
                       "#1D92A9", #ALI day17 IL4+13
                       "#902058", #ALI day3
                       "#EE8B64", #DM
                       #"#FAECDD", #EM
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
write.csv(allmarkers, "DEG_res1.tsv")

# --- UMAP of gene expression --- #
pdf("Tuftmarkers.pdf")
FeaturePlot(obj, c("NREP", "AVIL", "KRT13", "POU2F3"), combine = F, min.cutoff = 0, coord.fixed = T, pt.size = 5)
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
new.cluster.ids <- c("NREP+ tuft", #0
                     "KRT13+ tuft", #1
                     "KRT13+ tuft", #2
                     "KRT13+ tuft", #3
                     "BMX+ tuft") #4
Idents(obj) <- "seurat_clusters"
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# --- UMAP of new clusters --- #
pdf("UMAP_clusters_newID.pdf")
DimPlot(obj, reduction = "umap", cols = c("#59B045", "#BCEFB0", "#187000"), pt.size = 5)
DimPlot(obj, reduction = "umap", cols = c("#59B045", "#BCEFB0", "#187000"), pt.size = 5) +
  NoLegend()
dev.off()


#########################################################################
# Recalculate DE genes per cluster and save it
#########################################################################

allmarkers <- FindAllMarkers(obj, logfc.threshold = 0)
write.table(allmarkers, file = paste0("DEgenes_res1_allClusters_newIDs.tsv"))


#########################################################################
# Plot top20 DEGs
#########################################################################
# ---- Reorder identities ---- #              
obj@active.ident <- factor(obj@active.ident, 
                           levels=(c("BMX+ tuft", 
                                     "NREP+ tuft",
                                     "KRT13+ tuft")))

# ---- Dot plot of DEGs ---- #  
top20 <- c((allmarkers %>% filter(p_val_adj < 0.01, cluster == "KRT13+ tuft") %>% top_n(20, avg_log2FC))$gene,
           (allmarkers %>% filter(p_val_adj < 0.01, cluster == "NREP+ tuft") %>% top_n(20, avg_log2FC))$gene,
           (allmarkers %>% filter(p_val_adj < 0.01, cluster == "BMX+ tuft") %>% top_n(20, avg_log2FC))$gene)
           
pdf("Dot_DEG_per_cluster.pdf", height = 2.1, width = 16)
DotPlot(obj,
        features = unique(top20), assay = "RNA",
        dot.scale=6, scale = T, scale.min = 0) +
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +
  scale_color_gradient(low = "white", high = "black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()

# ---- Violin plot of selected markers ---- #  
pdf("Vln_selected.pdf")
VlnPlot(obj, c("KRT5", "KRT13", "BMX", "AVIL", "NREP", "TP63", "POU2F3", "SOX4", "SPDEF", "TUBB3"), 
        cols = c("#187000", "#59B045","#BCEFB0"), 
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
  
  pdf(paste0("GOterms_", celltype, ".pdf"), width = 6, height = 4)
  print(dotplot(enrich_go, font.size = 8, showCategory = 10))
  #, showCategory = selected_pathways))
  dev.off()
}

# --- Plot GO terms --- #
plot_go_subcl(allmarkers, "BMX+ tuft")
plot_go_subcl(allmarkers, "NREP+ tuft")
plot_go_subcl(allmarkers, "KRT13+ tuft")


