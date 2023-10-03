## ======================================================================= ##
##         Seurat analysis of human conjunctival basal cells 10X scSeq     ##
## ======================================================================= ##

#########################################################################
# Load libraries
#########################################################################

library(viridisLite) 
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
library(viridis) 
library(pheatmap) 
library(rlang) 


#########################################################################
# Open basal cells object
#########################################################################
obj <- readRDS("basal.rds")


#########################################################################
# Recluster with higher resolution (without scaling)
#########################################################################

# ---- Create directory to save basal cell data ---- #
dir.create("BasalCells/")
setwd("BasalCells/")

# ---- Clustering ---- #
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:15)
obj <- FindNeighbors(obj, dims = 1:15)
DefaultAssay(obj) <- "integrated"
obj <- FindClusters(obj, resolution = 0.7)

# ---- Create directory ---- #
dir.create("Res0.7/")
setwd("Res0.7/")


#########################################################################
# Data analysis
#########################################################################
# ----  Colors for new clusters ---- #
colors20 <- c("#EB9700", "#D93B8D", "#7C60A4", "#218A9E", "#0E5245",
              "#FFB8B8", "#A02E21", "#E02F28", "#556BAC", "#51D8E8",
              "#FF4F30", "#E856C9", "#0994E6", "#E8460C", "#FF0653",
              "#FF8A5E", "#E478FF", "#FFD536", "#E86556", "#7933E8")

# ----  Plot UMAPs of clusters ---- #
pdf("UMAP_clusters.pdf")
DimPlot(obj, reduction = "umap", cols = colors20, pt.size = 0.5)
DimPlot(obj, reduction = "umap", cols = colors20, pt.size = 0.5) +
  NoLegend()
dev.off()

# ----  Plot UMAPs of conditions ---- #
pdf("CellOrigin.pdf")
DimPlot(obj, group.by = "orig.ident", cols = c("#FFF58F", #ALI day0
                                               "#BE97E8", #ALI day17 M04
                                               "#BE97E8", #ALI day17 M16
                                               "#0DA2FF", #ALI day17 IL4+13
                                               "#EBC300", #ALI day3
                                               "#EB8200", #DM
                                               "#a9ded1", #EM
                                               "#CDFFD2"), pt.size = 0.5) #tissue
DimPlot(obj, group.by = "orig.ident", cols = c("#FFF58F", #ALI day0
                                               "#BE97E8", #ALI day17 M04
                                               "#BE97E8", #ALI day17 M16
                                               "#0DA2FF", #ALI day17 IL4+13
                                               "#EBC300", #ALI day3
                                               "#EB8200", #DM
                                               "#a9ded1", #EM
                                               "#CDFFD2"),#tissue
        pt.size = 0.5) + 
  NoLegend()
dev.off()

# ----  Plot UMAPs of basal cell markers ---- #
pdf("BasalCellMarkersUMAP.pdf", width = 10, height = 15)
FeaturePlot(obj, features = c("TP63", "COL17A1", "MKI67", "ROBO1", "ROBO2", 
                              "IGFBP2", "KRT5", "KRT14", "KRT23", "S100A9",
                              "EGFR", "LGR5", "LGR4", "LGR6", "TNFRSF19", 
                              "NGFR", "KRT6A"), 
            order = T, pt.size = 0.5, combine = FALSE, coord.fixed = TRUE, min.cutoff = 0)
dev.off()


#########################################################################
# get DE genes per cluster and save it
#########################################################################

allmarkers <- FindAllMarkers(obj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0)
write.csv(allmarkers, file = paste0("DEgenes_res0.7_allClusters.csv"))


#########################################################################
# UMAP of all DEG per cluster
#########################################################################

# ---- Identify significant DEGs ---- #
df <- allmarkers[allmarkers$p_val_adj < 0.01 & 
                 allmarkers$avg_logFC > 0.5,]

genes_to_plot <- df$gene

# ---- Plot significant DEGs ---- #
for (i in genes_to_plot) {
    pdf(paste0("UMAP_", i,".pdf"))
    FeaturePlot(obj, features = i, 
                order = T, pt.size = 1, min.cutoff = 0)
    dev.off()
    }  


#########################################################################
# Cell proportions per cell type in clusters
#########################################################################

write.csv(table(Idents(obj)), "cell_nb.csv")
write.csv(prop.table(table(Idents(obj))), "cell_prop.csv")


#########################################################################
# Save workspace
#########################################################################

save.image("envir_basalsubcl.RData")
