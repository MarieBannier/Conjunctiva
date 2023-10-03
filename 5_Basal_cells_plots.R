## ======================================================================= ##
##         Seurat analysis of human conjunctival goblet 10X scSeq          ##
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

directory <- "/hpc/hub_clevers/MB/conjunctiva/Conjunctiva_10X/Final_analysis/NewIdentities/BasalCells/Res0.7/"
setwd(directory)


#########################################################################
# Load workspace
#########################################################################

load("envir_basalsubcl.RData")


#########################################################################
# Color palettes
#########################################################################
# ----  Colors for conditions ---- #
colors8 <- c("#faa7d6", #ALI day0
             "#87235c", #ALI day17 M04
             "#87235c", #ALI day17 M16
             "#f5a511", #ALI day17 IL4+13
             "#c95195", #ALI day3
             "#add6ff", #DM
             "#a9ded1", #EM
             "#e3ca2b") #tissue

# ----  Colors for new clusters ---- #
colors20 <- c("#EB9700", "#D93B8D", "#7C60A4", "#218A9E", "#0E5245",
              "#FFB8B8", "#A02E21", "#E02F28", "#556BAC", "#51D8E8",
              "#FF4F30", "#E856C9", "#0994E6", "#E8460C", "#FF0653",
              "#FF8A5E", "#E478FF", "#FFD536", "#E86556", "#7933E8")


#########################################################################
# Rename clusters
#########################################################################

# ----  Save data in a new folder ---- #
dir.create("NewIdentities/")
setwd("NewIdentities/")

# ----  Rename clusters ---- # 
new.cluster.ids <- c("SERPINF1+ basal", #0
                      "SERPINF1+ basal", #1
                      "SERPINF1+ basal", #2
                      "SERPINF1+ basal", #3
                      "KRT17+ basal", #4
                      "SERPINF1+ basal", #5
                      "DPP6+ basal", #6
                      "HELLS+ cycling", #7
                      "KRT23+ early keratinocytes", #8
                      "PLAU+ basal", #9
                      "PTTG1+ cycling", #10
                      "PCLAF+ cycling", #11
                      "CXCL17+ early keratinocytes", #12
                      "PTPRJ+ early keratinocytes", #13
                      "SERPINF1+ basal", #14
                      "LGR4+ basal", #15
                      "SERPINB1+ early keratinocytes") #16

markers.per.cl <- c("TP63", "COL17A1", "KRT5", "KRT14",
                    "PTTG1", "MKI67", "PCLAF", "RRM2", "PCNA", "HELLS", "DEFB1", "FABP5",
                    "SERPINF1", "ROBO2", "ROBO1", "KRT17",
                    "DPP6", "LGR4", "LGR5", "LGR6", "EGFR", "NGFR",
                    "PLAU", "PLAUR", "CCL20", "CXCL1", "CXCL8", "AREG", "FGFBP1", "PLAT", "CXCL2",
                    "LYPD2", "LCN2", "SLPI",
                    "F3", "PTPRG", "PTPRJ", "KRT4", "EGF", "GCNT3", "FGFR2", "RUNX1",
                    "SERPINB1",
                    "KRT23","S100A4", "PSCA", "UPK1B")
                    
Idents(obj) <- "seurat_clusters"

names(new.cluster.ids) <- levels(obj)

obj <- RenameIdents(obj, new.cluster.ids)
message("Name of clusters:")
print(unique(Idents(obj)))

# ---- Reorder identities ---- #              
obj@active.ident <- factor(obj@active.ident, 
                           levels=rev(c("PTTG1+ cycling", 
                                        "PCLAF+ cycling",
                                        "HELLS+ cycling", 
                                        "LGR4+ basal", 
                                        "SERPINF1+ basal", 
                                        "DPP6+ basal", 
                                        "PLAU+ basal", 
                                        "KRT17+ basal", 
                                        "SERPINB1+ early keratinocytes", 
                                        "PTPRJ+ early keratinocytes",
                                        "CXCL17+ early keratinocytes",
                                        "KRT23+ early keratinocytes")))
                                                     

# ---- Plot UMAPs ---- #      
pdf("UMAP_clusters.pdf")
DimPlot(obj, reduction = "umap", cols = colors20, pt.size = 0.5, label = TRUE)
DimPlot(obj, reduction = "umap", cols = colors20, pt.size = 0.5) +
  NoLegend()
dev.off()

pdf("CellOrigin.pdf")
DimPlot(obj, group.by = "orig.ident", cols = colors8, pt.size = 0.5) #tissue
DimPlot(obj, group.by = "orig.ident", cols = colors8,
        pt.size = 0.5) + 
  NoLegend()
dev.off()


#########################################################################
# Recalculate DE genes per cluster and save it (or re-open it
#########################################################################

# ---- Calculate and save DEGs ---- # 
allmarkers <- FindAllMarkers(obj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.3)
write.csv(allmarkers, file = paste0("DEgenes_res0.7_newIDs_allClusters.csv"))

# ---- Re-load DEGs ---- # 
#allmarkers <- read.csv("DEgenes_res0.7_newIDs_allClusters.csv")


#########################################################################
# Plot selected DEGs
#########################################################################

# ---- List of markers ---- # 
markers <- rev(c("PTTG1", "MKI67", "CDC20", #PTTG1+ cycling
                 "PCLAF", "RRM2", "HIST1H4C", "STMN1", #PCLAF+ cycling
                 "HELLS", "KRT6A", "FABP5", #HELLS+ cycling
                 "LGR4", "IGF1R", "FOXP2", "FOXP1", "ROBO2", "ROBO1", "COL17A1", #LGR4+ basal
                 "SERPINF1", "KRT15", "IGFBP2", #SERPINF1+ basal
                 "DPP6", "KRT31", #DPP6+ basal
                 "PLAU", "PLAUR", "PLAT", "CXCL1", "CXCL2", "CXCL8", "CXCL10", "CCL20", "LAMC2", "AREG", "FGFBP1",#PLAU+ basal
                 "KRT17",  "TNFRSF12A", #KRT17+ basal
                 "SERPINB1", "S100A9", "S100A8", "LYPD2", "AQP5", #SERPINB1+ early keratinocytes
                 "PTPRJ", "PTPRK", "PTPRG", "STK39", "ADARB2", "MUC16", "ATP10B", "IKZF2",  #PTPRJ+ early keratinocytes
                 "S100A4", "HES4", "CXCL17", "F3", "LCN2",  #CXCL17+ early keratinocytes
                 "KRT23", "SLPI", "UPK1B")) #KRT23+ early keratinocytes

# ---- Plot markers per cluster ---- # 
pdf("Dot_markers_per_cluster.pdf", height = 4, width = 18.5)
DotPlot(obj, 
        features = markers, 
        assay="RNA",
        dot.scale=6) +
  RotatedAxis() +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +
        scale_color_gradient(low = "white", high = "black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
DotPlot(obj, 
        features = markers,
        assay="RNA", 
        dot.scale=6) +
  RotatedAxis() + NoLegend() +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +
        scale_color_gradient(low = "white", high = "black") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()
