## ======================================================================= ##
##              First clustering of human conjunctiva datasets             ##
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

directory <- "/hpc/hub_clevers/MB/conjunctiva/Conjunctiva_10X/Final_analysis/Analysis/Dim40/"
setwd(directory)


#########################################################################
# Load  environment
#########################################################################
load("envir_with_Resolutions.RData")
print("Environment loaded")


#########################################################################
# Plot UMAP with each resolution showing clusters and MUC5AC expression
#########################################################################

# ----  Set resolution at 2 ---- # 
Idents(dataset_combined) <- dataset_combined$integrated_snn_res.2

# ----  Calculate DEGs ---- #
dataset_combined.markers <- FindAllMarkers(dataset_combined, 
                                           only.pos = TRUE, 
                                           min.pct = 0.25, 
                                           logfc.threshold = 0.5)
write.csv(dataset_combined.markers, "markers_AllClusters.csv")

# ----  50 colors palette ---- #
colors50 <- c("#4287f5", "#F5D033", "#84C3BE", "#4287f5", "#CF3476",
              "#FFA420", "#A98307", "#FFFF00", "#31372B", "#5D9B9B",
              "#E6D690", "#1D1E33", "#4D5645", "#3E3B32", "#1D334A",
              "#C2B078", "#F54021", "#D84B20", "#922B3E", "#025669",
              "#82898F", "#C1876B", "#C51D34", "#3F888F", "#FF7514",
              "#412227", "#9D9101", "#A65E2E", "#00BB2D", "#D36E70",
              "#57A639", "#E7EBDA", "#1B5583", "#F5D033", "#7D8471",
              "#E55137", "#DC9D00", "#256D7B", "#CDA434", "#015D52",
              "#E6D690", "#78858B", "#606E8C", "#008F39", "#20603D",
              "#F39F18", "#763C28", "#E5BE01", "#2D572C", "#49678D")

# ----  Plot primary UMAP clusters ---- #
pdf("UMAP_clusters_res2.pdf") 
DimPlot(dataset_combined, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE,
        cols = colors50) +
  ggtitle('Res_2')
dev.off()

# ---- Identify top markers per cluster ---- # 
dataset_combined.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_logFC) -> top10
print("top10 markers identified")

# ---- Plot top markers per cluster ---- # 
pdf("Dot_top10_DEG_per_cluster.pdf", width = 60, height = 20)
DotPlot(dataset_combined, 
        features = unique(top10$gene), 
        dot.scale=6, 
        group.by = "integrated_snn_res.2") +
  RotatedAxis()
DotPlot(dataset_combined, 
        features = unique(top10$gene), 
        dot.scale=6, 
        group.by = "integrated_snn_res.2") +
  RotatedAxis() + NoLegend()
dev.off()
print("dotplot of top10 markers: done")

# ----  Colors for condition ---- #
colors8 <- c("#faa7d6", #ALI day0
             "#87235c", #ALI day17 M04
             "#87235c", #ALI day17 M16
             "#f5a511", #ALI day17 IL4+13
             "#c95195", #ALI day3
             "#add6ff", #DM
             "#a9ded1", #EM
             "#e3ca2b") #tissue

# ----  Plot condition ---- #
pdf("UMAP_condition.pdf")
DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident",
        cols = colors8)
dev.off()

pdf("UMAP_condition_splitted.pdf", width = 15, height = 2.5)
DimPlot(dataset_combined, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident",
        cols = colors8)
dev.off()

# ----  Plot each condition in black and white with grey background ---- #
plot.list <- list()
for (i in unique(x = dataset_combined$orig.ident)) {
  plot.list[[i]] <- DimPlot(
    object = dataset_combined, 
    cells.highlight = WhichCells(object = dataset_combined, expression = orig.ident == i),
    cols = "lightgrey", cols.highlight = "black", sizes.highlight = 0.2, pt.size = 0.2
  ) + NoLegend() + ggtitle(i)
}

for (i in 1:length(plot.list)) {
  png(paste0("UMAP_conditions_BW_", i, ".png"), height = 3000, width = 3000, res = 400)
  print(plot.list[i])
  dev.off()
  }


#########################################################################
# Redefine identities
#########################################################################

# ----  Save data in a new folder ---- #
dir.create(paste0(directory, "/NewIdentities/"))
setwd(paste0(directory,"/NewIdentities/"))

# ----  Rename clusters ---- # 
new.cluster.ids <- c("basal cells", #0
                     "basal cells", #1
                     "keratinocytes", #2
                     "keratinocytes", #3
                     "keratinocytes", #4
                     "keratinocytes", #5
                     "keratinocytes", #6
                     "basal cells", #7
                     "keratinocytes", #8
                     "keratinocytes", #9
                     "keratinocytes", #10
                     "keratinocytes", #11
                     "keratinocytes", #12
                     "keratinocytes", #13
                     "basal cells", #14
                     "basal cells", #15
                     "basal cells", #16
                     "keratinocytes", #17
                     "basal cells", #18
                     "basal cells", #19
                     "basal cells", #20
                     "keratinocytes", #21
                     "keratinocytes", #22
                     "keratinocytes", #23
                     "basal cells", #24
                     "keratinocytes", #25
                     "keratinocytes", #26
                     "basal cells", #27
                     "keratinocytes", #28
                     "basal cells", #29 PLAT+ tissue-specific
                     "basal cells", #30
                     "melanocytes", #31
                     "T cells", #32
                     "basal cells", #33
                     "basal cells", #34
                     "basal cells", #35 very neuronal
                     "tuft cells", #36
                     "basal cells", #37
                     "macrophages", #38
                     "basal cells", #39
                     "keratinocytes", #40
                     "goblet cells", #41
                     "endothelial cells", #42
                     "fibroblasts", #43
                     "dendritic cells") #44

names(new.cluster.ids) <- levels(dataset_combined)
dataset_combined <- RenameIdents(dataset_combined, 
                                 new.cluster.ids)

# ----  Colors for new clusters ---- #
colors20 <- c("#EB9700", "#D93B8D", "#7C60A4", "#218A9E", "#0E5245",
              "#FFB8B8", "#A02E21", "#E02F28", "#556BAC", "#51D8E8",
              "#FF4F30", "#E856C9", "#0994E6", "#E8460C", "#FF0653",
              "#FF8A5E", "#E478FF", "#FFD536", "#E86556", "#7933E8")

# ----  Plot primary UMAP clusters ---- #
pdf("UMAP_clusters_newID.pdf") 
DimPlot(dataset_combined, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE,
        cols = colors20)
DimPlot(dataset_combined, 
        reduction = "umap", 
        cols = colors20)
dev.off()

png("UMAP_clusters_newID_labelled.png", width = 2000, height = 2000, res = 500) 
DimPlot(dataset_combined, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE,
        cols = colors20)
dev.off()

png("UMAP_clusters_newID_nolabels.png", width = 2000, height = 2000, res = 500) 
DimPlot(dataset_combined, 
        reduction = "umap",
        cols = colors20) + NoLegend()
dev.off()


#########################################################################
# Plot patient genotype
#########################################################################

# ---- UMAP with genotype ---- #
print(head(dataset_combined$genotype))
print(unique(x = dataset_combined$genotype))

# Next, switch the identity class of all cells to show replicate ID
Idents(dataset_combined) <- "genotype"

dataset_combined@active.ident <- factor(dataset_combined@active.ident, 
                                        levels=c("NA", "0", "1", "2", "3"))

png("UMAP_genotype.png", width = 2500, height = 2000, res = 500)
DimPlot(subset(dataset_combined, idents = "NA", invert = TRUE), 
        reduction = "umap", 
        pt.size = 0.1, 
        cols = c("white", "#FED1D5","#FFADB3", "#FF626C", "#CC4E56"))
dev.off()

pdf("UMAP_genotype.pdf")
DimPlot(subset(dataset_combined, idents = "NA", invert = TRUE), 
        reduction = "umap", 
        pt.size = 0.1, 
        cols = c("white", "#FED1D5","#FFADB3", "#FF626C", "#CC4E56"))
dev.off()

# ---- obtain numbers of cells per cell type ---- # 
write.csv(table(Idents(dataset_combined), dataset_combined$orig.ident), "cells_per_cluster_renamed.csv")
write.csv(prop.table(table(Idents(dataset_combined), dataset_combined$orig.ident)), "cells_per_cluster_renamed_proportions.csv")


#########################################################################
# Extract objects for Goblet, tuft, and basal cells
#########################################################################

goblet <- subset(dataset_combined, idents = "goblet cells")
saveRDS(goblet, "goblet.rds")

tuft <- subset(dataset_combined, idents = "tuft cells")
saveRDS(tuft, "tuft.rds")

basal <- subset(dataset_combined, idents = "basal cells")
saveRDS(basal, "basal.rds")

print("individual objects saved...")


#########################################################################
# Save environment
#########################################################################
save.image("envir_with_newIDs.RData")


#########################################################################
# Save downsampled object (for cross-tissue comparison)
#########################################################################
downsampled_conjunctiva <- subset(dataset_combined, downsample = 200)
saveRDS(downsampled_conjunctiva, "conjunctivaDownsampled.rds")
