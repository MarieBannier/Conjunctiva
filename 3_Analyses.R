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


#########################################################################
# Load  environment
#########################################################################
load("envir_with_newIDs.RData")
print("Environment loaded")


#########################################################################
# Calculate DEGs for new identities (or load)
#########################################################################
# ----  Calculate DEGs for new identities---- #
dataset_combined.markers_newID <- FindAllMarkers(dataset_combined, 
                                                 only.pos = TRUE, 
                                                 min.pct = 0.1, 
                                                 logfc.threshold = 0.5)
write.csv(dataset_combined.markers_newID, "markers_AllClusters_newID.csv")

## ----  Load DEGs for new identities---- #
#dataset_combined.markers_newID <- read.csv("markers_AllClusters_newID.csv")

# ---- Identify top markers per cell type ---- # 
dataset_combined.markers_newID %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC) -> top10

dataset_combined.markers_newID %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_logFC) -> top30

print(head(top30))
print("top markers identified")

# ---- Plot top markers per cluster ---- # 
pdf("Dot_top10_DEG_per_cluster.pdf", width = 30, height = 7)
DotPlot(dataset_combined, 
        features = unique(top10$gene), 
        dot.scale=8) +
  RotatedAxis()
DotPlot(dataset_combined, 
        features = unique(top10$gene), 
        dot.scale=8) +
  RotatedAxis() + NoLegend()
dev.off()
print("dotplot of top10 markers: done")

# ---- List of markers ---- # 
markers <- rev(c("TP63", "KRT14", "COL17A1", #basal cells
                 "AQP5", "LCN2", "FABP5", "MUC1", #keratinocytes
                 "TFF1", "SPDEF", "MUC5AC", #Goblet cells
                 "POU2F3", "SOX4", "AVIL", "TRPM5", #tuft
                 "COL1A2", "COL3A1", "PDGFRA", #fibroblasts
                 "MLANA", "MITF", "TYRP1", #melanocytes
                 "CDH5", "TIE1", "CD36", #endothelial
                 "CD7", "CD247", "RUNX3", #T cells
                 "LIFR", "SPARCL1", "MCTP1", #DCs
                 "CD14", "CD74", "HLA-DRA")) #macrophages

# ---- Reorder identities ---- #              
dataset_combined@active.ident <- factor(dataset_combined@active.ident, 
                                        levels=rev(c("basal cells", 
                                                     "keratinocytes",
                                                     "goblet cells", 
                                                     "tuft cells", 
                                                     "fibroblasts", 
                                                     "melanocytes", 
                                                     "endothelial cells", 
                                                     "T cells", 
                                                     "dendritic cells", 
                                                     "macrophages")))

# ---- Plot markers per cluster ---- # 
pdf("Dot_markers_per_cluster.pdf", height = 3.5, width = 10)
DotPlot(dataset_combined, 
        features = markers, 
        assay="RNA", col.min = 0,
        dot.scale=6) +
  RotatedAxis() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
DotPlot(dataset_combined, 
        features = markers,
        assay="RNA", col.min = 0,
        dot.scale=6) +
  RotatedAxis() + NoLegend() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()

# ---- UMAP expression of some markers ---- #
dir.create(paste0(directory, "/GeneExpression/"))

markers_UMAP <- unique(c(markers, top30$gene))
print("markers identified")

for (i in markers_UMAP) {
  png(paste0(directory, "/GeneExpression/UMAP_",i,".png"), width = 2000, height = 2000, res = 500) 
  print(FeaturePlot(dataset_combined, 
                    features = i, 
                    reduction = "umap", 
                    min.cutoff = 0,
                    order = TRUE))
  dev.off()
  print(i)
}

# ---- Heatmap of DEG for epithelial cells ---- #
colors10 <- rev(c("#EB9700", "#D93B8D", "#7C60A4", "#218A9E", "#0E5245",
                  "#FFB8B8", "#A02E21", "#E02F28", "#556BAC", "#51D8E8"))

png("Heatmap_EpitheliumDEG.png", height = 2000, width = 2000, res = 500) 
print(DoHeatmap(subset(dataset_combined, downsample = 2000), 
                features = dataset_combined.markers_newID[dataset_combined.markers_newID$cluster %in% 
                           c("basal cells","keratinocytes", "goblet cells", "tuft cells"), "gene"], 
                label = FALSE,
                group.colors = colors10) + scale_fill_gradientn(colors = c("blue", "white", "red")) +
        theme(text = element_text(size = 2)))
dev.off()

pdf("Heatmap_EpitheliumDEG.pdf") 
print(DoHeatmap(subset(dataset_combined, downsample = 2000), 
                features = dataset_combined.markers_newID[dataset_combined.markers_newID$cluster %in% 
                           c("basal cells","keratinocytes", "goblet cells", "tuft cells"), "gene"], 
                label = FALSE,
                group.colors = colors10) + scale_fill_gradientn(colors = c("blue", "white", "red")) +
        theme(text = element_text(size = 2)))
dev.off()
print("heatmap plotted...")


#########################################################################
# Tuft cell markers
#########################################################################

# ---- Tuft DEG ---- #
tuft_DEG <- dataset_combined.markers_newID$gene[dataset_combined.markers_newID$cluster == "tuft cells" & 
                                                dataset_combined.markers_newID$p_val_adj < 0.01 &
                                                dataset_combined.markers_newID$avg_logFC > 0.5]

pdf("DotPlot_TuftDEG_percluster.pdf", width = 20, height = 4.5) 
print(DotPlot(dataset_combined, 
              features = tuft_DEG) + 
              scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
        RotatedAxis())
dev.off()

# ---- Tuft specific markers ---- #
tuft_markers <- c("SUCNR1", "TAS1R3", "TAS1R2", "TAS2R",  "P2RY2", "GNAT3", "PLCB2", "PLCG2", 
                  "SOX9", "COX", "GFI1B", "POU2F3", "TRPM5", "CHAT", "IL25", 
                  "P2X", "BMX", "AVIL", "SPIB", "ALOX5", "ALOX5AP", "DCLK1", "DCLK2",
                  "IL17RB",  "IL4R", "IL13RA1", "IL13RA2", "PTGS1", 
                  "CD300LF", "TSLP")

pdf("DotPlot_TuftMarkers_percluster.pdf", width = 10, height = 4.5) 
print(DotPlot(dataset_combined, 
              features = rev(tuft_markers), assay = "RNA", dot.scale = 7) + 
              scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()                  

#########################################################################
# Goblet cell markers
#########################################################################

# ---- Goblet DEG ---- #
goblet_DEG <- dataset_combined.markers_newID$gene[dataset_combined.markers_newID$cluster == "goblet cells" & 
                                                  dataset_combined.markers_newID$p_val_adj < 0.01 &
                                                  dataset_combined.markers_newID$avg_logFC > 0.5]

pdf("DotPlot_GobletDEG_percluster.pdf", width = 27, height = 4.5) 
print(DotPlot(dataset_combined, 
              features = goblet_DEG) + 
              scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
        RotatedAxis())
dev.off()


#########################################################################
# Effect of IL4/13 treatment
#########################################################################

# ---- Subset late ALI conditions with and without IL4/13 ---- #
IL_dataset <- subset(dataset_combined, subset = orig.ident %in% c("ALId17_hCjM04", "ALId17_hCjM16", "ALId17_IL4+13"))

# ---- Create a treatment metadata column ---- #
IL_dataset$treatment <- IL_dataset$orig.ident
IL_dataset$treatment[IL_dataset$orig.ident %in% c("ALId17_hCjM04", "ALId17_hCjM16")] <- "EM"
IL_dataset$treatment[IL_dataset$orig.ident == "ALId17_IL4+13"] <- "IL4.13"

# ---- Create a metadata column that contains both treatment and cell type information ---- #
IL_dataset$status <- paste0(IL_dataset@active.ident, "-", IL_dataset$treatment)
IL_dataset$CellType <- IL_dataset@active.ident
message("Metadata added")

# ---- Obtain DEG based on IL4/13 treatment in each cell type ---- #
message("Calculating DEGs based on IL4/13 treatment...")
Idents(IL_dataset) <- IL_dataset$status
for (i in unique(IL_dataset$CellType)) {
  print(i)
  IL_DEG <- FindMarkers(IL_dataset, ident.1 = paste0(i, "_IL4/13"), ident.2 = paste0(i, "_EM"))
  write.csv(IL_DEG, paste0(i, ".csv"))
  }

# ---- Plot Dotplot with seperate conditions ---- #
tuft_IL <- read.csv("tuft_cells.csv", row.names = 1, header = T)
goblet_IL <- read.csv("goblet_cells.csv", row.names = 1, header = T)
keratinocytes_IL <- read.csv("keratinocytes.csv", row.names = 1, header = T)
basal_IL <- read.csv("basal_cells.csv", row.names = 1, header = T)

tuft_IL <- filter(tuft_IL, p_val_adj < 0.01)
goblet_IL <- filter(goblet_IL, p_val_adj < 0.01)
keratinocytes_IL <- filter(keratinocytes_IL, p_val_adj < 0.01)
basal_IL <- filter(basal_IL, p_val_adj < 0.01)

DEG_IL <- c(rownames(tuft_IL), 
            rownames(goblet_IL), 
            rownames(keratinocytes_IL), 
            rownames(basal_IL))

Idents(IL_dataset) <- IL_dataset$CellType

pdf("ILeffect_DotPlot_allcells_p0,01.pdf", width = 23, height = 3.5)
print(DotPlot(IL_dataset, 
              features = unique(DEG_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              cols = c("#B6EAD7", "#0F99B2")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_tuft_p0,01.pdf", width = 6, height = 2)
print(DotPlot(subset(IL_dataset, idents = "tuft-cells"), 
              features = rownames(tuft_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              cols = c("#B6EAD7", "#0F99B2")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_goblet_p0,01.pdf", width = 11, height = 2)
print(DotPlot(subset(IL_dataset, idents = "goblet-cells"), 
              features = rownames(goblet_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              cols = c("#B6EAD7", "#0F99B2")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_basal_p0,01.pdf", width = 9, height = 2)
print(DotPlot(subset(IL_dataset, idents = "basal-cells"), 
              features = rownames(basal_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              cols = c("#B6EAD7", "#0F99B2")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_keratinocytes_p0,01.pdf", width = 10, height = 2)
print(DotPlot(subset(IL_dataset, idents = "keratinocytes"), 
              features = rownames(keratinocytes_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              cols = c("#B6EAD7", "#0F99B2")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

tuft_IL <- filter(tuft_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5) %>% arrange(desc(avg_logFC))
goblet_IL <- filter(goblet_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5) %>% arrange(desc(avg_logFC))
keratinocytes_IL <- filter(keratinocytes_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5) %>% arrange(desc(avg_logFC))
basal_IL <- filter(basal_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5)  %>% arrange(desc(avg_logFC))

DEG_IL <- c(rownames(tuft_IL), 
            rownames(goblet_IL), 
            rownames(keratinocytes_IL), 
            rownames(basal_IL))

pdf("ILeffect_DotPlot_fc0,5.pdf", width = 11, height = 3.5)
print(DotPlot(IL_dataset, 
              features = unique(DEG_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              cols = c("#B6EAD7", "#0F99B2")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_tuft_p0,01_fc05.pdf", width = 5, height = 2)
print(DotPlot(subset(IL_dataset, idents = "tuft-cells"), 
              features = rownames(tuft_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              scale.min = 0, scale.max= 100,
              #cols = c("#B6EAD7", "#0F99B2")) + 
              cols = c("black", "black")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_goblet_p0,01_fc05.pdf", width = 10, height = 2)
print(DotPlot(subset(IL_dataset, idents = "goblet-cells"), 
              features = rownames(goblet_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              scale.min = 0, scale.max= 100,
              #cols = c("#B6EAD7", "#0F99B2")) + 
              cols = c("black", "black")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_basal_p0,01_fc05.pdf", width = 5, height = 2)
print(DotPlot(subset(IL_dataset, idents = "basal-cells"), 
              features = rownames(basal_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              scale.min = 0, scale.max= 100,
              #cols = c("#B6EAD7", "#0F99B2")) + 
              cols = c("black", "black")) + 
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_DotPlot_keratinocytes_p0,01_fc05.pdf", width = 5, height = 2)
print(DotPlot(subset(IL_dataset, idents = "keratinocytes"), 
              features = rownames(keratinocytes_IL),
              split.by = "treatment",
              #group.by = "CellType",
              dot.scale = 7,
              scale.min = 0, scale.max= 100,
              #cols = c("#B6EAD7", "#0F99B2")) + 
              cols = c("black", "black")) +  
              #scale_color_gradient(low = "white", high = "black") +
              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
        RotatedAxis())
dev.off()

pdf("ILeffect_Vln_goblet_0,01_0,5_group.pdf", width = 3, height = 3)
VlnPlot(object = subset(IL_dataset, idents = "goblet-cells"), 
              features = rownames(goblet_IL),
              group.by = "treatment",
              cols = c("#B6EAD7", "#0F99B2"),
              combine = FALSE,
              assay = "RNA", pt.size = 1)
dev.off()

pdf("ILeffect_Vln_keratinocytes_0,01_0,5_group.pdf", width = 3, height = 3)
VlnPlot(object = subset(IL_dataset, idents = "keratinocytes"), 
              features = rownames(keratinocytes_IL),
              pt.size = 0.05, 
              group.by = "treatment",
              cols = c("#B6EAD7", "#0F99B2"),
              combine = FALSE,
              assay = "RNA")
dev.off()

pdf("ILeffect_Vln_basal_0,01_0,5_group.pdf", width = 3, height = 3)
VlnPlot(object = subset(IL_dataset, idents = "basal-cells"), 
              features = rownames(basal_IL),
              group.by = "treatment",
              cols = c("#B6EAD7", "#0F99B2"),
              combine = FALSE,
              assay = "RNA", pt.size = 0.05)
dev.off()

pdf("ILeffect_Vln_tuft_0,01_0,5_group.pdf", width = 3, height = 3)
VlnPlot(object = subset(IL_dataset, idents = "tuft-cells"), 
              features = rownames(tuft_IL),
              group.by = "treatment",
              cols = c("#B6EAD7", "#0F99B2"),
              combine = FALSE,
              assay = "RNA", pt.size = 1)
dev.off()


#########################################################################
# Plot UMAP of selected gene expression
#########################################################################

# ---- UMAPs of DEGs ---- #
for (i in DEG_IL_up) {
  pdf(paste0("UMAP_",i,".pdf")) 
  print(FeaturePlot(dataset_combined, 
                    features = i, 
                    reduction = "umap", 
                    min.cutoff = 0,
                    order = TRUE))
  print(FeaturePlot(subset(dataset_combined, idents = "tuft-cells"), 
                    features = i, 
                    reduction = "umap", 
                    min.cutoff = 0,
                    order = TRUE))
  print(FeaturePlot(subset(dataset_combined, idents = "goblet-cells"), 
                    features = i, 
                   reduction = "umap", 
                    min.cutoff = 0,
                    order = TRUE))
  dev.off()
  print(i)}

# ---- UMAPs of conjunctival Goblet cells markers ---- #
conj_goblet_markers <- c("S100A11", "CLU", "F3", "AQP3", 
                         "S100A4", "SNHG29", "KRT6A", 
                         "KRT14", "LY6D", "LYPD2", 
                         "KRT13", "S100A8", "S100A9")

png("UMAP_conjunctiva_Goblet_markers.png", res = 600,width = 2000, height = 2000)
FeaturePlot(dataset_combined, 
            features = conj_goblet_markers, 
            reduction = "umap", 
            min.cutoff = 0,
            order = TRUE,
            combine = F)
dev.off()



