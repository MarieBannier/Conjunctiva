## ======================================================================= ##
##              First clustering of human conjunctiva datasets             ##
## ======================================================================= ##


#########################################################################
# Load libraries
#########################################################################

library(viridisLite, lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/") 
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
library(viridis, lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/") 
library(pheatmap, lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/") 
library(rlang, lib.loc="/hpc/local/CentOS7/hub_clevers/R_libs/4.1.0/") 


#directory <- "/hpc/hub_clevers/MB/conjunctiva/Conjunctiva_10X/Final_analysis/Analysis/Dim40/"
directory <- "/hpc/hub_clevers/MB/conjunctiva/Conjunctiva_10X/Final_analysis/NewIdentities/"
setwd(directory)

#########################################################################
# Load  environment
#########################################################################
#load("envir_with_Resolutions.RData")
load("envir_with_newIDs.RData")
print("Environment loaded")

write.csv(rownames(dataset_combined), "universe.csv")

#########################################################################
# Plot UMAP with each resolution showing clusters and MUC5AC expression
#########################################################################

## ----  Set resolution at 2 ---- # 
#Idents(dataset_combined) <- dataset_combined$integrated_snn_res.2
#print(head(dataset_combined$integrated_snn_res.2))

# ----  Calculate DEGs ---- #
#dataset_combined.markers <- FindAllMarkers(dataset_combined, 
#                                           only.pos = TRUE, 
#                                           min.pct = 0.25, 
#                                           logfc.threshold = 0.5)
#write.csv(dataset_combined.markers, "markers_AllClusters.csv")

## ----  50 colors palette ---- #
#colors50 <- c("#4287f5", "#F5D033", "#84C3BE", "#4287f5", "#CF3476",
#              "#FFA420", "#A98307", "#FFFF00", "#31372B", "#5D9B9B",
#              "#E6D690", "#1D1E33", "#4D5645", "#3E3B32", "#1D334A",
#              "#C2B078", "#F54021", "#D84B20", "#922B3E", "#025669",
#              "#82898F", "#C1876B", "#C51D34", "#3F888F", "#FF7514",
#              "#412227", "#9D9101", "#A65E2E", "#00BB2D", "#D36E70",
#              "#57A639", "#E7EBDA", "#1B5583", "#F5D033", "#7D8471",
#              "#E55137", "#DC9D00", "#256D7B", "#CDA434", "#015D52",
#              "#E6D690", "#78858B", "#606E8C", "#008F39", "#20603D",
#              "#F39F18", "#763C28", "#E5BE01", "#2D572C", "#49678D")

## ----  Plot primary UMAP clusters ---- #
#pdf("UMAP_clusters_res2.pdf") 
#DimPlot(dataset_combined, 
#        reduction = "umap", 
#        label = TRUE, 
#        repel = TRUE,
#        cols = colors50) +
#  ggtitle('Res_2')
#dev.off()
#
#png("UMAP_clusters_res2.png", width = 2000, height = 2000) 
#DimPlot(dataset_combined, 
#        reduction = "umap", 
#        label = TRUE, 
#        repel = TRUE,
#        cols = colors50) +
#  ggtitle('Res_2')
#dev.off()
#
## ---- Reload dataset markers ---- #
#dataset_combined.markers <- read.csv("markers_AllClusters.csv")
#print("markers loaded")
#
## ---- Identify top markers per cluster ---- # 
#dataset_combined.markers %>%
#    group_by(cluster) %>%
#    top_n(n = 10, wt = avg_logFC) -> top10
#print("top10 markers identified")
#
## ---- Plot top markers per cluster ---- # 
#pdf("Dot_top10_DEG_per_cluster.pdf", width = 60, height = 20)
#DotPlot(dataset_combined, 
#        features = unique(top10$gene), 
#        dot.scale=6, 
#        group.by = "integrated_snn_res.2") +
#  RotatedAxis()
#DotPlot(dataset_combined, 
#        features = unique(top10$gene), 
#        dot.scale=6, 
#        group.by = "integrated_snn_res.2") +
#  RotatedAxis() + NoLegend()
#dev.off()
#print("dotplot of top10 markers: done")

## ----  Colors for condition ---- #
#colors8 <- c("#faa7d6", #ALI day0
#             "#87235c", #ALI day17 M04
#             "#87235c", #ALI day17 M16
#             "#f5a511", #ALI day17 IL4+13
#             "#c95195", #ALI day3
#             "#add6ff", #DM
#             "#a9ded1", #EM
#             "#e3ca2b") #tissue
#
## ----  Plot condition ---- #
#pdf("UMAP_condition.pdf")
#DimPlot(dataset_combined, reduction = "umap", group.by = "orig.ident",
#        cols = colors8)
#dev.off()
#
#pdf("UMAP_condition_splitted.pdf", width = 15, height = 2.5)
#DimPlot(dataset_combined, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident",
#        cols = colors8)
#dev.off()
#
## ----  Plot each condition in black and white with grey background ---- #
#plot.list <- list()
#for (i in unique(x = dataset_combined$orig.ident)) {
#  plot.list[[i]] <- DimPlot(
#    object = dataset_combined, 
#    cells.highlight = WhichCells(object = dataset_combined, expression = orig.ident == i),
#    cols = "lightgrey", cols.highlight = "black", sizes.highlight = 0.2, pt.size = 0.2
#  ) + NoLegend() + ggtitle(i)
#}
#
#for (i in 1:length(plot.list)) {
#  png(paste0("UMAP_conditions_BW_", i, ".png"), height = 3000, width = 3000, res = 400)
#  print(plot.list[i])
#  dev.off()
#  }



#########################################################################
# Redefine identities
#########################################################################

# ----  Save data in a new folder ---- #
dir.create(paste0(directory, "/NewIdentities/"))
setwd(paste0(directory,"/NewIdentities/"))

## ----  Rename clusters ---- # 
#refined.cluster.ids <- c("IGFBP2+ basal cells", #0
#                         "ROBO1+ basal cells", #1
#                         "LCN2+ keratinocytes", #2
#                         "LCN2+ keratinocytes", #3
#                         "early keratinocytes", #4
#                         "LCN2+ keratinocytes", #5
#                         "LCN2+ keratinocytes", #6
#                         "IGFBP2+ basal cells", #7
#                         "early keratinocytes", #8
#                         "early keratinocytes", #9
#                         "LCN2+ keratinocytes", #10
#                         "LCN2+ keratinocytes", #11
#                         "SERPINB1+ keratinocytes", #12
#                         "early keratinocytes", #13
#                         "basal cells - S/G2 phase", #14
#                         "ROBO1+ basal cells", #15
#                         "ROBO1+ basal cells", #16
#                         "early keratinocytes", #17
#                         "ROBO1+ basal cells", #18
#                         "IGFBP2+ basal cells", #19
#                         "KRT23+ basal cells", #20
#                         "early keratinocytes", #21
#                         "early keratinocytes", #22
#                         "SERPINB1+ keratinocytes", #23
#                         "IGFBP2+ basal cells", #24
#                         "SERPINB1+ keratinocytes", #25
#                         "early keratinocytes", #26
#                         "FGFBP1+ basal cells", #27
#                         "LCN2+ keratinocytes", #28
#                         "PLAT+ basal cells?", #29
#                         "basal cells - M phase", #30
#                         "melanocytes", #31
#                         "T cells", #32
#                         "basal cells - S/G2 phase", #33
#                         "basal cells - M phase", #34
#                         "neuronal?", #35
#                         "tuft cells", #36
#                         "ROBO1+ basal cells", #37
#                         "macrophages", #38
#                         "S100A9+ basal cells", #39
#                         "LCN2+ keratinocytes", #40
#                         "goblet cells", #41
#                         "endothelial cells", #42
#                         "fibroblasts", #43
#                         "dendritic cells") #44
#
#new.cluster.ids <- c("basal cells", #0
#                     "basal cells", #1
#                     "keratinocytes", #2
#                     "keratinocytes", #3
#                     "keratinocytes", #4
#                     "keratinocytes", #5
#                     "keratinocytes", #6
#                     "basal cells", #7
#                     "keratinocytes", #8
#                     "keratinocytes", #9
#                     "keratinocytes", #10
#                     "keratinocytes", #11
#                     "keratinocytes", #12
#                     "keratinocytes", #13
#                     "basal cells", #14
#                     "basal cells", #15
#                     "basal cells", #16
#                     "keratinocytes", #17
#                     "basal cells", #18
#                     "basal cells", #19
#                     "basal cells", #20
#                     "keratinocytes", #21
#                     "keratinocytes", #22
#                     "keratinocytes", #23
#                     "basal cells", #24
#                     "keratinocytes", #25
#                     "keratinocytes", #26
#                     "basal cells", #27
#                     "keratinocytes", #28
#                     "basal cells", #29 PLAT+ tissue-specific
#                     "basal cells", #30
#                     "melanocytes", #31
#                     "T cells", #32
#                     "basal cells", #33
#                     "basal cells", #34
#                     "basal cells", #35 very neuronal
#                     "tuft cells", #36
#                     "basal cells", #37
#                     "macrophages", #38
#                     "basal cells", #39
#                     "keratinocytes", #40
#                     "goblet cells", #41
#                     "endothelial cells", #42
#                     "fibroblasts", #43
#                     "dendritic cells") #44
#
#names(new.cluster.ids) <- levels(dataset_combined)
#dataset_combined <- RenameIdents(dataset_combined, 
#                                 new.cluster.ids)

# ----  Colors for new clusters ---- #
colors20 <- c("#EB9700", "#D93B8D", "#7C60A4", "#218A9E", "#0E5245",
              "#FFB8B8", "#A02E21", "#E02F28", "#556BAC", "#51D8E8",
              "#FF4F30", "#E856C9", "#0994E6", "#E8460C", "#FF0653",
              "#FF8A5E", "#E478FF", "#FFD536", "#E86556", "#7933E8")

## ----  Plot primary UMAP clusters ---- #
#pdf("UMAP_clusters_newID.pdf") 
#DimPlot(dataset_combined, 
#        reduction = "umap", 
#        label = TRUE, 
#        repel = TRUE,
#        cols = colors20)
#DimPlot(dataset_combined, 
#        reduction = "umap", 
#        cols = colors20)
#dev.off()
#
#png("UMAP_clusters_newID_labelled.png", width = 2000, height = 2000, res = 500) 
#DimPlot(dataset_combined, 
#        reduction = "umap", 
#        label = TRUE, 
#        repel = TRUE,
#        cols = colors20)
#dev.off()
#
#png("UMAP_clusters_newID_nolabels.png", width = 2000, height = 2000, res = 500) 
#DimPlot(dataset_combined, 
#        reduction = "umap",
#        cols = colors20) + NoLegend()
#dev.off()

# ----  Calculate DEGs for new identities---- #
#dataset_combined.markers_newID <- FindAllMarkers(dataset_combined, 
#                                                 only.pos = TRUE, 
#                                                 min.pct = 0.1, 
#                                                 logfc.threshold = 0.5)
#write.csv(dataset_combined.markers_newID, "markers_AllClusters_newID.csv")

# ----  Load DEGs for new identities---- #
dataset_combined.markers_newID <- read.csv("markers_AllClusters_newID.csv")

## ---- Identify top markers per cell type ---- # 
#dataset_combined.markers_newID %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_logFC) -> top10
#
#dataset_combined.markers_newID %>%
#  group_by(cluster) %>%
#  top_n(n = 30, wt = avg_logFC) -> top30
#
#print(head(top30))
#print("top markers identified")
#
#
## ---- Plot top markers per cluster ---- # 
#pdf("Dot_top10_DEG_per_cluster.pdf", width = 30, height = 7)
#DotPlot(dataset_combined, 
#        features = unique(top10$gene), 
#        dot.scale=8) +
#  RotatedAxis()
#DotPlot(dataset_combined, 
#        features = unique(top10$gene), 
#        dot.scale=8) +
#  RotatedAxis() + NoLegend()
#dev.off()
#print("dotplot of top10 markers: done")

# ---- List of markers ---- # 
#markers <- rev(c("TP63", "KRT14", "COL17A1", #basal cells
#                 "AQP5", "LCN2", "FABP5", "MUC1", #keratinocytes
#                 "TFF1", "SPDEF", "MUC5AC", #Goblet cells
#                 "POU2F3", "SOX4", "AVIL", "TRPM5", #tuft
#                 "COL1A2", "COL3A1", "PDGFRA", #fibroblasts
#                 "MLANA", "MITF", "TYRP1", #melanocytes
#                 "CDH5", "TIE1", "CD36", #endothelial
#                 "CD7", "CD247", "RUNX3", #T cells
#                 "LIFR", "SPARCL1", "MCTP1", #DCs
#                 "CD14", "CD74", "HLA-DRA")) #macrophages
#
## ---- Reorder identities ---- #              
#dataset_combined@active.ident <- factor(dataset_combined@active.ident, 
#                                        levels=rev(c("basal cells", 
#                                                     "keratinocytes",
#                                                     "goblet cells", 
#                                                     "tuft cells", 
#                                                     "fibroblasts", 
#                                                     "melanocytes", 
#                                                     "endothelial cells", 
#                                                     "T cells", 
#                                                     "dendritic cells", 
#                                                     "macrophages")))

# ---- Plot markers per cluster ---- # 
#pdf("Dot_markers_per_cluster.pdf", height = 3.5, width = 10)
#DotPlot(dataset_combined, 
#        features = markers, 
#        assay="RNA", col.min = 0,
#        dot.scale=6) +
#  RotatedAxis() +
#  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
#DotPlot(dataset_combined, 
#        features = markers,
#        assay="RNA", col.min = 0,
#        dot.scale=6) +
#  RotatedAxis() + NoLegend() +
#  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
#dev.off()
#
## ---- UMAP expression of some markers ---- #
#dir.create(paste0(directory, "/GeneExpression/"))

#markers_UMAP <- unique(c(markers, top30$gene))
#print("markers identified")
#
#for (i in markers_UMAP) {
#  png(paste0(directory, "/GeneExpression/UMAP_",i,".png"), width = 2000, height = 2000, res = 500) 
#  print(FeaturePlot(dataset_combined, 
#                    features = i, 
#                    reduction = "umap", 
#                    min.cutoff = 0,
#                    order = TRUE))
#  dev.off()
#  print(i)
#}

#png(paste0(directory, "/GeneExpression/UMAP_POU2F3.png"), width = 2000, height = 2000, res = 500) 
#print(FeaturePlot(dataset_combined, 
#                  features = "POU2F3", 
#                  reduction = "umap", 
#                  min.cutoff = 0,
#                  coord.fixed = T,
#                  order = TRUE))
#dev.off()
#
## ---- Add tuft cell score ---- #
##based on the gut cell atlas Figure S5
#tuft_markers <- list(c("SOX4", "SPTLC2", "PTPN18", "CCSER1", "OGDHL", "KIT", "COL27A1",
#                       "FYB1", "ALOX5", "GRASP", "AVIL", "PTGS1", "RGS13", "LRPM",
#                       "PSTPIP2", "TRPM5", "POU2F3", "HTR3E", "SH2D6", "BMX", "HCK",
#                       "SH2D7", "CCDC129", "HPGDS", "GNG13", "MATK", "PIK3CG", "AZGP1", 
#                       "GRK5", "PLCG2", "RUNX1", "HOTAIRM1", "IL17RB", "ZFHX3", "HIP1R",
#                       "ITPR2", "ADGRG6", "TLE4", "MYO1B", "TPM1", "CCDC14", "RASSF6",
#                       "ATP2A3", "IL13RA1"))
#
#tuft_markers <- list(c("SOX4", "CCSER1",  "KIT", "PTGS1", "LRPM", "TRPM5", "POU2F3",
#                       "AZGP1", "GRK5", "PLCG2", "TPM1", "RASSF6"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = tuft_markers, name = "tuft_score")
#print("module tuft added...")
#
#png("TuftScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "tuft_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
#png("TuftScore_zoom.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "tuft_score1", 
#                  min.cutoff = 0, coord.fixed = T,
#                  order = TRUE, 
#                  cells = WhichCells(dataset_combined, idents = "tuft cells")) + 
#    scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add goblet cell score ---- #
#goblet_markers <- list(c("TFF1", "TFF3", "MUC5AC", "SPDEF", "S100P", "MSMB", "DEFB1"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = goblet_markers, name = "goblet_score")
#print("module goblet added...")
#
#png("GobletScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "goblet_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add basal cell score ---- #
#basal_markers <- list(c("TP63", "KRT5", "KRT14"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = basal_markers, name = "basal_score")
#print("module basal cells added...")
#
#png("BasalScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "basal_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add keratinocyte score ---- #
#keratinocyte_markers <- list(c("MUC1", "MUC16", "MUC4", "MUC20","AQP5", "LCN2"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = keratinocyte_markers, name = "keratinocyte_score")
#print("module keratinocyte added...")
#
#png("KeratinocyteScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "keratinocyte_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add DC score ---- #
#DC_markers <- list(c("LIFR", "SPARCL1", "MCTP1"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = DC_markers, name = "DC_score")
#print("module DCs added...")
#
#png("DendriticScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "DC_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add macrophage score ---- #
#macrophage_markers <- list(c("CD14", "CD74", "HLA-DRA"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = macrophage_markers, name = "macrophage_score")
#print("module macrophages added...")
#
#png("MacrophageScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "macrophage_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add T cell score ---- #
#Tcell_markers <- list(c("CD7", "CD247", "RUNX3"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = Tcell_markers, name = "Tcell_score")
#print("module T cells added...")
#
#png("TcellScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "Tcell_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add emdothelial score ---- #
#endothelial_markers <- list(c("CDH5", "TIE1", "CD36"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = endothelial_markers, name = "endothelial_score")
#print("module endothelial cells added...")
#
#png("EndothelialScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "endothelial_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add fibroblast score ---- #
#fibroblast_markers <- list(c("COL1A2", "COL3A1", "PDGFRA"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = fibroblast_markers, name = "fibroblast_score")
#print("module fibroblasts added...")
#
#png("FibroblastScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "fibroblast_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
## ---- Add melanocyte score ---- #
#melanocyte_markers <- list(c("MLANA", "MITF", "TYRP1"))
#
#dataset_combined <- AddModuleScore(dataset_combined, features = melanocyte_markers, name = "melanocyte_score")
#print("module melanocytes added...")
#
#png("MelanocyteScore.png", height = 2000, width = 2000, res = 500)
#print(FeaturePlot(dataset_combined, 
#                  features = "melanocyte_score1",  coord.fixed = T,
#                  min.cutoff = 0,
#                  order = TRUE) + scale_colour_viridis(option="rocket", direction = -1))
#dev.off()
#
#
# ---- Dot Plot of all cell types and their respective scores ---- #
#pdf("DotplotScores.pdf", height = 4, width = 6.3)
#print(DotPlot(dataset_combined, 
#              features = c("basal_score1", "keratinocyte_score1",
#                           "goblet_score1", "tuft_score1", 
#                           "fibroblast_score1", "melanocyte_score1",
#                           "endothelial_score1", "Tcell_score1",
#                           "DC_score1", "macrophage_score1"), col.min = 0) + 
#        scale_color_gradient(low = "beige", high = "darkred") +
#        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
#        RotatedAxis())
#dev.off()
#
#
## ---- Heatmap of DEG for epithelial cells ---- #
colors10 <- rev(c("#EB9700", "#D93B8D", "#7C60A4", "#218A9E", "#0E5245",
                  "#FFB8B8", "#A02E21", "#E02F28", "#556BAC", "#51D8E8"))
#
#png("Heatmap_EpitheliumDEG.png", height = 2000, width = 2000, res = 500) 
#print(DoHeatmap(subset(dataset_combined, downsample = 2000), 
#                features = dataset_combined.markers_newID[dataset_combined.markers_newID$cluster %in% 
#                           c("basal cells","keratinocytes", "goblet cells", "tuft cells"), "gene"], 
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("blue", "white", "red")) +
#        theme(text = element_text(size = 2)))
#dev.off()
#
#pdf("Heatmap_EpitheliumDEG.pdf") 
#print(DoHeatmap(subset(dataset_combined, downsample = 2000), 
#                features = dataset_combined.markers_newID[dataset_combined.markers_newID$cluster %in% 
#                           c("basal cells","keratinocytes", "goblet cells", "tuft cells"), "gene"], 
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("blue", "white", "red")) +
#        theme(text = element_text(size = 2)))
#dev.off()
#print("heatmap plotted...")


## ---- Calculate cluster averages ---- # 
#orig.levels <- levels(dataset_combined)
#Idents(dataset_combined) <- gsub(pattern = " ", replacement = "-", x = Idents(dataset_combined))
#orig.levels <- gsub(pattern = " ", replacement = "-", x = orig.levels)
#levels(dataset_combined) <- orig.levels
#
#av.exp <- AverageExpression(dataset_combined)$RNA
#message("Avg Exp calculated...")

## ---- Obtain correlation between clusters ---- # 
#cor.exp <- as.data.frame(cor(av.exp))
#cor.exp$x <- rownames(cor.exp)
#
#cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("basal-cells", 
#                                                          "keratinocytes", 
#                                                          "goblet-cells", 
#                                                          "tuft-cells", 
#                                                          "fibroblasts", 
#                                                          "melanocytes", 
#                                                          "endothelial-cells", 
#                                                          "T-cells", 
#                                                          "dendritic-cells", 
#                                                          "macrophages"))
#
#pdf("Heatmap_cluster_cor.pdf")
#print(ggplot(cor.df, aes(x, y, fill = correlation)) +
#  geom_tile() + 
#      scale_fill_continuous(limits=c(0, 1)) +
#      scale_fill_gradientn(colors = c("white", "black"), values = c(0,1))+
#   coord_fixed()+
#   theme_classic()
#)
#dev.off()
#
# ---- Tuft DEG ---- #
#tuft_DEG <- dataset_combined.markers_newID$gene[dataset_combined.markers_newID$cluster == "tuft cells" & 
#                                                dataset_combined.markers_newID$p_val_adj < 0.01 &
#                                                dataset_combined.markers_newID$avg_logFC > 0.5]
#
#pdf("DotPlot_TuftDEG_percluster.pdf", width = 20, height = 4.5) 
#print(DotPlot(dataset_combined, 
#              features = tuft_DEG) + 
#              scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
#        RotatedAxis())
#dev.off()
#
## ---- Tuft specific markers ---- #
#tuft_markers <- c("SUCNR1", "TAS1R3", "TAS1R2", "TAS2R",  "P2RY2", "GNAT3", "PLCB2", "PLCG2", 
#                  "SOX9", "COX", "GFI1B", "POU2F3", "TRPM5", "CHAT", "IL25", 
#                  "P2X", "BMX", "AVIL", "SPIB", "ALOX5", "ALOX5AP", "DCLK1", "DCLK2",
#                  "IL17RB",  "IL4R", "IL13RA1", "IL13RA2", "PTGS1", 
#                  "CD300LF", "TSLP")
#
#pdf("DotPlot_TuftMarkers_percluster.pdf", width = 10, height = 4.5) 
#print(DotPlot(dataset_combined, 
#              features = rev(tuft_markers), assay = "RNA", dot.scale = 7) + 
#              scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()                  
#
## ---- Tuft specific markers as heatmap across clusters ---- #
#av.exp_tuft <- av.exp[rownames(av.exp) %in% c(tuft_markers, tuft_DEG),]
#
#pdf("Heatmap_tuftmarkers.pdf", height = 10)
#pheatmap(av.exp_tuft, scale = "row", cellheight = 7, cellwidth = 8, fontsize_row = 7, fontsize_col = 7, cluster_rows = TRUE,
#         border_color = "black", color = colorRampPalette(c("navy", "white", "red"))(100))         
#dev.off()
##
### ---- Plot CDs, ILs, CCLs, CXCLs, across clusters ---- #
##pdf("DotPlot_ImmuneMarkers_percluster.pdf", width = 20, height = 4.5) 
##print(DotPlot(dataset_combined, 
##              features = c(grep("^CD", rownames(dataset_combined), value = T),
##                           grep("^CCL", rownames(dataset_combined), value = T),
##                           grep("^CXCL", rownames(dataset_combined), value = T),
##                           grep("^IL", rownames(dataset_combined), value = T)),
##              assay = "RNA", dot.scale = 7) + 
##              scale_color_gradient(low = "white", high = "black") +
##              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
##        RotatedAxis())
##dev.off()    
#
## ---- Goblet DEG ---- #
#goblet_DEG <- dataset_combined.markers_newID$gene[dataset_combined.markers_newID$cluster == "goblet cells" & 
#                                                  dataset_combined.markers_newID$p_val_adj < 0.01 &
#                                                  dataset_combined.markers_newID$avg_logFC > 0.5]
#
#pdf("DotPlot_GobletDEG_percluster.pdf", width = 27, height = 4.5) 
#print(DotPlot(dataset_combined, 
#              features = goblet_DEG) + 
#              scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +  
#        RotatedAxis())
#dev.off()
  

# ---- Goblet specific markers as heatmap across clusters ---- #
#av.exp_goblet <- av.exp[rownames(av.exp) %in% goblet_DEG,]
#
#pdf("Heatmap_gobletmarkers.pdf", height = 15)
#pheatmap(av.exp_goblet, scale = "row", cellheight = 7, cellwidth = 8, fontsize_row = 7, fontsize_col = 7, cluster_rows = TRUE,
#         border_color = "black", color = colorRampPalette(c("navy", "white", "red"))(100))         
#dev.off()


#########################################################################
# Effect of IL4/13 treatment
#########################################################################

## ---- Subset late ALI conditions with and without IL4/13 ---- #
#IL_dataset <- subset(dataset_combined, subset = orig.ident %in% c("ALId17_hCjM04", "ALId17_hCjM16", "ALId17_IL4+13"))
#
## ---- Create a treatment metadata column ---- #
#IL_dataset$treatment <- IL_dataset$orig.ident
#IL_dataset$treatment[IL_dataset$orig.ident %in% c("ALId17_hCjM04", "ALId17_hCjM16")] <- "EM"
#IL_dataset$treatment[IL_dataset$orig.ident == "ALId17_IL4+13"] <- "IL4.13"
#
## ---- Create a metadata column that contains both treatment and cell type information ---- #
#IL_dataset$status <- paste0(IL_dataset@active.ident, "-", IL_dataset$treatment)
#IL_dataset$CellType <- IL_dataset@active.ident
#message("Metadata added")

# ---- Obtain DEG based on IL4/13 treatment in each cell type ---- #
#message("Calculating DEGs based on IL4/13 treatment...")
#
#Idents(IL_dataset) <- IL_dataset$status
#for (i in unique(IL_dataset$CellType)) {
#  print(i)
#  IL_DEG <- FindMarkers(IL_dataset, ident.1 = paste0(i, "_IL4/13"), ident.2 = paste0(i, "_EM"))
#  write.csv(IL_DEG, paste0(i, ".csv"))
#  }
#
#
## ---- Plot Dotplot with seperate conditions ---- #
#tuft_IL <- read.csv("tuft_cells.csv", row.names = 1, header = T)
#goblet_IL <- read.csv("goblet_cells.csv", row.names = 1, header = T)
#keratinocytes_IL <- read.csv("keratinocytes.csv", row.names = 1, header = T)
#basal_IL <- read.csv("basal_cells.csv", row.names = 1, header = T)
#
#tuft_IL <- filter(tuft_IL, p_val_adj < 0.01)
#goblet_IL <- filter(goblet_IL, p_val_adj < 0.01)
#keratinocytes_IL <- filter(keratinocytes_IL, p_val_adj < 0.01)
#basal_IL <- filter(basal_IL, p_val_adj < 0.01)
#
#DEG_IL <- c(rownames(tuft_IL), 
#            rownames(goblet_IL), 
#            rownames(keratinocytes_IL), 
#            rownames(basal_IL))
#
#
#Idents(IL_dataset) <- IL_dataset$CellType
#
#pdf("ILeffect_DotPlot_allcells_p0,01.pdf", width = 23, height = 3.5)
#print(DotPlot(IL_dataset, 
#              features = unique(DEG_IL),
#              split.by = "treatment",
#              #group.by = "CellType",
#              dot.scale = 7,
#              cols = c("#B6EAD7", "#0F99B2")) + 
#              #scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()
#
#pdf("ILeffect_DotPlot_tuft_p0,01.pdf", width = 6, height = 2)
#print(DotPlot(subset(IL_dataset, idents = "tuft-cells"), 
#              features = rownames(tuft_IL),
#              split.by = "treatment",
#              #group.by = "CellType",
#              dot.scale = 7,
#              cols = c("#B6EAD7", "#0F99B2")) + 
#              #scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()
#
#pdf("ILeffect_DotPlot_goblet_p0,01.pdf", width = 11, height = 2)
#print(DotPlot(subset(IL_dataset, idents = "goblet-cells"), 
#              features = rownames(goblet_IL),
#              split.by = "treatment",
#              #group.by = "CellType",
#              dot.scale = 7,
#              cols = c("#B6EAD7", "#0F99B2")) + 
#              #scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()
#
#pdf("ILeffect_DotPlot_basal_p0,01.pdf", width = 9, height = 2)
#print(DotPlot(subset(IL_dataset, idents = "basal-cells"), 
#              features = rownames(basal_IL),
#              split.by = "treatment",
#              #group.by = "CellType",
#              dot.scale = 7,
#              cols = c("#B6EAD7", "#0F99B2")) + 
#              #scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()
#
#pdf("ILeffect_DotPlot_keratinocytes_p0,01.pdf", width = 10, height = 2)
#print(DotPlot(subset(IL_dataset, idents = "keratinocytes"), 
#              features = rownames(keratinocytes_IL),
#              split.by = "treatment",
#              #group.by = "CellType",
#              dot.scale = 7,
#              cols = c("#B6EAD7", "#0F99B2")) + 
#              #scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()
#
#tuft_IL <- filter(tuft_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5) %>% arrange(desc(avg_logFC))
#goblet_IL <- filter(goblet_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5) %>% arrange(desc(avg_logFC))
#keratinocytes_IL <- filter(keratinocytes_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5) %>% arrange(desc(avg_logFC))
#basal_IL <- filter(basal_IL, p_val_adj < 0.01, abs(avg_logFC) > 0.5)  %>% arrange(desc(avg_logFC))
#
#DEG_IL <- c(rownames(tuft_IL), 
#            rownames(goblet_IL), 
#            rownames(keratinocytes_IL), 
#            rownames(basal_IL))
#
#pdf("ILeffect_DotPlot_fc0,5.pdf", width = 11, height = 3.5)
#print(DotPlot(IL_dataset, 
#              features = unique(DEG_IL),
#              split.by = "treatment",
#              #group.by = "CellType",
#              dot.scale = 7,
#              cols = c("#B6EAD7", "#0F99B2")) + 
#              #scale_color_gradient(low = "white", high = "black") +
#              geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke = 0.5) +  
#        RotatedAxis())
#dev.off()
#
#
#av.exp_IL <- av.exp[rownames(av.exp) %in% DEG_IL,]
#
#Idents(IL_dataset) <- IL_dataset$CellType
#
#pdf("Heatmap_ILtreatment_allcells.pdf", height = 4)
#print(DoHeatmap(subset(IL_dataset, downsample = 300), 
#                features = unique(DEG_IL), 
#                group.by = "CellType",
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("#3C9AB2", "#E8C520", "#F22300")) +
#        theme(text = element_text(size = 6)))
#dev.off()
#
#pdf("Heatmap_ILtreatment_basalcells.pdf", height = 1)
#print(DoHeatmap(subset(IL_dataset, idents = "basal-cells"), 
#                features = rownames(basal_IL), 
#                group.by = "treatment",
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("#3C9AB2", "#E8C520", "#F22300")) +
#        theme(text = element_text(size = 6)))
#dev.off()
#
#pdf("Heatmap_ILtreatment_goblet.pdf", height = 2)
#print(DoHeatmap(subset(IL_dataset, idents = "goblet-cells"), 
#                features = rownames(goblet_IL), 
#                group.by = "treatment",
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("#3C9AB2", "#E8C520", "#F22300")) +
#        theme(text = element_text(size = 6)))
#dev.off()
#
#pdf("Heatmap_ILtreatment_tuft.pdf", height = 2)
#print(DoHeatmap(subset(IL_dataset, idents = "tuft-cells"), 
#                features = rownames(tuft_IL), 
#                group.by = "treatment",
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("#3C9AB2", "#E8C520", "#F22300")) +
#        theme(text = element_text(size = 6)))
#dev.off()
#
#pdf("Heatmap_ILtreatment_keratinocytes.pdf", height = 1)
#print(DoHeatmap(subset(IL_dataset, idents = "keratinocytes"), 
#                features = rownames(keratinocytes_IL), 
#                group.by = "treatment",
#                label = FALSE,
#                group.colors = colors10) + scale_fill_gradientn(colors = c("#3C9AB2", "#E8C520", "#F22300")) +
#        theme(text = element_text(size = 6)))
#dev.off()
#
## ---- Overlap with MS ---- #
#
## Load and reshape data
#MS_IL_up <- read.csv("VolcanoPlot.csv", header = TRUE, sep = ";")
#
#MS_IL_up <- MS_IL_up[!duplicated(MS_IL_up$X),]
#rownames(MS_IL_up) <- MS_IL_up$X
#MS_IL_up <- MS_IL_up[,-1]
#
## Filter gene names that are up with fc > 2
#MS_genes_IL_up <- rownames(MS_IL_up[MS_IL_up$fold.change > 2 & MS_IL_up$p < 0.05,])
#
#write.csv(MS_genes_IL_up, "MS_genes_IL_up.csv")
#
## Calculate overlap with scSeq data of all cell types
#tuft_IL <- read.csv("tuft_cells.csv", row.names = 1, header = T)
#goblet_IL <- read.csv("goblet_cells.csv", row.names = 1, header = T)
#keratinocytes_IL <- read.csv("keratinocytes.csv", row.names = 1, header = T)
#basal_IL <- read.csv("basal_cells.csv", row.names = 1, header = T)
#
#tuft_IL_up <- filter(tuft_IL, p_val_adj < 0.01,avg_logFC > 0)
#goblet_IL_up <- filter(goblet_IL, p_val_adj < 0.01,avg_logFC > 0)
#keratinocytes_IL_up <- filter(keratinocytes_IL, p_val_adj < 0.01,avg_logFC > 0)
#basal_IL_up <- filter(basal_IL, p_val_adj < 0.01,avg_logFC > 0)
#
#DEG_IL_up <- c(rownames(tuft_IL_up), 
#              rownames(goblet_IL_up), 
#              rownames(keratinocytes_IL_up), 
#              rownames(basal_IL_up))
#
#write.csv(DEG_IL_up, "DEG_IL_up.csv")
#
## ---- Plot UMAP of all DEGs ---- #
#
#for (i in DEG_IL_up) {
#  pdf(paste0("UMAP_",i,".pdf")) 
#  print(FeaturePlot(dataset_combined, 
#                    features = i, 
#                    reduction = "umap", 
#                    min.cutoff = 0,
#                    order = TRUE))
#  print(FeaturePlot(subset(dataset_combined, idents = "tuft-cells"), 
#                    features = i, 
#                    reduction = "umap", 
#                    min.cutoff = 0,
#                    order = TRUE))
#  print(FeaturePlot(subset(dataset_combined, idents = "goblet-cells"), 
#                    features = i, 
#                    reduction = "umap", 
#                    min.cutoff = 0,
#                    order = TRUE))
#  dev.off()
#  print(i)}



## ---- Cell proportions per cell type in EM vs IL4/13---- #
#write.csv(table(Idents(IL_dataset),IL_dataset$treatment), "cell_nb_ILtreatment.csv")
#write.csv(prop.table(table(Idents(IL_dataset),IL_dataset$treatment)), "cell_prop_ILtreatment.csv")
#
## ---- Cow plots of DEG upon IL4/13 ---- #
#theme_set(theme_cowplot())
#tuft_cells <- subset(IL_dataset, idents = "tuft-cells")
#Idents(tuft_cells) <- "treatment"
#avg.tuft_cells <- as.data.frame(log1p(AverageExpression(tuft_cells, verbose = FALSE)$RNA))
#avg.tuft_cells$gene <- rownames(avg.tuft_cells)
#
#goblet_cells <- subset(IL_dataset, idents = "goblet-cells")
#Idents(goblet_cells) <- "treatment"
#avg.goblet_cells <- as.data.frame(log1p(AverageExpression(goblet_cells, verbose = FALSE)$RNA))
#avg.goblet_cells$gene <- rownames(avg.goblet_cells)
#
#keratinocyte_cells <- subset(IL_dataset, idents = "keratinocytes")
#Idents(keratinocyte_cells) <- "treatment"
#avg.keratinocyte_cells <- as.data.frame(log1p(AverageExpression(keratinocyte_cells, verbose = FALSE)$RNA))
#avg.keratinocyte_cells$gene <- rownames(avg.keratinocyte_cells)
#
#basal_cells <- subset(IL_dataset, idents = "basal-cells")
#Idents(basal_cells) <- "treatment"
#avg.basal_cells <- as.data.frame(log1p(AverageExpression(basal_cells, verbose = FALSE)$RNA))
#avg.basal_cells$gene <- rownames(avg.basal_cells)
#
#p1 <- ggplot(avg.tuft_cells, aes(EM, IL4.13)) + geom_point() +
#      geom_text(hjust=0, vjust=0, label = rownames(avg.tuft_cells)) + 
#      ggtitle("tuft_cells")
##p1 <- LabelPoints(plot = p1, points = rownames(avg.tuft_cells)[avg.tuft_cells$EM > 1.5 | avg.tuft_cells$IL4.13 > 1.5] , repel = TRUE)
#p2 <- ggplot(avg.goblet_cells, aes(EM, IL4.13))  + geom_point() + 
#      geom_text(hjust=0, vjust=0, label = rownames(avg.goblet_cells)) + 
#      ggtitle("goblet_cells")
##p2 <- LabelPoints(plot = p2, points = rownames(avg.goblet_cells)[avg.goblet_cells$EM > 1.5 | avg.goblet_cells$IL4.13 > 1.5] , repel = TRUE)
#p3 <- ggplot(avg.keratinocyte_cells, aes(EM, IL4.13)) + geom_point() +
#      geom_text(hjust=0, vjust=0, label = rownames(avg.keratinocyte_cells)) + 
#      ggtitle("keratinocyte_cells")
##p3 <- LabelPoints(plot = p1, points = rownames(avg.keratinocyte_cells)[avg.keratinocyte_cells$EM > 1.5 | avg.keratinocyte_cells$IL4.13 > 1.5] , repel = TRUE)
#p4 <- ggplot(avg.basal_cells, aes(EM, IL4.13)) + geom_point() +
#      geom_text(hjust=0, vjust=0, label = rownames(avg.basal_cells)) + 
#      ggtitle("basal_cells")
##p4 <- LabelPoints(plot = p1, points = rownames(avg.basal_cells)[avg.basal_cells$EM > 1.5 | avg.basal_cells$IL4.13 > 1.5] , repel = TRUE)
#
#
#pdf("Cowplot.pdf")
#p1
#p2
#p3
#p4
#dev.off()

# ---- UMAPs of conjunctival Goblet cells markers ---- #

#conj_goblet_markers <- c("S100A11", "CLU", "F3", "AQP3", 
#                         "S100A4", "SNHG29", "KRT6A", 
#                         "KRT14", "LY6D", "LYPD2", 
#                         "KRT13", "S100A8", "S100A9")
#
#png("UMAP_conjunctiva_Goblet_markers.png", res = 600,width = 2000, height = 2000)
#FeaturePlot(dataset_combined, 
#            features = conj_goblet_markers, 
#            reduction = "umap", 
#            min.cutoff = 0,
#            order = TRUE,
#            combine = F)
#dev.off()
#
#
#
#downsampled_conjunctiva <- subset(dataset_combined, downsample = 200)
#saveRDS(downsampled_conjunctiva, "conjunctivaDownsampled.rds")

#########################################################################
# Plot cell cycle score
#########################################################################

## A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
## segregate this list into markers of G2/M phase and markers of S phase
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#
#dataset_combined <- CellCycleScoring(dataset_combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#
#pdf("UMAP_CellCycleScore.pdf")
#DimPlot(dataset_combined, reduction = "umap")
#dev.off()

#########################################################################
# Plot patient genotype
#########################################################################

# ---- UMAP with genotype ---- #
#print(head(dataset_combined$genotype))
#print(unique(x = dataset_combined$genotype))

# Next, switch the identity class of all cells to reflect replicate ID
#Idents(dataset_combined) <- "genotype"
#
#dataset_combined@active.ident <- factor(dataset_combined@active.ident, 
#                                        levels=c("NA", "0", "1", "2", "3"))
#
#png("UMAP_genotype.png", width = 2500, height = 2000, res = 500)
#DimPlot(subset(dataset_combined, idents = "NA", invert = TRUE), 
#        reduction = "umap", 
#        pt.size = 0.1, 
#        cols = c("white", "#FED1D5","#FFADB3", "#FF626C", "#CC4E56"))
#dev.off()
#
#pdf("UMAP_genotype.pdf")
#DimPlot(subset(dataset_combined, idents = "NA", invert = TRUE), 
#        reduction = "umap", 
#        pt.size = 0.1, 
#        cols = c("white", "#FED1D5","#FFADB3", "#FF626C", "#CC4E56"))
#dev.off()



## ---- obtain numbers of cells per cell type ---- # 
#write.csv(table(Idents(dataset_combined), dataset_combined$orig.ident), "cells_per_cluster_renamed.csv")
#write.csv(prop.table(table(Idents(dataset_combined), dataset_combined$orig.ident)), "cells_per_cluster_renamed_proportions.csv")


## ---- Extract objects for Goblet, tuft, and basal cells ---- # 
#goblet <- subset(dataset_combined, idents = "goblet cells")
#saveRDS(goblet, "goblet.rds")
#
#tuft <- subset(dataset_combined, idents = "tuft cells")
#saveRDS(tuft, "tuft.rds")
#
#basal <- subset(dataset_combined, idents = "basal cells")
#saveRDS(basal, "basal.rds")
#
#print("individual objects saved...")

# ---- TODO ---- #
#DONE#
#Heatmap of tuft markers in tuft cluster vs other clusters
# Compare all cell types in ALI IL4+13 vs w/o (table and then dotplot)
# Heatmap of goblet markers in goblet cluster vs other cluster


# Subcluster basal cells






#########################################################################
# Save environment
#########################################################################
#save.image("envir_with_newIDs.RData")
