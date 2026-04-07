library(ggplot2)
library(Seurat)
library(sctransform)
library(dplyr)
##########################################################################################
#"identify cell types"
#  map cluster numbers → cell types
#tell you exactly which cluster = beta, alpha, etc.
#guide marker gene validation
##########################################################################################
##########################################################################################
islet <- readRDS("data/processed/islet_clustered.rds")

markers <- FindAllMarkers(  islet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)  # These are marker genes used for cell-type identification
library(dplyr)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top_markers
write.csv(markers, "cluster_markers.csv", row.names = FALSE)
write.csv(top_markers, "top_markers_per_cluster.csv", row.names = FALSE)

islet_clean <- subset(islet, idents = c("1", "4", "6"))  #Keep only endocrine clusters:

#"confirm cell types with plots"
#1.Plot canonical markers on UMAP
p_markers <-FeaturePlot(islet, features = c("INS", "GCG", "SST", "PPY"), ncol = 2)
#2: Validate exocrine and stromal contamination
p_exocrine <- FeaturePlot(islet,features = c("PRSS1", "CPA1", "KRT19", "COL1A1"),ncol = 2)
#3: Quantitative validation (VlnPlot)
p_vln <- VlnPlot( islet, features = c("INS", "GCG", "SST", "PRSS1", "COL1A1"), ncol = 3)
#4: Check cluster identity labels
p_clusters <- DimPlot(islet, label = TRUE)
#5: Optional overlay (clean visualization)
p_ins <- FeaturePlot(islet,features = "INS",reduction = "umap",  cols = c("lightgrey", "red"))

#dir.create("plots/validation", recursive = TRUE, showWarnings = FALSE)

ggsave("plots/validation/endocrine.png", p_markers, width = 10, height = 8, dpi = 300)
ggsave("plots/validation/exocrine.png", p_exocrine, width = 10, height = 8, dpi = 300)
ggsave("plots/validation/vln.png", p_vln, width = 12, height = 8, dpi = 300)
ggsave("plots/validation/umap_clusters.png", p_clusters, width = 8, height = 6, dpi = 300)
ggsave("plots/validation/INS_overlay.png", p_ins, width = 6, height = 5, dpi = 300)


# Check cluster levels first

#1: Final cluster → cell type mapping
levels(islet)
new.cluster.ids <- c(
  "Alpha",        # 0
  "Beta",         # 1
  "Acinar",       # 2
  "Ductal",       # 3
  "Beta",         # 4
  "Alpha",        # 5
  "Delta",        # 6
  "Alpha-like",   # 7 (optional subtype)
  "Endocrine",    # 8 (mixed/unclear)
  "Fibroblast",   # 9
  "Beta-like",    # 10 (optional subtype)
  "PP",           # 11
  "PP-like",      # 12
  "Ductal-like",  # 13
  "PP-like",      # 14
  "Alpha-like"    # 15
)

names(new.cluster.ids) <- levels(islet)
islet <- RenameIdents(islet, new.cluster.ids)
p_final <- DimPlot(islet, label = TRUE, repel = TRUE)
ggsave("plots/validation/UMAP_celltype_annotation.png",plot = p_final, width = 10, height = 8, dpi = 300)

p_pub <- DimPlot(  islet,label = TRUE, repel = TRUE, pt.size = 1) + NoLegend()
ggsave("plots/validation/UMAP_celltype_clean_annotation.png",plot = p_pub,width = 8,height = 6,dpi = 300)

#4: Store annotation in metadata
islet$cell_type <- Idents(islet)
#5: Quick sanity check
table(islet$cell_type)

#OPTIONAL (merge subclusters into major types)
islet$cell_type_simple <- recode(
  islet$cell_type,
  "Beta-like" = "Beta",
  "Alpha-like" = "Alpha",
  "PP-like" = "PP",
  "Ductal-like" = "Ductal",
  "Endocrine" = "Endocrine"
)

saveRDS(islet, "islet_annotated.rds")
write.csv(islet@meta.data, "cell_metadata.csv")

#STEP: Filter only endocrine cells
endocrine <- subset( islet,subset = cell_type_simple %in% c("Beta", "Alpha", "Delta", "PP"))
table(endocrine$cell_type_simple)
cc1 <- DimPlot(endocrine, group.by = "cell_type_simple", label = TRUE, repel = TRUE)
ggsave("plots/validation/UMAP_celltype_endocrine_only.png",plot = cc1,width = 8,height = 6,dpi = 300)

saveRDS(endocrine, "islet_endocrine_only.rds")


endocrine <- NormalizeData(endocrine)
endocrine <- FindVariableFeatures(endocrine)
endocrine <- ScaleData(endocrine)

endocrine <- RunPCA(endocrine)
p_elbow <-ElbowPlot(endocrine)

endocrine <- RunUMAP(endocrine, dims = 1:15)
endocrine <- FindNeighbors(endocrine, dims = 1:15)
endocrine <- FindClusters(endocrine, resolution = 0.4)
p_umap <- DimPlot(endocrine, label = TRUE)
p_umap_clean <- DimPlot(endocrine, label = TRUE, repel = TRUE, pt.size = 1) + NoLegend()
p_pca <- DimPlot(endocrine, reduction = "pca")

dir.create("plots/endocrine_final", recursive = TRUE, showWarnings = FALSE)

ggsave("plots/endocrine_final/umap.png", p_umap, width = 8, height = 6, dpi = 300)
ggsave("plots/endocrine_final/umap_clean.png", p_umap_clean, width = 8, height = 6, dpi = 300)
ggsave("plots/endocrine_final/elbow.png", p_elbow, width = 6, height = 5, dpi = 300)
ggsave( "plots/endocrine_final/endocrine_pca.png", plot = p_pca, width = 8, height = 6, dpi = 300)
write.csv(endocrine@meta.data, "plots/endocrine_final/endocrine_metadata.csv")



##########################################################################################
# Part 2 "analyze endocrine clusters"
##########################################################################################


#STEP 1: Recompute markers (on endocrine-only data)
endo_markers <- FindAllMarkers( endocrine,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
library(dplyr)
top_endo <- endo_markers %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC)
top_endo

# STEP 2: Identify cluster identities again
p1 <- FeaturePlot( endocrine,features = c("INS", "IAPP", "GCG", "SST"),ncol = 2)
#STEP 3: Look for beta cell states (VERY IMPORTANT)
p2 <- FeaturePlot(endocrine,features = c("INS", "IAPP", "ERO1B", "DDIT3", "HSPA5"),ncol = 3)
#STEP 4: Check alpha cell variation
p3 <- FeaturePlot(endocrine, features = c("GCG", "ARX", "MAFB"),ncol = 2)

#STEP 5: Quantify with VlnPlot
p4 <- VlnPlot( endocrine, features = c("INS", "DDIT3", "HSPA5"), ncol = 3)
#STEP 6: Assign functional labels
new_ids <- c(
  "Beta_healthy",    # cluster 0
  "Beta_stressed",   # cluster 1
  "Alpha",           # cluster 2
  "Alpha",           # cluster 3
  "Delta",           # cluster 4
  "Beta_healthy",    # cluster 5
  "Beta_stressed",   # cluster 6
  "Alpha",           # cluster 7
  "Delta",           # cluster 8
  "Beta",            # cluster 9
  "Alpha",           # cluster 10
  "Beta",            # cluster 11
  "PP"               # cluster 12
)

names(new_ids) <- levels(endocrine)
endocrine <- RenameIdents(endocrine, new_ids)
#STEP 7: Visualize final states
p5 <- DimPlot(endocrine, label = TRUE, repel = TRUE)
p6 <- DimPlot(endocrine, label = TRUE, repel = TRUE) + NoLegend()


dir.create("endocrine_analysis", showWarnings = FALSE)
dir.create("endocrine_analysis/plots", showWarnings = FALSE)
dir.create("endocrine_analysis/data", showWarnings = FALSE)

write.csv(endo_markers, "endocrine_analysis/data/endo_markers_all.csv", row.names = FALSE)

write.csv(top_endo, "endocrine_analysis/data/top_endo_markers.csv", row.names = FALSE)

saveRDS(endocrine, "endocrine_analysis/data/endocrine_annotated.rds")

ggsave("endocrine_analysis/plots/featureplot_endocrine_markers.png", plot = p1, width = 10, height = 8, dpi = 300)
ggsave("endocrine_analysis/plots/featureplot_beta_states.png", plot = p2, width = 12, height = 8, dpi = 300)
ggsave("endocrine_analysis/plots/featureplot_alpha_states.png", plot = p3, width = 10, height = 6, dpi = 300)
ggsave("endocrine_analysis/plots/vlnplot_beta_states.png", plot = p4, width = 12, height = 8, dpi = 300)
ggsave("endocrine_analysis/plots/umap_endocrine_states.png", plot = p5, width = 10, height = 8, dpi = 300)
ggsave("endocrine_analysis/plots/umap_clean.png", plot = p6, width = 8, height = 6, dpi = 300)

