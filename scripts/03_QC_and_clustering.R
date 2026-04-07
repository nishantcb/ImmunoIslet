library(ggplot2)
library(Seurat)
library(sctransform)
dir.create("plots", showWarnings = FALSE)

islet <- readRDS("data/processed/islet_merged.rds")
mt_genes <- grep("^MT_", rownames(islet), value = TRUE)
length(mt_genes)
# The dataset itself simply does not include mitochondrial genes Skip mitochondrial QC completely
#🎯 Use these QC metrics instead


#Sanity Check
head(rownames(islet))
tail(rownames(islet))
sum(rownames(islet) == "")
any(duplicated(rownames(islet)))
colnames(islet@meta.data)


packageVersion("Seurat")
packageVersion("SeuratObject")

islet <- NormalizeData(islet)

# ✅ Now plot QC metrics
# Create plot object
p <- VlnPlot(object = islet, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1, log = TRUE)
ggsave("QC_violin_plot.pdf", plot = p, width = 8, height = 5)
s <- FeatureScatter(islet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("scatter_plot.pdf", plot = s, width = 8, height = 5)

#QC THRESHOLDS
#islet <- subset(islet, subset = nFeature_RNA > 800 & nFeature_RNA < 4500 & nCount_RNA < 25000)
islet <- subset(islet, subset = nFeature_RNA > 800 & nFeature_RNA < 4200 & nCount_RNA > 1500 & nCount_RNA < 20000)

p1 = VlnPlot(islet, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
s1 = FeatureScatter(islet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave("QC_violin_plot_verify2.png", plot = p1, width = 8, height = 5)
ggsave("scatter_plot_verify2.png", plot = s1, width = 8, height = 5)

#######################################################################
# Normalization and Clustering
####################################################################
# SCTransform normalization
islet <- SCTransform(islet, verbose = FALSE)

#STEP 2: PCA (Dimensionality Reduction)
islet <- RunPCA(islet, verbose = FALSE)
p_elbow <- ElbowPlot(islet)
#dims = 1:20 
ggsave("plots/PCA_elbow_plot.png",plot = p_elbow,width = 6,height = 5, dpi = 300)
p_pca <- DimPlot(islet, reduction = "pca")
ggsave("plots/PCA_dimplot.png", plot = p_pca, width = 8, height = 6, dpi = 300)
p_load <- VizDimLoadings(islet, dims = 1:2, reduction = "pca")
ggsave("plots/PCA_loadings.png", plot = p_load, width = 8, height = 6, dpi = 300)

#STEP 3: UMAP (Visualization)
islet <- RunUMAP(islet, dims = 1:20)
#STEP 4: Neighbors + Clustering
islet <- FindNeighbors(islet, dims = 1:30)
islet <- FindClusters(islet, resolution = 0.5)
#STEP 5: Visualize clusters
p1 <- DimPlot(islet, reduction = "umap", label = TRUE)
p2 <- DimPlot(islet, reduction = "umap", group.by = "orig.ident")


ggsave("plots/UMAP_clusters.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("plots/UMAP_samples.png",plot = p2, width = 8, height = 6, dpi = 300)
saveRDS(islet, "data/processed/islet_clustered.rds")
##########################################################################################
##########################################################################################
##########################################################################################
# "identify cell types"
#  map cluster numbers → cell types
#tell you exactly which cluster = beta, alpha, etc.
#guide marker gene validation
##########################################################################################