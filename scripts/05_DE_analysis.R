##########################################################################################
# Part 3 Differential Expression (DE)
##########################################################################################


library(Seurat)
library(dplyr)
library(ggplot2)
#3: Read your saved Seurat object
endocrine <- readRDS("endocrine_analysis/data/endocrine_annotated.rds")
#STEP 5: Confirm labels
levels(Idents(endocrine))
table(Idents(endocrine))
#6: Run differential expression
de_beta <- FindMarkers(endocrine, ident.1 = "Beta_stressed", ident.2 = "Beta_healthy", min.pct = 0.25, logfc.threshold = 0.25)
#7: Sort and inspect
de_beta <- de_beta[order(de_beta$avg_log2FC, decreasing = TRUE), ]
head(de_beta, 15)
write.csv(de_beta, "endocrine_analysis/data/beta_state_DE.csv")

############Analysis of result##################
#dir.create("endocrine_analysis/results", showWarnings = FALSE)

ident.1 = "Beta_stressed"
ident.2 = "Beta_healthy"
#STEP 1: Filter significant genes
de_sig <- de_beta %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
write.csv(de_sig, "endocrine_analysis/results/DE_significant.csv", row.names = TRUE)
nrow(de_sig)

#STEP 2: Separate genes by direction
#Up in Beta_stressed
up_stressed <- de_sig %>% filter(avg_log2FC > 0.25) 
write.csv(up_stressed, "endocrine_analysis/results/DE_up_stressed.csv", row.names = TRUE)
#Up in Beta_healthy
up_healthy <- de_sig %>% filter(avg_log2FC < -0.25) 
write.csv(up_healthy, "endocrine_analysis/results/DE_up_healthy.csv", row.names = TRUE)

#STEP 3: Get top genes
head(up_stressed, 10); head(up_healthy, 10)

p1 = FeaturePlot(endocrine, features = c("INS", "IAPP", "ERO1B", "HSPA5", "DDIT3"), ncol = 3,  cols = c("lightgrey","red"))
ggsave("endocrine_analysis/results/FeaturePlot_beta_states.png", plot = p1, width = 12, height = 8, dpi = 300)

p_vln <- VlnPlot(endocrine, features = c("INS","ERO1B","HSPA5"), group.by = "cell_type")
ggsave("endocrine_analysis/results/VlnPlot_beta_markers.png", plot = p_vln, width = 10, height = 6, dpi = 300)

top_genes <- rownames(up_stressed)[1:20]
p_heat <- DoHeatmap(endocrine, features = top_genes, group.by = "cell_type")
ggsave("endocrine_analysis/results/Heatmap_top_stressed_genes.png", plot = p_heat, width = 12, height = 8, dpi = 300)

de_beta$gene <- rownames(de_beta)
volcano_plot <- ggplot(de_beta, aes(x = avg_log2FC, y = -log10(p_val_adj))) + geom_point(aes(color = avg_log2FC), alpha = 0.6) + scale_color_gradient2(low = "blue", mid = "grey", high = "red") + theme_classic()
ggsave("endocrine_analysis/results/Volcano_beta_states.png",  plot = volcano_plot, width = 8, height = 6,dpi = 300)

# Pathway Analysis

library(clusterProfiler); library(org.Hs.eg.db)
genes_up <- rownames(up_stressed)
ego <- enrichGO(gene = genes_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
write.csv(as.data.frame(ego), "endocrine_analysis/results/GO_enrichment.csv", row.names = FALSE)

head(ego)
ggsave("endocrine_analysis/results/GO_enrichment_plot.png", plot = dotplot(ego) + theme_classic(), width = 8, height = 6, dpi = 300)
