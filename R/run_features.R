

FeaturePlot(cml_seurat, features = c("GZMB", "B3GAT1"), order = T, label = T, repel = T, cols = c("lightgrey", "salmon"))

FeaturePlot(cd8_cells, features = c("GZMB", "B3GAT1", "CD28", "CD27"), order = T, label = T, repel = T, cols = c("lightgrey", "salmon"))

DotPlot(cd8_cells, features = c("GZMB", "B3GAT1", "CD28", "CD27")) #, order = T, label = T, repel = T, cols = c("lightgrey", "salmon"))
ggsave("results/cluster_markers/latent_umap_cd8_mature_features.png", width = 8, height = 6)
