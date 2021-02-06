

## Init output options
output_dir       <- paste0("results/cluster_markers/")
dir.create(output_dir, showWarnings = F)

## No need to modify
#### =====================

DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T, cols = getPalette(nClusters)) + theme_void() + theme(legend.position = "none") 
ggsave(paste0(output_dir, "latent_umap_cluster.png"), width = 8, height = 6)

p <- DotPlot(cml_seurat, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("")
ggsave(plot = p, "results/cluster_markers/dotplot_big.pdf", width = 12, height = 5)

p <- DotPlot(cml_seurat, features = (van_galen_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_vangalen.pdf", width = 9, height = 5)

celltype = Idents(cml_seurat) %>% unique() %>% extractCoarsePhenotype() 
cluster  = Idents(cml_seurat) %>% unique()
data.frame(celltype,cluster) %>% filter(celltype == "B-cell") 
data.frame(celltype,cluster) %>% filter(celltype == "CD8+") 
data.frame(celltype,cluster) %>% filter(celltype == "CD4+") 
data.frame(celltype,cluster) %>% filter(celltype == "NK") 
data.frame(celltype,cluster) %>% filter(celltype == "Monocytes") 

## Immune cells
nk_idents        <- c(1,2,12,14,15) %>% as.factor() %>% getClusterPhenotypes()
cd4_idents       <- c(3,5,8) %>% as.factor() %>% getClusterPhenotypes()
cd8_idents       <- c(0,7,11) %>% as.factor() %>% getClusterPhenotypes()
b_idents         <- c(6,10,13) %>% as.factor() %>% getClusterPhenotypes()
monocytes_idents <- c(4,9,16,17,20) %>% as.factor() %>% getClusterPhenotypes()

t_idents        <- c(as.character(cd4_idents), as.character(cd8_idents)) %>% as.factor()
tnk_cell_idents <- c(as.character(cd4_idents), as.character(cd8_idents), as.character(nk_idents)) %>% as.factor()


## Dotplots
guo_markers_genes       <- guo_markers       %>% unlist %>% unique %>% rev
van_galen_markers_genes <- van_galen_markers %>% unlist %>% unique %>% rev
zhang_cd4_markers_genes <- zhang_cd4_markers %>% unlist %>% unique %>% rev
zhang_cd8_markers_genes <- zhang_cd8_markers %>% unlist %>% unique %>% rev
nk_markers_genes        <- nk_mark           %>% unlist %>% unique %>% rev


## T/NK cells
tnk_cells <- subset(cml_seurat, idents = tnk_cell_idents)
tnk_cells <- getLatentUMAP(tnk_cells)
tnk_cells <- fixSeurat(tnk_cells)
nColors   <- length(unique(Idents(tnk_cells)))

p <- DimPlot(tnk_cells, reduction = "latent_umap", label = T, cols = getPalette(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_tnk.png", width = 12, height = 8)

p <- DotPlot(tnk_cells, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_tnk_big.pdf", width = 9, height = 5)

p <- DotPlot(tnk_cells, features = guo_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_tnk_guo.pdf", width = 9, height = 5)

p <- DotPlot(tnk_cells, features = rev(zhang_cd4_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_tnk_zheng_cd4.pdf", width = 9, height = 5)

p <- DotPlot(tnk_cells, features = rev(zhang_cd8_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_tnk_zheng_cd8.pdf", width = 9, height = 5)


tnk_markers     <- FindAllMarkers(tnk_cells, test.use = "t", verbose = T, max.cells.per.ident = 5e3, random.seed = 123)
tnk_top_markers <- all_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(tnk_markers, "results/cluster_markers/tnk_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(tnk_top_markers, "results/cluster_markers/tnk_top_markers.txt", sep = "\t", quote = F, row.names = F)
tnk_top_markers <- tnk_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(tnk_cells, features = rev(unique(tnk_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_tnk_de_markers.png", width = 9, height = 5)







## NK cells
nk_cells <- subset(cml_seurat, idents = nk_idents)
nk_cells <- getLatentUMAP(nk_cells)
nk_cells <- fixSeurat(nk_cells)

nk_cells$cluster    <- Idents(nk_cells)
nk_cells$patient    <- extractName(nk_cells$orig.ident)
nk_cells$timepoint  <- extractTimepoint(nk_cells$orig.ident)

nColors  <- length(unique(Idents(nk_cells)))

p <- nk_cells@meta.data %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(nColors)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_nk_per_patient_per_cluster.png", width = 6, height = 4)

p <- nk_cells@meta.data %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(cluster, freq, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(4)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_nk_per_cluster_per_patient.png", width = 6, height = 4)

p <- DimPlot(nk_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_nk.png", width = 12, height = 8)

p <- DotPlot(nk_cells, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_nk_big.pdf", width = 9, height = 5)

p <- DotPlot(nk_cells, features = (nk_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_nk_pferfferle.pdf", width = 9, height = 5)

nk_markers     <- FindAllMarkers(nk_cells, test.use = "t", verbose = T, max.cells.per.ident = 1e3, random.seed = 123)
nk_top_markers <- nk_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(nk_markers, "results/cluster_markers/nk_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(nk_top_markers, "results/cluster_markers/nk_top_markers.txt", sep = "\t", quote = F, row.names = F)
nk_top_markers <- nk_top_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(nk_cells, features = rev(unique(nk_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_nk_de_markers.png", width = 9, height = 5)







## cluster 5
cd4_tem_cells <- subset(cml_seurat, idents = "5 CD4+ Tem")
cd4_tem_cells <- getLatentUMAP(cd4_tem_cells)
cd4_tem_cells <- fixSeurat(cd4_tem_cells)

DimPlot(cd4_tem_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
FeaturePlot(cd4_tem_cells, reduction = "latent_umap", features = c("FOXP3", "LYZ", "CD8B", "nCount_RNA"), order = T)




## t cells
t_cells <- subset(cml_seurat, idents = t_idents)
t_cells <- getLatentUMAP(t_cells)
t_cells <- fixSeurat(t_cells)

t_cells$cluster   <- Idents(t_cells)
t_cells$patient   <- extractName(t_cells$orig.ident)
t_cells$timepoint <- extractTimepoint(t_cells$orig.ident)

nColors  <- length(unique(Idents(t_cells)))

p <- t_cells@meta.data %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(nColors)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_t_per_patient_per_cluster.png", width = 6, height = 4)

p <- t_cells@meta.data %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(cluster, freq, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(5)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_t_per_cluster_per_patient.png", width = 6, height = 4)


p <- DimPlot(t_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_t.png", width = 12, height = 8)

p <- DotPlot(t_cells, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_t_big.pdf", width = 9, height = 5)

p <- DotPlot(t_cells, features = guo_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_t_guo.pdf", width = 9, height = 5)

p <- DotPlot(t_cells, features = rev(zhang_cd4_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_t_zheng_cd4.pdf", width = 9, height = 5)

p <- DotPlot(t_cells, features = rev(zhang_cd8_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_t_zheng_cd8.pdf", width = 9, height = 5)


t_markers     <- FindAllMarkers(t_cells, test.use = "t", verbose = T, max.cells.per.ident = 1e3, random.seed = 123)
t_top_markers <- t_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(t_markers, "results/cluster_markers/t_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(t_top_markers, "results/cluster_markers/t_top_markers.txt", sep = "\t", quote = F, row.names = F)
t_top_markers <- t_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(t_cells, features = rev(unique(t_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_t_de_markers.png", width = 9, height = 5)











## cd8 cells
cd8_cells <- subset(cml_seurat, idents = cd8_idents)
cd8_cells <- getLatentUMAP(cd8_cells)
cd8_cells <- fixSeurat(cd8_cells)

cd8_cells$cluster   <- Idents(cd8_cells)
cd8_cells$patient   <- extractName(cd8_cells$orig.ident)
cd8_cells$timepoint <- extractTimepoint(cd8_cells$orig.ident)

nColors  <- length(unique(Idents(cd8_cells)))

p <- cd8_cells@meta.data %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(nColors)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_cd8_per_patient_per_cluster.png", width = 6, height = 4)

p <- cd8_cells@meta.data %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(cluster, freq, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(5)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_cd8_per_cluster_per_patient.png", width = 6, height = 4)




p <- DimPlot(cd8_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_cd8.png", width = 12, height = 8)

p <- DotPlot(cd8_cells, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd8_big.pdf", width = 9, height = 5)

p <- DotPlot(cd8_cells, features = guo_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd8_guo.pdf", width = 13, height = 5)

p <- DotPlot(cd8_cells, features = (zhang_cd4_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd8_zheng_cd4.pdf", width = 12, height = 5)

p <- DotPlot(cd8_cells, features = (zhang_cd8_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd8_zheng_cd8.pdf", width = 15, height = 5)

cd8_markers     <- FindAllMarkers(cd8_cells, test.use = "t", verbose = T, max.cells.per.ident = 1e3, random.seed = 123)
cd8_top_markers <- cd8_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(cd8_markers, "results/cluster_markers/cd8_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(cd8_top_markers, "results/cluster_markers/cd8_top_markers.txt", sep = "\t", quote = F, row.names = F)
cd8_top_markers <- cd8_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(cd8_cells, features = rev(unique(cd8_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd8_de_markers.png", width = 9, height = 5)








## cd4 cells
cd4_cells <- subset(cml_seurat, idents = cd4_idents)
cd4_cells <- getLatentUMAP(cd4_cells)
cd4_cells <- fixSeurat(cd4_cells)

cd4_cells$cluster   <- Idents(cd4_cells)
cd4_cells$patient   <- extractName(cd4_cells$orig.ident)
cd4_cells$timepoint <- extractTimepoint(cd4_cells$orig.ident)

nColors  <- length(unique(Idents(cd4_cells)))

p <- cd4_cells@meta.data %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(nColors)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_cd4_per_patient_per_cluster.png", width = 6, height = 4)

p <- cd4_cells@meta.data %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(cluster, freq, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(5)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_cd4_per_cluster_per_patient.png", width = 6, height = 4)




p <- DimPlot(cd4_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_cd4.png", width = 12, height = 8)

p <- DotPlot(cd4_cells, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd4_big.pdf", width = 9, height = 5)

p <- DotPlot(cd4_cells, features = guo_markers_genes, cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd4_guo.pdf", width = 9, height = 5)

p <- DotPlot(cd4_cells, features = rev(zhang_cd4_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd4_zheng_cd4.pdf", width = 9, height = 5)

p <- DotPlot(cd4_cells, features = rev(zhang_cd4_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd4_zheng_cd4.pdf", width = 9, height = 5)

cd4_markers     <- FindAllMarkers(cd4_cells, test.use = "t", verbose = T, max.cells.per.ident = 1e3, random.seed = 123)
cd4_top_markers <- cd4_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(cd4_markers, "results/cluster_markers/cd4_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(cd4_top_markers, "results/cluster_markers/cd4_top_markers.txt", sep = "\t", quote = F, row.names = F)
cd4_top_markers <- cd4_top_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(cd4_cells, features = rev(unique(cd4_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_cd4_de_markers.png", width = 9, height = 5)







## b cells
b_cells <- subset(cml_seurat, idents = b_idents)
b_cells <- getLatentUMAP(b_cells)
b_cells <- fixSeurat(b_cells)

b_cells$cluster   <- Idents(b_cells)
b_cells$patient   <- extractName(b_cells$orig.ident)
b_cells$timepoint <- extractTimepoint(b_cells$orig.ident)

nColors  <- length(unique(Idents(b_cells)))

p <- b_cells@meta.data %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(nColors)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_b_per_patient_per_cluster.png", width = 6, height = 4)

p <- b_cells@meta.data %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(cluster, freq, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(5)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_b_per_cluster_per_patient.png", width = 6, height = 4)



p <- DimPlot(b_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_b.png", width = 12, height = 8)

p <- DotPlot(b_cells, features = rev(big_markers), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_b_big.pdf", width = 9, height = 5)

p <- DotPlot(b_cells, features = rev(b_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_b_vangalen.pdf", width = 9, height = 5)



b_markers     <- FindAllMarkers(b_cells, test.use = "t", verbose = T, max.cells.per.ident = 1e3, random.seed = 123)
b_top_markers <- b_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(b_markers, "results/cluster_markers/b_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(b_top_markers, "results/cluster_markers/b_top_markers.txt", sep = "\t", quote = F, row.names = F)
b_top_markers <- b_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(b_cells, features = rev(unique(b_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_b_de_markers.png", width = 9, height = 5)













## monocytes cells
monocytes_cells <- subset(cml_seurat, idents = monocytes_idents)
monocytes_cells <- getLatentUMAP(monocytes_cells)
monocytes_cells <- fixSeurat(monocytes_cells)

monocytes_cells$cluster   <- Idents(monocytes_cells)
monocytes_cells$patient   <- extractName(monocytes_cells$orig.ident)
monocytes_cells$timepoint <- extractTimepoint(monocytes_cells$orig.ident)

nColors  <- length(unique(Idents(monocytes_cells)))

p <- monocytes_cells@meta.data %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(nColors)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_monocytes_per_patient_per_cluster.png", width = 6, height = 4)

p <- monocytes_cells@meta.data %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(cluster, freq, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(5)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/cluster_markers/bar_monocytes_per_cluster_per_patient.png", width = 6, height = 4)

p <- DimPlot(monocytes_cells, reduction = "latent_umap", label = T, cols = getPalette3(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_monocytes.png", width = 12, height = 8)

p <- DotPlot(monocytes_cells, features = (van_galen_markers_genes), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_monocytes_van_galen.pdf", width = 12, height = 5)

p <- DotPlot(monocytes_cells, features = rev(unique(unlist(c(undiff_mark, myeloid_mark)))), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_monocytes_van_galen_supp.pdf", width = 9, height = 5)

monocytes_markers     <- FindAllMarkers(monocytes_cells, test.use = "t", verbose = T, max.cells.per.ident = 1e3, random.seed = 123)
monocytes_top_markers <- monocytes_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
fwrite(monocytes_markers, "results/cluster_markers/monocytes_markers.txt", sep = "\t", quote = F, row.names = F)
fwrite(monocytes_top_markers, "results/cluster_markers/monocytes_top_markers.txt", sep = "\t", quote = F, row.names = F)
monocytes_top_markers <- monocytes_top_markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

p <- DotPlot(monocytes_cells, features = rev(unique(monocytes_top_markers$gene)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot = p, "results/cluster_markers/dotplot_monocytes_de_markers.png", width = 9, height = 5)

