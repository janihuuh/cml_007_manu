
dir.create("results/manuscript/", showWarnings = F)
dir.create("results/manuscript/overall/", showWarnings = F)

cml_seurat <- readRDS("results/cml_seurat_new.rds")
Idents(cml_seurat) <- Idents(cml_seurat) %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()
cml_seurat$cluster <- Idents(cml_seurat)

cml_seurat@meta.data %>% ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(23)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/overall/vln_nCount_rna.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(23)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/overall/vln_nFeature_RNA.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(23)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/overall/vln_percent.mt.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(23)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/overall/vln_percent.ribo.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, percent.cycle, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(23)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/overall/vln_percent.cycle.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, hybrid_doublet_score, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(23)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/overall/vln_hybrid_doublet_score.pdf", width = 7, height = 4)

## Plot most notable markers
big_markers           <- c("CD3E", "TRAC",            ## T cell
                           "SELL", "IL7R", "LEF1", "TCF7", "CD4", "IL2",
                           "CD8A", "CD8B", "GZMA", "PRF1", "GZMB", "GNLY", "CXCR3",     ## CD8+ cell
                           "IFNG",
                           "GZMK", "CD27", "CD28", "IL2RA", "CCL5",             ## memory
                           "FOXP3", "PDCD1", "LAG3", "HAVCR2", "CTLA4", ## Inhibitory genes
                           "NKG7", "FCGR3A","KLRG1", "KLRB1", "KLRD1", "NCAM1", ## NK-cell
                           "LYZ", "CD14", "CST3",              ## monocytes
                           "FCER1A", "CLEC10A",                         ## cDC
                           "CLEC4C", "PTPRS", "TCF4",          ## pDC
                           "MS4A1", "CD19", "IL4R",                   ## b cells
                           "TNFRSF17", "JCHAIN",       ## plasma cells
                           "MKI67" )

DimPlot(cml_seurat, reduction = "umap", cols = getPalette3(20), label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave("results/manuscript/overall/latent_umap.png", width = 7, height = 6)

DimPlot(cml_seurat, reduction = "umap", cols = getPalette3(20), label = T, repel = T, label.size = 4) + theme_classic(base_size = 20) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave("results/manuscript/overall/latent_umap2.png", width = 7, height = 6)

DotPlot(cml_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") + theme(legend.position = "top")
ggsave("results/manuscript/overall/dotplot_big_markers.pdf", width = 11, height = 7)

DotPlot(cml_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave("results/manuscript/overall/dotplot_dufva_markers.pdf", width = 9, height = 6)

guo_genes <- guo_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()

p <- DotPlot(cml_seurat, features = guo_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/overall/dotplot_guo_markers.pdf", width = 12, height = 6)


cd4_cells  <- cml_seurat@meta.data %>% filter(cluster %in% c("5 CD4", "7 CD4 Treg", "14 CD4", "15 CD4", "16 CD4 exhausted")) %>% pull(barcode)
cd4_seurat <- subset(cml_seurat, cells = cd4_cells)
zhang_cd4_genes <- zhang_cd4_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(cd4_seurat, features = zhang_cd4_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/overall/dotplot_zhang_cd4_markers.pdf", width = 16, height = 3)


cd8_cells  <- cml_seurat@meta.data %>% filter(cluster %in% c("1 CD8 EM", "2 CD8 EM", "6 CD8 IFNg", "12 CD8")) %>% pull(barcode)
cd8_seurat <- subset(cml_seurat, cells = cd8_cells)
zhang_cd8_genes <- zhang_cd8_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(cd8_seurat, features = zhang_cd4_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/overall/dotplot_zhang_cd8_markers.pdf", width = 16, height = 3)



## Add singler
singler <- readRDS('results/singler/cml_singler.rds')
singler_metadata <- data.frame(hpca_pred      = singler$singler[[1]]$SingleR.single$labels,
                               blueprint_pred = singler$singler[[2]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

## Add into Seurat
singler_metadata                  <- singler_metadata[match(colnames(cml_seurat), singler_metadata$barcode), ]
cml_seurat$singler_hpca_pred      <- singler_metadata$hpca_pred
cml_seurat$singler_blueprint_pred <- singler_metadata$blueprint_pred



p <- cml_seurat@meta.data %>% group_by(cluster, singler_hpca_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_hpca_pred, prop, fill = singler_hpca_pred, label = singler_hpca_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(20)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel()
ggsave(plot = p, "results/manuscript/overall/bar_predictions_hpca.png", width = 12, height = 8)

p <- cml_seurat@meta.data %>% group_by(cluster, singler_blueprint_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_blueprint_pred, prop, fill = singler_blueprint_pred, label = singler_blueprint_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(18)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel() + ylim(values = c(0,1))
ggsave(plot = p, "results/manuscript/overall/bar_predictions_blueprint.png", width = 12, height = 8)



clusters <- cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>% group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(`!is.na(new_clonotypes_id)` == T) %>% filter(prop > 0.4) %>% pull(clusters) %>% as.character() %>% unique()

cml_seurat@meta.data %>% 
  mutate(cluster = Idents(cml_seurat)) %>%
  filter(cluster %in% clusters) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%
  
  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/manuscript/overall/scatter_gini_clonality.pdf", width = 5, height = 4)



cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>%
  group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%
  
  ggplot(aes(reorder(clusters, prop), prop)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1))
ggsave("results/manuscript/overall/bar_tcrab_cluster.pdf", width = 5, height = 4)


cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>%
  group_by(orig.ident, clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%
  
  ggplot(aes(reorder(clusters, prop), prop)) + geom_boxplot(outlier.shape = NA, fill = "lightgrey") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1)) + geom_jitter(size = 0.5)
ggsave("results/manuscript/overall/box_tcrab_cluster.pdf", width = 5, height = 4)





cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>%
  group_by(clusters) %>% summarise(n = n()) %>%
  
  ggplot(aes(reorder(clusters, n), n, fill = clusters, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(24)) + theme(legend.position = "none") + labs(x = "", y = "nCells") +
  geom_text()
ggsave("results/manuscript/overall/bar_n_cluster.pdf", width = 5, height = 4)




## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(cml_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript/overall/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
all_markers <- fread("results/manuscript/overall/all_markers_1e3.txt")#  %>% mutate(cluster %>% extractClusterNumber() %>% as.factor() %>% getClusterPhenotypes())

set.seed(123)
top10  <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

cells <- cml_seurat@meta.data %>% mutate(cluster = Idents(cml_seurat)) %>% group_by(cluster) %>% sample_n(27) %>% pull(barcode)
p <- DoHeatmap(cml_seurat, cells = cells, features = top10$gene, angle = 90, group.colors = getPalette(24)) + theme(legend.position = "none")
ggsave(plot = p, "results/manuscript/overall/heatmap_top10_markers.png", width = 16, height = 20)



##  Per time point
cml_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% 
  ggplot(aes(timepoint, prop, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(23)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/manuscript/overall/bar_cluster.pdf", width = 7, height = 4)


cml_seurat@meta.data %>% 
  #  filter(!cluster %in% c("0 Monocytes low quality")) %>% 
  group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% 
  mutate(celltype = extractCoarsePhenotype(cluster)) %>% 
  mutate(celltype = ifelse(celltype %in% c("cycling", "Cytotoxic", "pDC", "T"), "other", celltype)) %>% 
  
  ggplot(aes(timepoint, prop, color = cluster, group = cluster)) + geom_point() + geom_path() + facet_wrap(~celltype, scales = "free_y") + scale_color_manual(values = getPalette3(23)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45) +
  facets_nice
ggsave("results/manuscript/overall/line_cluster.pdf", width = 10, height = 6)


cml_seurat$patient.x
patient_df <- cml_seurat@meta.data %>% group_by(orig.ident, timepoint, patient.x) %>% summarise(n = n()) %>% dplyr::select(-n)

cml_seurat@meta.data %>% 
  group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>% 
  ggplot(aes(timepoint,prop, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic(base_size = 12) +
  facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + 
  ggpubr::stat_compare_means(label = "p.format") + scale_fill_manual(values = getPalette3(4)) +
  theme(legend.position = "none") + facets_nice + geom_jitter(size = 0.5) + labs(x = "")
ggsave("results/manuscript/overall/box_cluster.pdf", width = 15, height = 7)

cml_seurat@meta.data %>% 
  group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>% 
  filter(timepoint %in% c("dg", "3mo")) %>% 
  ggplot(aes(timepoint,prop, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + ggpubr::stat_compare_means(label = "p.format") + scale_fill_manual(values = getPalette3(4)) +
  theme(legend.position = "none") + facets_nice + geom_jitter(size = 0.5)
ggsave("results/manuscript/overall/box_cluster_2v1.pdf", width = 12, height = 7)

cml_seurat@meta.data %>% 
  group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>% 
  filter(timepoint %in% c("3mo", "12mo")) %>% 
  ggplot(aes(timepoint,prop, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + ggpubr::stat_compare_means(label = "p.format") + scale_fill_manual(values = getPalette3(4)[c(3:4)]) +
  theme(legend.position = "none") + facets_nice + geom_jitter(size = 0.5)
ggsave("results/manuscript/overall/box_cluster_3v2.pdf", width = 12, height = 7)








df        <- cml_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
median_df <- df %>% group_by(timepoint, cluster) %>% summarise(median = median(prop)) %>% top_n(n = 5, wt = median) %>% arrange(desc(median))
# median_df <- df %>% group_by(timepoint, cluster) %>% summarise(median = median(prop)) %>% group_by(cluster) %>% top_n(n = 1, wt = median) %>% arrange(desc(median))

# umap_means_df <- umap_means_df %>% left_join(median_df)
umap_means_df <- cml_seurat@reductions$umap@cell.embeddings %>% as.data.frame %>% bind_cols(cml_seurat@meta.data) %>% group_by(cluster) %>% summarise(umap1 = median(UMAP_1), umap2 = median(UMAP_2)) 
umap_means_df <- median_df %>% left_join(umap_means_df)

cml_seurat@reductions$umap@cell.embeddings %>% as.data.frame %>% bind_cols(cml_seurat@meta.data) %>%
  #  filter(!is.na(orig.ident)) %>% filter(orig.ident != "NA") %>% mutate(orig.ident = droplevels(as.factor(orig.ident))) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = timepoint)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~timepoint, ncol = 7) + 
  scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  ggrepel::geom_text_repel(data = umap_means_df, aes(umap1,umap2,label=cluster), color = "black", size = 3.5) +
  scale_color_manual(values = getPalette3(7)) + guides(color=FALSE) + theme(legend.position = "none")
ggsave("results/manuscript/overall/umap_dens.png", width = 7, height = 4)




## Get DEGs
DEG_cluster_df <- lapply(levels(Idents(cml_seurat)), getDEGbyCluster, seurat_object = cml_seurat) %>% rbindlist()
write.table(DEG_cluster_df, "results/manuscript/overall/deg_total.txt", sep = "\t", quote = F, row.names = F)
DEG_cluster_df <- fread("results/manuscript/overall/deg_total.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.factor() %>% getClusterPhenotypes())

df <- DEG_cluster_df %>% mutate(celltype = extractCoarsePhenotype(cluster)) %>% mutate(cluster = cluster %>% reorderClusters()) %>% 
  group_by(celltype,cluster, timepoint) %>% summarise(n=n()) %>% ungroup() %>% 
  mutate(celltype = ifelse(celltype %in% c("cycling", "Cytotoxic", "pDC"), "other", celltype)) %>% 
  filter(timepoint %in% c("3movdg", "12mov3mo")) %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("3movdg", "12mov3mo")))
  
ggplot(df, aes(timepoint,n,group=cluster,color=cluster)) + geom_point() + geom_path() + facet_wrap(~celltype) + scale_color_manual(values=getPalette3(23)) + facets_nice + labs(x = "", y = "DEGs") + ggpubr::rotate_x_text(45) +
  ggrepel::geom_text_repel(data = subset(df, timepoint == "3movdg"), aes(label=cluster), nudge_y = 100) + theme(legend.position = "none")
ggsave("results/manuscript/overall/line_deg.pdf", width = 8, height = 6)


