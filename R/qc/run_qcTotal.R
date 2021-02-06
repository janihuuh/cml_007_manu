
cml_seurat <- readRDS("results/cml_seurat_scvi.rds")

## Detect doublets
require(scds)

#- Annotate doublet using co-expression based doublet scoring:
cml_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cml_seurat@assays$RNA@counts), colData = cml_seurat@meta.data)
cml_sce <- cxds(cml_sce)
cml_sce <- bcds(cml_sce)
cml_sce <- cxds_bcds_hybrid(cml_sce)

## Add into Seurat
cml_seurat$cxds_doublet_score   <- SingleCellExperiment::colData(cml_sce)$cxds_score
cml_seurat$bcds_doublet_score   <- SingleCellExperiment::colData(cml_sce)$bcds_score
cml_seurat$hybrid_doublet_score <- SingleCellExperiment::colData(cml_sce)$hybrid_score

dir.create("results/qc/", showWarnings = F)
qc_df <- cml_seurat@meta.data %>% as.data.frame()

###################

min_mito     <- 0
max_mito     <- 10

min_ribo     <- 5
max_ribo     <- 50

min_features <- 600
max_features <- 3000

min_counts   <- 1.75e3
max_counts   <- 10e3

min_pct50    <- 25
max_pct50    <- 60

###################

cml_seurat$temp <- cml_seurat$nCount_RNA > 2e3
DimPlot(cml_seurat, reduction = "latent_umap", group.by = "temp")

#### Violin plots
plotQcViolin(qc_df, var_to_plot = "nFeature_RNA", grouping = "orig.ident", min = min_features, max = max_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/qc/violin_nFeature_RNA.pdf", width = 5, height = 4)

plotQcViolin(qc_df, var_to_plot = "nCount_RNA", grouping = "orig.ident", min = min_counts, max = max_counts) + scale_y_log10() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/qc/violin_nCount_RNA.pdf", width = 5, height = 4)

plotQcViolin(qc_df, var_to_plot = "percent.mt", grouping = "orig.ident", min = min_mito, max = max_mito) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/qc/violin_percent_mt.pdf", width = 5, height = 4)

plotQcViolin(qc_df, var_to_plot = "percent.ribo", grouping = "orig.ident", min = min_ribo, max = max_ribo) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/qc/violin_percent_ribo.pdf", width = 5, height = 4)

plotQcViolin(qc_df, var_to_plot = "hybrid_doublet_score", grouping = "orig.ident", min = 0, max = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/qc/violin_percent_dbl.pdf", width = 5, height = 4)


## Scatter plots
p <- qc_df %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = min_counts, linetype = "dotted") +
  geom_vline(xintercept = max_counts, linetype = "dotted") +
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") + theme(legend.position = "none")
ggsave(plot = p, "results/qc/scatter_counts_vs_genes.png", width = 8, height = 6)

p <- qc_df %>%
  ggplot(aes(percent.ribo, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  scale_y_log10() +
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") +
  geom_vline(xintercept = min_ribo, linetype = "dotted") +
  geom_vline(xintercept = max_ribo, linetype = "dotted") + theme(legend.position = "none")
ggsave(plot = p, "results/qc/scatter_ribo_vs_genes.png", width = 8, height = 6)

p <- qc_df %>%
  ggplot(aes(percent.mt, nFeature_RNA, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  scale_y_log10() +
  geom_hline(yintercept = min_features, linetype = "dotted") +
  geom_hline(yintercept = max_features, linetype = "dotted") +
  geom_vline(xintercept = max_mito, linetype = "dotted") + theme(legend.position = "none")
ggsave(plot = p, "results/qc/scatter_mito_vs_genes.png", width = 8, height = 6)

p <- qc_df %>%
  ggplot(aes(percent.ribo, percent.mt, color = orig.ident)) + geom_point(size = 0.3, alpha = 0.5) + add_guide +
  geom_hline(yintercept = min_mito, linetype = "dotted") +
  geom_hline(yintercept = max_mito, linetype = "dotted") +
  geom_vline(xintercept = min_ribo, linetype = "dotted") +
  geom_vline(xintercept = max_ribo, linetype = "dotted") + theme(legend.position = "none")
ggsave(plot = p, "results/qc/scatter_ribo_vs_mito.png", width = 8, height = 6)



##############################################################################

## In total, we remove with the following conditions:
percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()
dbls_outlier         <- qc_df %>% dplyr::filter(hybrid_doublet_score   > 1)   %>% pull(barcode) %>% as.character()

##############################################################################


outlier_cells        <- c(percent_mito_outlier,
                          percent_ribo_outlier,
                          features_outlier,
                          umis_outlier,
                          dbls_outlier)

reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                          rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                          rep("features_outlier",     length(features_outlier)),
                          rep("umis_outlier",         length(umis_outlier)),
                          rep("dbls_outlier",         length(dbls_outlier)))

outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))

outlier_df %>% group_by(from,reason) %>% summarise(n = n()) %>%
  ggplot(aes(reorder(reason,n),n,fill=from,label=n)) + geom_bar(stat = "identity") + labs(x = "") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()#+ ggrepel::geom_text_repel()
ggsave("results/qc/bar_outliers.png", width = 4, height = 3)




## As percentages per sample
outlier_df <- outlier_df[!duplicated(outlier_df$barcode), ]
tot_cells <- data.frame(total_cells = table(extractName(qc_df$barcode)))
colnames(tot_cells) <- c("from", "total")
# outlier_df <- merge(outlier_df, tot_cells)

outlier_df %>% group_by(from) %>% summarise(n = n()) %>%
  left_join(tot_cells, by = "from") %>% dplyr::mutate(freq = n/total) %>% 
  ggplot(aes(from, freq, fill = from)) + geom_bar(stat = "identity") + theme_bw()
ggsave("results/qc/bar_outliers_per_sample.png", width = 4, height = 3)

write.table(outlier_df, "results/qc/outliers.txt", sep = "\t", quote = F, row.names = F)
# outlier_df <- fread("results/qc/outliers.txt")




## Remove the cells from Seurat-object and save a new seurat-object
cells.to.use <- colnames(cml_seurat)[!colnames(cml_seurat) %in% outlier_df$barcode]
cml_seurat   <- subset(cml_seurat, cells = cells.to.use)

saveRDS(cml_seurat, "results/cml_seurat_qc.rds")

DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T) + theme(legend.position = "none")
