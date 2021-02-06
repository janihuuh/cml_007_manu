
## Init
output_dir       <- paste0("results/dimred/")
dir.create(output_dir, showWarnings = F)


## No need to modify
#### =====================

## Init visualisation object
viz_df <- cbind(cml_seurat@meta.data,
                # cml_seurat[["pca"]]@cell.embeddings[,1:6],
                cml_seurat[["latent_umap"]]@cell.embeddings,
                cluster = Idents(cml_seurat)) # %>% dplyr::rename(latent_umap_1 = latent_1, latent_umap_2 = latent_2)
write.table(viz_df, paste0(output_dir, "viz_df.txt"), sep = "\t", quote = F, row.names = F)

latent_umap_mean <- data.frame(aggregate(latent_umap_1 ~ cluster, viz_df, median), latent_umap_2 = aggregate(latent_umap_2 ~ cluster, viz_df, median)[,2])
nClusters <- levels(viz_df$cluster) %>% length


### ===== Visualize PCA

ElbowPlot(cml_seurat, ndims = 50)
ggsave(paste0(output_dir, "pca_elbow.png"), width = 6, height = 4)

VizDimLoadings(cml_seurat, dims = 1:6, reduction = "pca", col = "salmon", balanced = T, ncol = 3, nfeatures = 10)
ggsave(paste0(output_dir, "pca_loadings.png"), width = 14, height = 6)


### ===== Visualize latent_umap

DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T, cols = getPalette(nClusters)) + theme(legend.position = "none")
ggsave(paste0(output_dir, "latent_umap_cluster.png"), width = 8, height = 6)

# viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = cluster)) + geom_point(size = 0.3) + scale_color_manual(values = getPalette(nClusters)) + add_guide
# ggsave(paste0(output_dir, "latent_umap_cluster.png"), width = 8, height = 6)

ggplot() + 
  geom_point(data = subset(viz_df, is.na(new_clonotypes_id)), aes(latent_umap_1, latent_umap_2), size = 0.3, alpha = 0.5, color = "lightgrey") +
  geom_point(data = subset(viz_df, !is.na(new_clonotypes_id)), aes(latent_umap_1, latent_umap_2), size = 0.3, alpha = 0.8, color = "salmon") +
  labs(color = "TCRab")
ggsave(paste0(output_dir, "latent_umap_tcarb.png"), width = 8, height = 6)

viz_df %>%
  ggplot(aes(latent_umap_1, latent_umap_2, color = ifelse(percent.mt < 5 & nFeature_RNA < 500, "suspect", "nonsuspect"))) + geom_point(size = 0.3, alpha = 0.5) +
  labs(color = "percent.mt < 5 & nFeature_RNA < 500") + add_guide
ggsave(paste0(output_dir, "latent_umap_mito_small.png"), width = 8, height = 6)


## Viz QC
viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = log10(nFeature_RNA))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "latent_umap_nFeature.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = ifelse(nFeature_RNA < 500, "suspect", "nonsuspect"))) + geom_point(size = 0.3, alpha = 0.5) + labs(color = "nFeature_RNA < 500") + add_guide
ggsave(paste0(output_dir, "latent_umap_nFeature500.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = log10(nCount_RNA))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "latent_umap_nCount.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = ifelse(nCount_RNA < 1000, "suspect", "nonsuspect"))) + geom_point(size = 0.3, alpha = 0.5) + labs(color = "nCount_RNA < 1000") + add_guide
ggsave(paste0(output_dir, "latent_umap_nCount1000.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = percent.mt)) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "latent_umap_percmito.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = percent.ribo)) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "latent_umap_percribo.png"), width = 8, height = 6)

viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = log10(percent.cycle))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(paste0(output_dir, "latent_umap_perccycle.png"), width = 8, height = 6)

## Plot individual clusters
p <- NULL
i <- 1
colors = getPalette(nClusters)

for(cluster_temp in levels(viz_df$cluster)){
  
  p[[i]] <- ggplot() +
    geom_point(data = subset(viz_df, cluster != cluster_temp), aes(x = latent_umap_1, y = latent_umap_2), color = "lightgrey", size = 0.8) +
    geom_point(data = subset(viz_df, cluster == cluster_temp), aes(x = latent_umap_1, y = latent_umap_2), color = colors[i], size = 0.8) +
    theme_void() + theme(legend.position = "none")
  
  i <- i + 1
}

png(paste0(output_dir, "latent_umap_per_cluster.png"), width = 8, height = 6)
p
dev.off()















## How much cells with TCRab per cluster?
viz_df %>% group_by(cluster) %>% summarise(has_tcrab = sum(!is.na(new_clonotypes_id)), n = n()) %>% mutate(freq = has_tcrab / n) %>%
  ggplot(aes(cluster, freq, fill = cluster)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(nClusters)) +
  ylim(values = c(0,1)) + theme_bw() + theme(legend.position = "none") + labs(y = "cells with TCRa | TCRb | TCRab") + 
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_hline(yintercept = 0.25, linetype = "dotted")
ggsave(paste0(output_dir, "bar_cluster_with_tcrab.png"), width = 8, height = 6)


## How much cells per cluster?
p <- melt(table(Idents(cml_seurat))) %>%
  ggplot(aes(as.factor(reorder(Var1, value)), value, label = value, fill = as.factor(Var1), label = value)) +
  geom_bar(stat = "identity", color = "lightgrey") + geom_text() +
  scale_fill_manual(values = getPalette(nClusters)) +
  theme_minimal() + theme_bw() + theme(legend.position = "none") + coord_flip() + labs(x = "", y = "nCells")
ggsave(plot = p, paste0(output_dir, "bar_cluster.png"), width = 8, height = 8)


colnames(viz_df) <- colnames(viz_df) %>% make.unique()

## Qc on clusters?
p <- viz_df %>%
  ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + scale_y_log10() + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_nCount_RNA.png"), width = 16, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + scale_y_log10() + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_nFeature_RNA.png"), width = 16, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_percent_mt.png"), width = 16, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_percent.ribo.png"), width = 16, height = 6)

p <- viz_df %>%
  ggplot(aes(cluster, percent.cycle, fill = cluster)) + geom_violin(alpha = 0.3) + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "violin_percent.cycle.png"), width = 16, height = 6)


## How patient specific are the clusters?
p <- viz_df %>% group_by(patient, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(patient), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_project_per_cluster.png"), width = 6, height = 4)

p <- viz_df %>% group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(patient))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_project.png"), width = 12, height = 4)

p <- viz_df %>% group_by(cluster, effusion) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(effusion))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted")
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_effusion.png"), width = 12, height = 4)

p <- viz_df %>% group_by(effusion, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(effusion), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_effusion_per_cluster.png"), width = 6, height = 4)

p <- viz_df %>% group_by(cluster, patient, timepoint) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(timepoint))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + facet_wrap(~patient, ncol=1) + facets_nice
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_patient_per_time.png"), width = 6, height = 12)


p <- viz_df %>% 
  filter(timepoint == "dg") %>%   
  group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(patient))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_project_baseline.png"), width = 12, height = 4)

p <- viz_df %>% 
  filter(timepoint == "dg") %>%   
  group_by(cluster, effusion) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(effusion))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted")
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_effusion_baseline.png"), width = 12, height = 4)



## Per time
viz_df %>% 
  group_by(patient, timepoint, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(timepoint, freq, fill = cluster, group = cluster)) + geom_point(shape = 21, size = 3) + geom_path() + facet_wrap(~patient, scales = "free") + theme_bw() +
  scale_fill_manual(values = getPalette(nClusters)) + facets_nice
ggsave(paste0(output_dir, "line_clusters_per_timepoint.png"), width = 8, height = 4)



viz_df %>% calculateFoldchange() %>% filter(variable == "X1") %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_2v1_all.png"), width = 8, height = 4)

viz_df %>% filter(effusion == "effusion") %>% calculateFoldchange() %>% filter(variable == "X1") %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_2v1_effusion.png"), width = 8, height = 4)

viz_df %>% filter(effusion == "no effusion") %>% calculateFoldchange() %>% filter(variable == "X1") %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_2v1_no_effusion.png"), width = 8, height = 4)




viz_df %>% calculateFoldchange() %>% filter(variable == "X3") %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_3v2_all.png"), width = 8, height = 4)

viz_df %>% filter(effusion == "effusion") %>% calculateFoldchange() %>% filter(variable == "X3") %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_3v2_effusion.png"), width = 8, height = 4)

viz_df %>% filter(effusion == "no effusion") %>% calculateFoldchange() %>% filter(variable == "X3") %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_3v2_no_effusion.png"), width = 8, height = 4)






df_temp <- viz_df %>% filter(timepoint == "dg") %>% 
  tidyr::complete(cluster, effusion, fill = list(z = 0)) %>%
  
  group_by(timepoint, cluster, effusion) %>%
  summarise(n = n()) %>% mutate(prop = n / sum(n))

df_temp1 <- df_temp %>% filter(effusion == "effusion")
df_temp2 <- df_temp %>% filter(effusion == "no effusion")

df_temp <- data.frame(cluster = df_temp1$cluster,
                      "1" = log2(df_temp2$prop / df_temp1$prop)) %>% dplyr::rename("value" = "X1")

df_temp %>% plotBarFoldchange()
ggsave(paste0(output_dir, "bar_folchange_baseline_effusion.png"), width = 8, height = 4)

