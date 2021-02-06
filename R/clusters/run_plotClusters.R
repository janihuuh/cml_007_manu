


aml_seurat$from <- extractName(aml_seurat$orig.ident)

aml_seurat$project <- aml_seurat$orig.ident %>% as.character()
aml_seurat$project[extractName(aml_seurat$orig.ident) %in% c("001", "002", "008", "009", "012", "014", "017")] <- "tim3 aml"
aml_seurat$project[extractName(aml_seurat$orig.ident) %in% c("FHRB3667", "FHRB5897", "FHRB6386")] <- "hemap aml"
aml_seurat$project[extractName(aml_seurat$orig.ident) == "aa"] <- "aplastic anemia"
aml_seurat$project[extractName(aml_seurat$orig.ident) == "petti"] <- "petti aml"
aml_seurat$project[aml_seurat$orig.ident %in% c("petti_Normal_sorted_170607", "petti_Normal_sorted_170531")] <- "petti normal"

aml_seurat$celltype <- Idents(aml_seurat) %>% extractCoarsePhenotype()


## No need to modify
#### =====================

## Init output options
output_dir <- "results/clusters/"
dir.create(output_dir, showWarnings = F)

## Init viz object
viz_df <- cbind(aml_seurat@meta.data,
                aml_seurat[["latent_umap"]]@cell.embeddings,
                cluster = Idents(aml_seurat)) #%>% dplyr::rename(latent_umap_1 = latent_1, latent_umap_2 = latent_2) 
fwrite(viz_df, paste0(output_dir, "viz_df.txt"), sep = "\t", quote = F, row.names = F)

latent_umap_mean <- data.frame(aggregate(latent_umap_1 ~ cluster, viz_df, median), latent_umap_2 = aggregate(latent_umap_2 ~ cluster, viz_df, median)[,2])
nClusters <- levels(viz_df$cluster) %>% length

### ===== Visualize QC on latent_umap

p <- viz_df %>%
  ggplot(aes(latent_umap_1, latent_umap_2, color = as.factor(timepoint))) + geom_point(size = 0.3) + scale_color_manual(values = brewer.pal(5, "Pastel1")) +
  add_guide + theme_void() + labs(fill = "time point") + labs(fill = "time point")
ggsave(plot = p, paste0(output_dir, "latent_umap_timepoint.png"), width = 8, height = 6)

p <- viz_df %>% mutate(prepost = ifelse(timepoint == 1, "pre", "post")) %>%
  ggplot(aes(latent_umap_1, latent_umap_2, color = prepost)) + geom_point(size = 0.3) + scale_color_manual(values = brewer.pal(5, "Pastel1")) +
  add_guide + theme_void() + labs(fill = "time point")
ggsave(plot = p, paste0(output_dir, "latent_umap_prepost.png"), width = 8, height = 6)

p <- DimPlot(aml_seurat, reduction = "latent_umap", cols = getPalette(nClusters), label = T)
ggsave(plot = p, paste0(output_dir, "latent_umap_cluster.png"), width = 14, height = 12)

p <- ggplot() +
  geom_point(data = filter(viz_df, is.na(new_clonotypes_id)), aes(latent_umap_1, latent_umap_2), color = "lightgrey", size = 0.3, alpha = 0.5) +
  geom_point(data = filter(viz_df, !is.na(new_clonotypes_id)), aes(latent_umap_1, latent_umap_2), color = "salmon", size = 0.3, alpha = 0.5)
ggsave(plot = p, paste0(output_dir, "latent_umap_tcarb.png"), width = 8, height = 6)


## Viz QC
p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = log10(nFeature_RNA))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(plot = p, paste0(output_dir, "latent_umap_nFeature.png"), width = 8, height = 6)

p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = ifelse(nFeature_RNA < 500, "suspect", "nonsuspect"))) + geom_point(size = 0.3, alpha = 0.5) + labs(color = "nFeature_RNA < 500") + add_guide
ggsave(plot = p, paste0(output_dir, "latent_umap_nFeature500.png"), width = 8, height = 6)

p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = log10(nCount_RNA))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(plot = p, paste0(output_dir, "latent_umap_nCount.png"), width = 8, height = 6)

p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = ifelse(nCount_RNA < 1000, "suspect", "nonsuspect"))) + geom_point(size = 0.3, alpha = 0.5) + labs(color = "nCount_RNA < 1000") + add_guide
ggsave(plot = p, paste0(output_dir, "latent_umap_nCount1000.png"), width = 8, height = 6)

p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = percent.mt)) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(plot = p, paste0(output_dir, "latent_umap_percmito.png"), width = 8, height = 6)

p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = percent.ribo)) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(plot = p, paste0(output_dir, "latent_umap_percribo.png"), width = 8, height = 6)

p <- viz_df %>% ggplot(aes(latent_umap_1, latent_umap_2, color = log10(percent.cycle))) + geom_point(size = 0.3, alpha = 0.5) + scale_color_viridis_c()
ggsave(plot = p, paste0(output_dir, "latent_umap_perccycle.png"), width = 8, height = 6)




## Plot individual clusters
p <- NULL
i <- 1
colors = getPalette(length(unique(viz_df$cluster)))

for(cluster_temp in levels(viz_df$cluster)){
  
  message(cluster_temp)
  
  p[[i]] <- ggplot() +
    geom_point(data = subset(viz_df, cluster != cluster_temp), aes(x = latent_umap_1, y = latent_umap_2), color = "lightgrey", size = 0.8) +
    geom_point(data = subset(viz_df, cluster == cluster_temp), aes(x = latent_umap_1, y = latent_umap_2), color = colors[i], size = 0.8) +
    theme_void() + theme(legend.position = "none") + labs(title = paste("cluster", cluster_temp))
  
  i <- i + 1
}

png(paste0(output_dir, "latent_latent_umap_per_cluster.png"), width = 1024, height = 768)
do.call(grid.arrange, c(p, ncol = 6))
dev.off()



## Plot individual time points
p <- NULL
i <- 1
colors = getPalette(length(unique(viz_df$timepoint)))

for(cluster_temp in unique(viz_df$timepoint)){
  
  message(cluster_temp)
  
  p[[i]] <- ggplot() +
    geom_point(data = subset(viz_df, patient != cluster_temp), aes(x = latent_umap_1, y = latent_umap_2), color = "lightgrey", size = 0.8) +
    geom_point(data = subset(viz_df, patient == cluster_temp), aes(x = latent_umap_1, y = latent_umap_2), color = colors[i], size = 0.8) +
    theme_void() + theme(legend.position = "none") + labs(title = cluster_temp)
  
  i <- i + 1
}

png(paste0(output_dir, "latent_latent_umap_per_timepoint.png"), width = 1024, height = 333)
do.call(grid.arrange, c(p, ncol = 3))
dev.off()















## Basic clusters stuff



## How much cells with TCRab per cluster?
viz_df %>% group_by(cluster) %>% summarise(has_tcrab = sum(!is.na(new_clonotypes_id)), n = n()) %>% mutate(freq = has_tcrab / n) %>%
  ggplot(aes(cluster, freq, fill = cluster)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(nClusters)) +
  ylim(values = c(0,1)) + theme_bw() + theme(legend.position = "none") + labs(y = "cells with TCRa | TCRb | TCRab") + geom_hline(yintercept = 0.5, linetype = "dotted")
ggsave(paste0(output_dir, "bar_cluster_with_tcrab.png"), width = 8, height = 6)



## How much cells per cluster?
p <- melt(table(Idents(aml_seurat))) %>%
  ggplot(aes(as.factor(reorder(Var1, value)), value, label = value, fill = as.factor(Var1), label = value)) +
  geom_bar(stat = "identity", color = "lightgrey") + geom_text() +
  scale_fill_manual(values = getPalette(nClusters)) +
  theme_minimal() + theme_bw() + theme(legend.position = "none") + coord_flip() + labs(x = "", y = "nCells")
ggsave(plot = p, paste0(output_dir, "bar_cluster.png"), width = 8, height = 8)


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
p <- viz_df %>% group_by(project, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(project), freq, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(nClusters)) + theme_bw() + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_project_per_cluster.png"), width = 6, height = 4)

p <- viz_df %>% group_by(cluster, project) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(as.factor(cluster), freq, fill = as.factor(project))) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_project.png"), width = 12, height = 4)




## Follow up index patient 001
extractTime <- function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*_", "", str1)
}


viz_df$name <- extractName(viz_df$orig.ident)

viz_001_df <- viz_df %>% filter(name == "001") 
viz_001_df$timepoint <- gsub("001_", "", viz_001_df$orig.ident)
viz_001_df$timepoint <- factor(as.character(viz_001_df$timepoint), c("scr", "161017", "201117", "C22D1"))

cluster_df <- select(viz_001_df, cluster, celltype)
cluster_df <- cluster_df[!duplicated(cluster_df), ]

df <- viz_001_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% left_join(cluster_df)

ggplot(df, aes(timepoint, freq, fill = cluster, group = cluster)) + geom_point(shape = 21) + geom_path() + 
  ggrepel::geom_label_repel(data = subset(df, timepoint == "C22D1"), aes(timepoint, freq, label = cluster, fill = cluster), nudge_x = 5,  direction = "y", segment.size = 0.5) +
  theme(legend.position = "none") + facet_wrap(~celltype) + theme_bw() + scale_fill_manual(values = getPalette3(45)) + theme(legend.position = "none")
ggsave(paste0(output_dir, "evolution_per_cluster_001_1.png"), width = 12, height = 8)

ggplot(df, aes(timepoint, freq, fill = cluster, group = cluster)) + geom_point(shape = 21) + geom_path() + 
  ggrepel::geom_label_repel(data = subset(df, timepoint == "C22D1"), aes(timepoint, freq, label = cluster, fill = cluster), nudge_x = 5,  direction = "y", segment.size = 0.5) +
  theme(legend.position = "none") + facet_wrap(~celltype, scales = "free") + theme_bw() + scale_fill_manual(values = getPalette3(45)) + theme(legend.position = "none")
ggsave(paste0(output_dir, "evolution_per_cluster_001_2.png"), width = 12, height = 8)



geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(8)) + theme_bw() + labs(x = "", fill = "") + ggpubr::rotate_x_text(45)
ggsave(plot = p, paste0(output_dir, "bar_per_cluster_per_project.png"), width = 12, height = 4)



