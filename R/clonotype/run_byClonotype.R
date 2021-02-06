


## Init output options
output_dir       <- paste0("results/byClonotype/")
dir.create(output_dir, showWarnings = F)



## No need to modify
#### =====================



########
plotLatentUmapClonotype <- function(seurat_object, clonotype, output_dir, name = ""){
  
  message(clonotype)
  
  if(nchar(name) == 0) name = clonotype
  
  ## Init visualisation object
  viz_df <- cbind(seurat_object@meta.data,
                  seurat_object[["pca"]]@cell.embeddings[,1:6],
                  # seurat_object[["umap"]]@cell.embeddings,
                  seurat_object[["latent_umap"]]@cell.embeddings,
                  cluster = Idents(seurat_object)) %>%
    mutate(cluster = factor(cluster, levels = levels(cluster)[order(levels(cluster))])) %>% dplyr::rename(umap_1 = latent_umap_1, umap_2 = latent_umap_2)
  viz_df$cluster = factor(viz_df$cluster, levels = levels(viz_df$cluster)[order(as.numeric(levels(viz_df$cluster)))])
  nClusters <- viz_df %>% pull(cluster) %>% unique %>% length()
  
  umap_mean <- data.frame(aggregate(umap_1 ~ cluster, viz_df, median),
                          umap_2 = aggregate(umap_2 ~ cluster, viz_df, median)[,2])
  
  ## Visualise
  colors     <- getPalette(nClusters)
  timepoints <- viz_df %>% subset(new_clonotypes_id %in% clonotype) %>% pull(timepoint) %>% unique %>% as.numeric()
  
  ## Plot UMAPs with time points highlighted
  p <- ggplot() +
    geom_point(data = subset(viz_df, !new_clonotypes_id %in% clonotype), aes(x = umap_1, y = umap_2), color = "lightgrey", size = 0.8, alpha = 0.8) +
    geom_smooth(data = subset(viz_df, new_clonotypes_id %in% clonotype), aes(x = umap_1, y = umap_2), color = "darkred", fill = NA) +
    geom_point(data = subset(viz_df, new_clonotypes_id %in% clonotype), aes(x = umap_1, y = umap_2, fill = timepoint), size = 1.5, shape = 21) +
    ggrepel::geom_text_repel(data = umap_mean, aes(x = umap_1, y = umap_2, label = cluster), size = 5, color = "black") +
    
    theme_void() + add_guide +
    scale_fill_manual(values = c("salmon", "dodgerblue", "darkolivegreen4"))
  
  ggsave(plot = p, paste0(output_dir, name, ".png"), width = 8, height = 6)
  
}


## Analyze clonotypes that have at least n cells in total
tcrab_df <- viz_df %>% filter(celltype %in% c("CD4+", "CD8+") & !is.na(new_clonotypes_id))

clonotypes_df <- tcrab_df %>% group_by(new_clonotypes_id) %>% summarise(n = n()) %>% filter(n > 19)
cells.to.use  <- tcrab_df %>% filter(new_clonotypes_id %in% clonotypes_df$new_clonotypes_id) %>% pull(barcode)

t_cell <- subset(cml_seurat, cells = cells.to.use)
t_cell <- getLatentUMAP(t_cell)
nColors <- Idents(t_cell) %>% unique() %>%  length()

p <- DimPlot(t_cell, reduction = "latent_umap", label = T, cols = getPalette(nColors)) + theme_void() + add_guide
ggsave(plot = p, "results/cluster_markers/latent_umap_tcr.png", width = 12, height = 8)

over_n = clonotypes_df %>% pull(new_clonotypes_id)
message(paste("Analyzing", length(over_n), "clonotypes that are over", 20, "cells"))

clonotype_dir <- paste0("results/byClonotype/clonotype_umap/")
dir.create(clonotype_dir, showWarnings = F)

################## Visualize

lapply(over_n, FUN = plotLatentUmapClonotype, seurat_object = t_cell, output_dir = clonotype_dir)


## Specific clonotypes
tcrab_min_df %>% filter(direction_2v1 == "2v1_down") %>% pull(new_clonotypes_id) %>% unique() %>%
  plotLatentUmapClonotype(seurat_object = cd8, output_dir = paste0("results/scrnaseq/", patient, "/byClonotype/", "specific/"), name = "2v1_down_clonotypes")

tcrab_min_df %>% filter(direction_2v1 == "2v1_unsigf") %>% pull(new_clonotypes_id) %>% unique() %>%
  plotLatentUmapClonotype(seurat_object = cd8, output_dir = paste0("results/scrnaseq/", patient, "/byClonotype/", "specific/"), name = "2v1_unsigf")

tcrab_min_df %>% filter(direction_2v1 == "2v1_up") %>% pull(new_clonotypes_id) %>% unique() %>%
  plotLatentUmapClonotype(seurat_object = cd8, output_dir = paste0("results/scrnaseq/", patient, "/byClonotype/", "specific/"), name = "2v1_up")




tcrab_min_df %>% filter(direction_3v2 == "3v2_down") %>% pull(new_clonotypes_id) %>% unique() %>%
  plotLatentUmapClonotype(seurat_object = cd8, output_dir = paste0("results/scrnaseq/", patient, "/byClonotype/", "specific/"), name = "3v2_down_clonotypes")

tcrab_min_df %>% filter(direction_3v2 == "3v2_unsigf") %>% pull(new_clonotypes_id) %>% unique() %>%
  plotLatentUmapClonotype(seurat_object = cd8, output_dir = paste0("results/scrnaseq/", patient, "/byClonotype/", "specific/"), name = "3v2_unsigf")

tcrab_min_df %>% filter(direction_3v2 == "3v2_up") %>% pull(new_clonotypes_id) %>% unique() %>%
  plotLatentUmapClonotype(seurat_object = cd8, output_dir = paste0("results/scrnaseq/", patient, "/byClonotype/", "specific/"), name = "3v2_up")





################## Get DEG




i <- 1
p <- NULL

for(clonotype in unique(tcrab_min_df$new_clonotypes_id)){
  message(i)
  p[[i]] <-  getDEGbyClonotype(seurat_object = cd8, clonotypename = clonotype)
  i <- i + 1
}


DEG_cd8_clonotype <- p %>% rbindlist()
fwrite(DEG_cd8_clonotype, paste0("results/scrnaseq/", patient, "/effectOnClonotypes/cd8/clonotype_deg.txt"), sep = "\t", quote = F, row.names = F)





##############################################################################################################################


## Visualize


## Effect overtime
p <- DEG_cd8_clonotype %>% filter(timepoint %in% c("2v1", "3v2")) %>%
  group_by(new_clonotypes_id,timepoint,direction) %>% summarise(n = n()) %>%
  
  ggplot(aes(timepoint, n, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
  scale_fill_manual(values = getPalette(8)) + theme_bw() + theme(legend.position = "none") + labs(x = "", y = "# DEG after therapy initiation") + facet_wrap(~direction) + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("2v1", "3v2")))
ggsave(plot = p, paste0("results/scrnaseq/", patient, "/effectOnClonotypes/cd8/box_nDEG_timepoints.png"), width = 6, height = 4)




p <- DEG_cd8_clonotype %>% filter(timepoint %in% c("2v1", "3v2")) %>%
  group_by(new_clonotypes_id,timepoint) %>% summarise(n = n()) %>%
  
  ggplot(aes(timepoint, n, fill = timepoint, label=substr(new_clonotypes_id, 10, nchar(new_clonotypes_id)))) + geom_boxplot(outlier.shape = NA, alpha=0.1) + geom_path(aes(group = new_clonotypes_id)) + ggrepel::geom_text_repel() +
  scale_fill_manual(values = getPalette(8)) + theme_bw() + theme(legend.position = "none") + labs(x = "", y = "# DEG after therapy initiation") +
  ggsignif::geom_signif(comparisons = list(c("2v1", "3v2")))
ggsave(plot = p, paste0("results/scrnaseq/", patient, "/effectOnClonotypes/cd8/path_nDEG_timepoints.png"), width = 6, height = 4)




## Expanded vs suppressed clonotypes over time
df <- DEG_cd8_clonotype %>% filter(timepoint == "2v1") %>%
  group_by(new_clonotypes_id,timepoint,direction) %>% summarise(n = n()) %>%
  left_join(select(tcrab_min_df,new_clonotypes_id, direction_2v1), by = "new_clonotypes_id")
df <- df[!duplicated(df), ]

df %>%
  ggplot(aes(direction_2v1, n, fill = direction_2v1)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
  scale_fill_manual(values = getPalette3(8)) + theme(legend.position = "none") + theme_bw() + facet_wrap(~direction) + facets_nice +
  labs(x = "", y = "# DEG after therapy initiation")
ggsave(plot = p, paste0("results/scrnaseq/", patient, "/effectOnClonotypes/cd8/box_2v1_expansion.png"), width = 4, height = 3)




df <- DEG_cd8_clonotype %>% filter(timepoint == "3v2") %>%
  group_by(new_clonotypes_id,timepoint,direction) %>% summarise(n = n()) %>%
  left_join(select(tcrab_min_df,new_clonotypes_id, direction_3v2), by = "new_clonotypes_id")
df <- df[!duplicated(df), ]

df %>%
  ggplot(aes(direction_3v2, n, fill = direction_3v2)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
  scale_fill_manual(values = getPalette3(8)) + theme(legend.position = "none") + theme_bw() + facet_wrap(~direction) + facets_nice +
  labs(x = "", y = "# DEG after therapy initiation")
ggsave(plot = p, paste0("results/scrnaseq/", patient, "/effectOnClonotypes/cd8/box_3v2_expansion.png"), width = 4, height = 3)



