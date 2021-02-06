

## No need to modify
#### =====================

## Init output options
output_dir2 <- "results/clusters/"
dir.create(output_dir2, showWarnings = F)

output_dir <- "results/clusters/scores"
dir.create(output_dir, showWarnings = F)

## Calculate Van Galen
i <- 1
for(markers in van_galen_markers){
  
  name <- names(van_galen_markers)[i]
  message(name)
  cml_seurat <- AddModuleScore(cml_seurat, features = list(markers), name = name)
  i <- i + 1
  
}


## Plot Van Galen
cml_seurat <- fixSeurat(cml_seurat)

i <- 1
for(markers in van_galen_markers){
  
  name <- names(van_galen_markers)[i]
  message(name)
  p <- FeaturePlot(cml_seurat, reduction = "latent_umap", feature = paste0(name, "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T)+ theme_void()
  ggsave(plot = p, paste0("results/clusters/scores/vangalen_", name, ".png"), width = 6, height = 4)
  i <- i + 1
  
}



## Violin plots
scores <- lapply(van_galen_markers, function(x){
  
  message(x)
  score <- Matrix::colMeans(cml_seurat@assays$RNA[rownames(cml_seurat@assays$RNA) %in% x, ])
  return(score)
}
)

scores_df <- do.call(cbind, scores)
colnames(scores_df) <- make.unique(colnames(scores_df))
scores_df <- scores_df %>% as.data.frame() %>% mutate(cluster = Idents(cml_seurat), barcode = colnames(cml_seurat))
scores_molten <- melt(scores_df, id = c("barcode", "cluster"))

p <- scores_molten %>% # filter(variable == "pDC") %>%
  ggplot(aes(cluster,value,fill=cluster)) + geom_violin(alpha = 0.4) + geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = getPalette(nClusters)) +
  scale_fill_manual(values = getPalette(nClusters)) + facet_wrap(~variable,ncol=2, scales = "free") + scale_y_log10() # + facets_nice # + scale_y_log10()
ggsave(plot = p, paste0(output_dir, "violin_van_galen_scores.png"), width = 24, height = 30)






##################################################




## Calculate zheng_markers
i <- 1
for(markers in zheng_markers){
  
  name <- names(zheng_markers)[i]
  message(name)
  cml_seurat <- AddModuleScore(cml_seurat, features = list(markers), name = name)
  i <- i + 1
  
  
}



## Plot zheng_markers
i <- 1
for(markers in zheng_markers){
  
  name <- names(zheng_markers)[i]
  message(name)
  p <- FeaturePlot(cml_seurat, reduction = "latent_umap", feature = paste0(name, "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T)+ theme_void()
  ggsave(plot = p, paste0("results/clusters/scores/zheng_", name, ".png"), width = 6, height = 4)
  i <- i + 1
  
}






## Calculate zhang_cd4_markers
i <- 1
for(markers in zhang_cd4_markers){
  
  name <- names(zhang_cd4_markers)[i]
  message(name)
  cml_seurat <- AddModuleScore(cml_seurat, features = list(markers), name = name)
  i <- i + 1
  
  
}


## Plot zhang_cd4_markers
i <- 1
for(markers in zhang_cd4_markers){
  
  name <- names(zhang_cd4_markers)[i]
  message(name)
  p <- FeaturePlot(tnk_cells, reduction = "latent_umap", feature = paste0(name, "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()
  ggsave(plot = p, paste0("results/clusters/scores/zhang_cd4_", name, ".png"), width = 6, height = 4)
  i <- i + 1
  
}




#############


## Calculate zhang_cd8_markers
i <- 1
for(markers in zhang_cd8_markers){
  
  name <- names(zhang_cd8_markers)[i]
  message(name)
  cml_seurat <- AddModuleScore(cml_seurat, features = list(markers), name = name)
  i <- i + 1
  
  
}


## Plot zhang_cd8_markers
i <- 1
for(markers in zhang_cd8_markers){
  
  name <- names(zhang_cd8_markers)[i]
  message(name)
  p <- FeaturePlot(tnk_cells, reduction = "latent_umap", feature = paste0(name, "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()
  ggsave(plot = p, paste0("results/clusters/scores/zhang_cd8_", name, ".png"), width = 6, height = 4)
  i <- i + 1
  
}


## Plot gammadelta genes
gd_genes <- c(grep("^TRD", rownames(fhrb1641), value = T), grep("^TRG", rownames(fhrb1641), value = T))
fhrb1641 <- AddModuleScore(tcell, features = list(gd_genes), name = "gammadelta")
p <- FeaturePlot(fhrb1641, reduction = "latent_umap", feature = paste0("gammadelta", "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()
ggsave(plot = p, paste0("results/clusters/scores/zhang_cd8_", "gammadelta", ".png"), width = 6, height = 4)



## Plot CD8 genes
cd8_genes <- c("CD8A", "CD8B")
cd8 <- AddModuleScore(tcell, features = list(cd8_genes), name = "cd8")
p <- FeaturePlot(cd8, reduction = "latent_umap", feature = paste0("cd8", "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()
ggsave(plot = p, paste0("results/clusters/scores/zhang_cd8_", "cd8", ".png"), width = 6, height = 4)




## Plot cytot_genes genes
p <- FeaturePlot(tcell, reduction = "latent_umap", feature = paste0("cytotoxicity", "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()
ggsave(plot = p, paste0("results/clusters/scores/zhang_cd8_", "cytotoxicity", ".png"), width = 6, height = 4)


p <- FeaturePlot(tcell, reduction = "latent_umap", feature = "percent.cycle", cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()



###### NK-cell markers

nk_cells <- fixSeurat(nk_cells)

## Calculate nk_pfefferle
i <- 1
for(markers in nk_mark){
  
  name <- names(nk_mark)[i]
  message(name)
  fhrb1641 <- AddModuleScore(fhrb1641, features = list(markers), name = name) #, nbin = 10)
  i <- i + 1
  
}


## Plot nk_pfefferle
i <- 1
for(markers in nk_mark){
  
  name <- names(nk_mark)[i]
  message(name)
  p <- FeaturePlot(fhrb1641, reduction = "latent_umap", feature = paste0(name, "1"), cols = c("gray91", "tomato"), order = T, min.cutoff = 0, label = T, repel = T) + theme_void()
  ggsave(plot = p, paste0("results/clusters/scores/nk_pfefferle", name, ".png"), width = 6, height = 4)
  i <- i + 1
  
}













###### Violin plots







## Zhang scores
zhang_cd4_scores <- lapply(zhang_cd4_markers, function(x){
  
  message(x)
  score <- Matrix::colMeans(tcell@assays$RNA[rownames(tcell@assays$RNA) %in% x, ])
  return(score)
}
)

zhang_cd4_scores_df <- do.call(cbind, zhang_cd4_scores)
colnames(zhang_cd4_scores_df) <- make.unique(colnames(zhang_cd4_scores_df))
zhang_cd4_scores_df <- zhang_cd4_scores_df %>% as.data.frame() %>% mutate(cluster = Idents(tcell), barcode = colnames(tcell))
zhang_cd4_scores_molten <- melt(zhang_cd4_scores_df, id = c("barcode", "cluster"))


p <- zhang_cd4_scores_molten %>%
  ggplot(aes(cluster,value,fill=cluster)) + geom_violin(alpha = 0.4) + geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = getPalette(nClusters)) +
  scale_fill_manual(values = getPalette(nClusters)) + facet_wrap(~variable,ncol=2, scales = "free") + facets_nice # + scale_y_log10()
ggsave(plot = p, paste0(output_dir, "violin_zhang_cd4_scores.png"), width = 15, height = 30)




zhang_cd8_scores <- lapply(zhang_cd8_markers, function(x){
  
  message(x)
  score <- psych::geometric.mean(tcell@assays$RNA[rownames(tcell@assays$RNA) %in% x, ])
  return(score)
  
}
)

zhang_cd8_scores_df <- do.call(cbind, zhang_cd8_scores)
colnames(zhang_cd8_scores_df) <- make.unique(colnames(zhang_cd8_scores_df))
zhang_cd8_scores_df <- zhang_cd8_scores_df %>% as.data.frame() %>% mutate(cluster = Idents(tcell), barcode = colnames(tcell))
zhang_cd8_scores_molten <- melt(zhang_cd8_scores_df, id = c("barcode", "cluster"))


p <- zhang_cd8_scores_molten %>% # filter(variable == "pDC") %>%
  ggplot(aes(cluster,value,fill=cluster)) + geom_violin(alpha = 0.4) + geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = getPalette(nClusters)) +
  scale_fill_manual(values = getPalette(nClusters)) + facet_wrap(~variable,ncol=2, scales = "free") + facets_nice # + scale_y_log10()
ggsave(plot = p, paste0(output_dir, "violin_zhang_cd8_scores.png"), width = 15, height = 30)




