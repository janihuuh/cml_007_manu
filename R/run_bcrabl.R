
## 
bcrabl_hallmark <- fread("results/bcr_abl/klein_targets_of_bcr_abl_fusion.txt")[-1,]
ray_bcrabl_hallmark <- fread("results/bcr_abl/ray_targets_of_bcr_abl1_fusion.txt")[-1,]

cml_seurat <- AddModuleScore(cml_seurat, list(bcrabl_hallmark$KLEIN_TARGETS_OF_BCR_ABL1_FUSION), name = "klein_bcr_abl")
cml_seurat <- AddModuleScore(cml_seurat, list(ray_bcrabl_hallmark$RAY_TARGETS_OF_P210_BCR_ABL_FUSION_UP), name = "ray_bcr_abl")

cml_seurat$cluster <- Idents(cml_seurat)

cml_seurat@meta.data %>% 
  group_by(cluster, klein_bcr_abl1) %>% 
  ggplot(aes(cluster, klein_bcr_abl1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters))  + ggpubr::rotate_x_text(45)
ggsave("results/bcr_abl/violin_timepoint_cluster_total.png", width = 8, height = 4)


cml_seurat@meta.data %>% filter(timepoint == "dg") %>% 
  group_by(cluster, klein_bcr_abl1) %>% 
  ggplot(aes(cluster, klein_bcr_abl1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) + ggpubr::rotate_x_text(45)
ggsave("results/bcr_abl/violin_timepoint_cluster_dg.png", width = 8, height = 4)

cml_seurat@meta.data %>% filter(timepoint == "3mo") %>% 
  group_by(cluster, klein_bcr_abl1) %>% 
  ggplot(aes(cluster, klein_bcr_abl1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) + ggpubr::rotate_x_text(45)
ggsave("results/bcr_abl/violin_timepoint_cluster_3mo.png", width = 8, height = 4)

cml_seurat@meta.data %>% filter(timepoint == "12mo") %>% 
  group_by(cluster, klein_bcr_abl1) %>% 
  ggplot(aes(cluster, klein_bcr_abl1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) + ggpubr::rotate_x_text(45)
ggsave("results/bcr_abl/violin_timepoint_cluster_12mo.png", width = 8, height = 4)




cml_seurat@meta.data %>% 
  group_by(timepoint, cluster, klein_bcr_abl1) %>% 
  ggplot(aes(timepoint, klein_bcr_abl1, fill = cluster)) + 
  geom_violin(alpha = 0.6) + geom_boxplot(outlier.shape = NA, width = 0.4) +
  theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) +
  facets_nice +  
  facet_wrap(~cluster, scales = "free_y") + 
  ggsignif::geom_signif(comparisons = list(c("dg", "3mo")), manual = F) + 
  # ggpubr::stat_compare_means(comparisons = list(c("dg", "3mo"))) +
  geom_hline(yintercept = 0, linetype = "dotted") 
ggsave("results/bcr_abl/violin_klein_timepoint_cluster.png", width = 12, height = 8)

cml_seurat@meta.data %>% 
  group_by(timepoint, cluster, ray_bcr_abl1) %>% 
  ggplot(aes(timepoint, ray_bcr_abl1, fill = cluster)) + 
    geom_violin(alpha = 0.6) + geom_boxplot(outlier.shape = NA, width = 0.4) +
    theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) +
    facets_nice +  
    facet_wrap(~cluster, scales = "free_y") + 
    ggsignif::geom_signif(comparisons = list(c("dg", "3mo")), manual = F) + 
    # ggpubr::stat_compare_means(comparisons = list(c("dg", "3mo"))) +
    geom_hline(yintercept = 0, linetype = "dotted") 
ggsave("results/bcr_abl/violin_ray_timepoint_cluster.png", width = 12, height = 8)





cml_seurat@meta.data %>% filter(timepoint == "dg") %>% 
  group_by(cluster, klein_bcr_abl1, effusion) %>% 
  ggplot(aes(effusion, klein_bcr_abl1, fill = effusion)) + geom_violin(alpha = 0.5) + geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme(legend.position = "none") + scale_fill_manual(values = c("salmon", "lightgrey")) + ggpubr::rotate_x_text(45) +
  facet_wrap(~cluster, scales = "free") + facets_nice + theme_bw() + theme(legend.position = "none") +
  ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion"))) 
ggsave("results/bcr_abl/violin_effusion_dg.png", width = 12, height = 12)

cml_seurat@meta.data %>% filter(timepoint == "3mo") %>% 
  group_by(cluster, klein_bcr_abl1, effusion) %>% 
  ggplot(aes(effusion, klein_bcr_abl1, fill = effusion)) + geom_violin(alpha = 0.5) + geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme(legend.position = "none") + scale_fill_manual(values = c("salmon", "lightgrey")) + ggpubr::rotate_x_text(45) +
  facet_wrap(~cluster, scales = "free") + facets_nice + theme_bw() + theme(legend.position = "none") +
  ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion"))) 
ggsave("results/bcr_abl/violin_effusion_3mo.png", width = 12, height = 12)

cml_seurat@meta.data %>% filter(timepoint == "12mo") %>% 
  group_by(cluster, klein_bcr_abl1, effusion) %>% 
  ggplot(aes(effusion, klein_bcr_abl1, fill = effusion)) + geom_violin(alpha = 0.5) + geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme(legend.position = "none") + scale_fill_manual(values = c("salmon", "lightgrey")) + ggpubr::rotate_x_text(45) +
  facet_wrap(~cluster, scales = "free") + facets_nice + theme_bw() + theme(legend.position = "none") +
  ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion"))) 
ggsave("results/bcr_abl/violin_effusion_12mo.png", width = 12, height = 12)



cml_seurat@meta.data %>% filter(effusion = "no effusion") %>% 
  group_by(timepoint, cluster, klein_bcr_abl1) %>% 
  ggplot(aes(timepoint, klein_bcr_abl1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) +
  facets_nice + 
  facet_wrap(~cluster, scales = "free_y") +
  ggsignif::geom_signif(comparisons = list(c("dg", "12mo"))) 
ggsave("results/bcr_abl/violin_timepoint_cluster_12mo_r.png", width = 12, height = 8)
