
## 
ifna_hallmark <- fread("results/interferon_alpha/hallmark_ifna.txt")[-1,]
cml_seurat <- AddModuleScore(cml_seurat, list(ifna_hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE), name = "hallmark_ifna")


ifna_hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE

cml_seurat$cluster <- Idents(cml_seurat)
 
cml_seurat@meta.data %>% 
  group_by(cluster, hallmark_ifna1) %>% 
  ggplot(aes(cluster, hallmark_ifna1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters))  + ggpubr::rotate_x_text(45)
ggsave("results/interferon_alpha/violin_timepoint_cluster_total.png", width = 8, height = 4)

cml_seurat@meta.data %>% filter(timepoint == "dg") %>% 
  group_by(cluster, hallmark_ifna1) %>% 
  ggplot(aes(cluster, hallmark_ifna1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) + ggpubr::rotate_x_text(45)
ggsave("results/interferon_alpha/violin_timepoint_cluster_dg.png", width = 8, height = 4)

cml_seurat@meta.data %>% filter(timepoint == "12mo") %>% 
  group_by(cluster, hallmark_ifna1) %>% 
  ggplot(aes(cluster, hallmark_ifna1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) + ggpubr::rotate_x_text(45)
ggsave("results/interferon_alpha/violin_timepoint_cluster_12mo.png", width = 8, height = 4)

cml_seurat@meta.data %>% 
  group_by(timepoint, cluster, hallmark_ifna1) %>% 
  ggplot(aes(timepoint, hallmark_ifna1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) +
  facets_nice + 
  facet_wrap(~cluster, scales = "free_y") +
  ggsignif::geom_signif(comparisons = list(c("dg", "12mo"))) 
ggsave("results/interferon_alpha/violin_timepoint_cluster.png", width = 12, height = 8)

cml_seurat@meta.data %>% filter(timepoint == "dg") %>% 
  group_by(cluster, hallmark_ifna1, effusion) %>% 
  ggplot(aes(effusion, hallmark_ifna1, fill = effusion)) + geom_violin(alpha = 0.5) + geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme(legend.position = "none") + scale_fill_manual(values = c("salmon", "lightgrey")) + ggpubr::rotate_x_text(45) +
  facet_wrap(~cluster, scales = "free") + facets_nice + theme_bw() + theme(legend.position = "none") +
  ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion"))) 
ggsave("results/interferon_alpha/violin_effusion_dg.png", width = 12, height = 12)

cml_seurat@meta.data %>% filter(timepoint == "12mo") %>% 
  group_by(cluster, hallmark_ifna1, effusion) %>% 
  ggplot(aes(effusion, hallmark_ifna1, fill = effusion)) + geom_violin(alpha = 0.5) + geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme(legend.position = "none") + scale_fill_manual(values = c("salmon", "lightgrey")) + ggpubr::rotate_x_text(45) +
  facet_wrap(~cluster, scales = "free") + facets_nice + theme_bw() + theme(legend.position = "none") +
  ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion"))) 
ggsave("results/interferon_alpha/violin_effusion_12mo.png", width = 12, height = 12)

cml_seurat@meta.data %>% filter(effusion = "no effusion") %>% 
  group_by(timepoint, cluster, hallmark_ifna1) %>% 
  ggplot(aes(timepoint, hallmark_ifna1, fill = cluster)) + geom_violin() + theme(legend.position = "none") + scale_fill_manual(values = getPalette(nClusters)) +
  facets_nice + 
  facet_wrap(~cluster, scales = "free_y") +
  ggsignif::geom_signif(comparisons = list(c("dg", "12mo"))) 
ggsave("results/interferon_alpha/violin_timepoint_cluster.png", width = 12, height = 8)





######





## IFNa at 3mo
df <- cml_seurat@meta.data %>% filter(timepoint == "3mo") %>% group_by(cluster,effusion) %>% summarise(median=median(hallmark_ifna1))
df_eff <- df %>% filter(effusion == "effusion") %>% dplyr::select(-effusion)
df_noeff <- df %>% filter(effusion == "no effusion") %>% dplyr::select(-effusion)

ord <- df_eff %>% left_join(df_noeff, by = "cluster") %>% mutate(log2fc = log2(median.x/median.y)) %>% arrange(desc(log2fc))
ord %>% ggplot(aes(reorder(cluster,log2fc),log2fc, fill = ifelse(abs(log2fc)>1, "y", "n"))) + geom_bar(stat="identity") + coord_flip() + ylim(values=c(-1.5,1.5)) + theme(legend.position = "none")

cml_seurat@meta.data %>% 
  filter(timepoint == "3mo") %>% left_join(ord, by ="cluster") %>% 
  ggplot(aes(effusion,hallmark_ifna1)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(~reorder(cluster,log2fc), scales = "free_y")



## IFNa at 12mo
cml_seurat$late_effusion <- ifelse(cml_seurat$patient == "720", "late effusion", "no late effusion")
cml_seurat$late_effusion <- ifelse(cml_seurat$patient == "716", "early effusion", cml_seurat$late_effusion)

df <- cml_seurat@meta.data %>% filter(timepoint == "12mo") %>% group_by(cluster,effusion) %>% summarise(median=median(hallmark_ifna1))
df_eff <- df %>% filter(effusion == "effusion") %>% dplyr::select(-effusion)
df_noeff <- df %>% filter(effusion == "no effusion") %>% dplyr::select(-effusion)

ord <- df_eff %>% left_join(df_noeff, by = "cluster") %>% mutate(log2fc = log2(median.x/median.y)) %>% arrange(desc(log2fc))
ord %>% ggplot(aes(reorder(cluster,log2fc),log2fc, fill = ifelse(abs(log2fc)>1, "y", "n"))) + geom_bar(stat="identity") + coord_flip() + ylim(values=c(-1.5,1.5)) + theme(legend.position = "none")

cml_seurat@meta.data %>% 
  filter(late_effusion != "early effusion") %>% 
  filter(timepoint == "12mo") %>% left_join(ord, by ="cluster") %>% 
  filter(celltype %in% c("CD8", "NK")) %>% 
  mutate(patient = factor(as.character(patient), levels = c("706", "730", "720"))) %>% 
  
  ggplot(aes(patient,hallmark_ifna1,fill=late_effusion)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(~reorder(cluster,log2fc), scales = "free_y")



cml_seurat@meta.data %>% 
  filter(late_effusion != "early effusion") %>% 
  filter(timepoint == "3mo") %>% left_join(ord, by ="cluster") %>% 
  filter(celltype %in% c("CD8", "NK")) %>% 
  mutate(patient = factor(as.character(patient), levels = c("706", "730", "720"))) %>% 
  
  ggplot(aes(patient,hallmark_ifna1,fill=late_effusion)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(~reorder(cluster,log2fc), scales = "free_y")





######

df_eff     <- cml_seurat@meta.data %>% filter(effusion == "effusion") %>% group_by(cluster,timepoint) %>% summarise(median=median(hallmark_ifna1))
df_eff_0m  <- df_eff %>% filter(timepoint == "dg") %>% dplyr::select(-timepoint)
df_eff_3m  <- df_eff %>% filter(timepoint == "3mo") %>% dplyr::select(-timepoint)
df_eff_12m <- df_eff %>% filter(timepoint == "12mo") %>% dplyr::select(-timepoint)

df_no_eff     <- cml_seurat@meta.data %>% filter(effusion == "no effusion") %>% group_by(cluster,timepoint) %>% summarise(median=median(hallmark_ifna1))
df_no_eff_0m  <- df_no_eff %>% filter(timepoint == "dg") %>% dplyr::select(-timepoint)
df_no_eff_3m  <- df_no_eff %>% filter(timepoint == "3mo") %>% dplyr::select(-timepoint)
df_no_eff_12m <- df_no_eff %>% filter(timepoint == "12mo") %>% dplyr::select(-timepoint)

df_eff_12v3   <- df_eff_3m %>% left_join(df_eff_12m, by = "cluster") %>% mutate(log2fc = log2(median.y/median.x)) %>% arrange(desc(log2fc))
df_noeff_12v3 <- df_no_eff_3m %>% left_join(df_no_eff_12m, by = "cluster") %>% mutate(log2fc = log2(median.y/median.x)) %>% arrange(desc(log2fc))

df_eff_12v3 %>% left_join(df_noeff_12v3, by = "cluster") %>% 
  
  ggplot(aes(log2fc.x,log2fc.y,label=cluster)) + geom_point() + xlim(c(-1,1)) + ylim(c(-1,1)) + theme_bw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_abline(slope=1) +
  labs(x = "log2fc in effusion", y = "log2fc in no effusion") + ggrepel::geom_text_repel()



cml_seurat@meta.data %>% 
  filter(celltype %in% c("CD8", "NK")) %>% 
  filter(timepoint %in% c("3mo", "12mo")) %>% 
  # filter(!cluster %in% c("12 CD8 MAIT-like CXCR3+", "13 NK CD56bright", "1 CD8 EM GZMK+ CXCR3+", "2 CD8 cytotoxic GZMH+", "8 CD8 naive TCF7+", "5 CD8 EMRA ZNF683+")) %>% 
  filter(!cluster %in% c("12 CD8 MAIT-like CXCR3+", "13 NK CD56bright")) %>% 
  
  ggplot(aes(timepoint,hallmark_ifna1,fill=effusion)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(cluster~effusion, ncol=2)


cml_seurat@meta.data %>% 
  filter(celltype %in% c("CD8", "NK")) %>% 
  # filter(timepoint %in% c("3mo", "12mo")) %>% 
  # filter(!cluster %in% c("12 CD8 MAIT-like CXCR3+", "13 NK CD56bright", "1 CD8 EM GZMK+ CXCR3+", "2 CD8 cytotoxic GZMH+", "8 CD8 naive TCF7+", "5 CD8 EMRA ZNF683+")) %>% 
  filter(!cluster %in% c("12 CD8 MAIT-like CXCR3+", "13 NK CD56bright")) %>% 
  
  ggplot(aes(timepoint,hallmark_ifna1,fill=effusion)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(~cluster) + ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion")))



