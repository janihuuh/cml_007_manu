

dir.create("results/tcr/", showWarnings = F)

viz_df <- cbind(cml_seurat@meta.data,
                # cml_seurat[["pca"]]@cell.embeddings[,1:6],
                cml_seurat[["latent_umap"]]@cell.embeddings,
                cluster = Idents(cml_seurat)) # %>% dplyr::rename(latent_umap_1 = latent_1, latent_umap_2 = latent_2)
viz_df$celltype <- extractCoarsePhenotype(viz_df$cluster)
tcrab_df <- viz_df %>% filter(celltype %in% c("CD4+", "CD8+") & !is.na(new_clonotypes_id))

clonality <- function(x){
  
  1 - ineq::entropy(x) / log(length(unique((x))))
  # 1 - ineq::entropy(x) / log(length(((x))))
  
}


#############################################

tcrab_df %>% group_by(cluster) %>% 
  summarise(clonality = 1 - ineq::Gini(as.factor(new_clonotypes_id))) %>% 
  # summarise(clonality = clonality(as.factor(new_clonotypes_id))) %>% 
  ggplot(aes(reorder(cluster, clonality), clonality, fill = cluster)) + geom_bar(stat = "identity") + labs(x = "") + theme_bw() + 
  geom_hline(yintercept = 0.5, linetype = "dotted") + ylim(values = c(0,1)) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + scale_fill_manual(values = getPalette(7))
ggsave("results/tcr/bar_clonality_cluster.png", width = 4, height = 4)

tcrab_df %>% group_by(cluster, patient) %>% 
  summarise(clonality = 1 - ineq::Gini(as.factor(new_clonotypes_id))) %>% 
  # summarise(clonality = clonality(as.factor(new_clonotypes_id))) %>% 
  ggplot(aes(reorder(cluster, clonality), clonality, fill = cluster)) + geom_boxplot(outlier.shape = NA) + labs(x = "") + theme_bw() + 
  geom_hline(yintercept = 0.5, linetype = "dotted") + ylim(values = c(0,1)) + scale_fill_manual(values = getPalette(7)) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none")
ggsave("results/tcr/box_clonality_cluster.png", width = 4, height = 4)

tcrab_df %>% group_by(patient, timepoint, effusion) %>% 
  summarise(clonality = 1 - ineq::Gini(as.factor(new_clonotypes_id))) %>% 
  # summarise(clonality = clonality(as.factor(new_clonotypes_id))) %>% 
  ggplot(aes(timepoint, clonality, group = patient, fill = effusion)) + geom_point(shape = 21, size = 3) + geom_path() + ylim(values = c(0.5,1)) + theme_bw() + scale_fill_manual(values = c("salmon", "lightgrey"))
ggsave("results/tcr/box_clonality_effusion.png", width = 4, height = 4)

tcrab_df %>% group_by(patient, timepoint, effusion) %>% 
  summarise(clonality = 1 - ineq::Gini(as.factor(new_clonotypes_id))) %>% 
  filter(timepoint == "dg") %>% 

  ggplot(aes(as.factor(effusion), clonality, fill = effusion)) + geom_boxplot(outlier.shape = NA) + ylim(values = c(0,1)) + 
    ggsignif::geom_signif(comparisons = list(c("effusion", "no effusion"))) + theme_bw() + labs(x = "") + theme(legend.position = "none") + geom_jitter(size = 0.3) + scale_fill_manual(values = c("salmon", "lightgrey"))
ggsave("results/tcr/box_clonality_effusion_baseline.png", width = 4, height = 4)

#############################################


tcrab_df %>% group_by(new_clonotypes_id) %>% summarise(n = n()) %>% arrange(desc(n)) %>% mutate(patient = extractName(new_clonotypes_id)) %>% 
  group_by(patient) %>% top_n(n = 10, wt = n) %>% 
  
  ggplot(aes(patient,n,fill=new_clonotypes_id,label=n)) + geom_bar(stat = "identity", position = "dodge", color = "black") + theme_bw() + theme(legend.position = "none") + geom_text(position = position_dodge(width = 1)) +
  labs(y = "nCells in clonotype")
ggsave("results/tcr/bar_top10_clonesize.png", width = 8, height = 4)



clonesize_df <- tcrab_df %>% group_by(new_clonotypes_id) %>% summarise(n = n()) %>% dplyr::rename(total_n = n)
top10_clonesize_df <- clonesize_df %>% mutate(patient = extractName(new_clonotypes_id)) %>% group_by(patient) %>% top_n(n = 10, wt = total_n)

top10_clonesize_df %>% 
  ggplot(aes(patient,total_n,fill=new_clonotypes_id,label=total_n)) + geom_bar(stat = "identity", position = "dodge", color = "black") + theme_bw() + theme(legend.position = "none") + geom_text(position = position_dodge(width = 1)) +
  labs(y = "nCells in clonotype")
ggsave("results/tcr/bar_top10_clonesize.png", width = 8, height = 4)



tcrab_df %>% filter(new_clonotypes_id %in% top10_clonesize_df$new_clonotypes_id) %>% 
  group_by(new_clonotypes_id, timepoint) %>% summarise(n = n()) %>% arrange(desc(n)) %>% left_join(top10_clonesize_df, by = "new_clonotypes_id") %>% 
  
  ggplot(aes(timepoint,n,group=new_clonotypes_id,fill=patient)) + geom_point(shape = 21, size = 3) + geom_path() + facet_wrap(~patient, scales = "free") +
  theme_bw() + theme(legend.position = "none") + 
  labs(y = "nCells in clonotype") + facets_nice
ggsave("results/tcr/bar_top10_clonesize_timepoint.png", width = 8, height = 6)



