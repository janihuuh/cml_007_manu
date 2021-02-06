
dir.create("results/manuscript/pe/")

## 
facs_24mo <- fread("results/manuscript/facs_df.txt")

pe_patients <- c("1444", "1463", "1421", "1545", "1582", "1478")
pe_patients <- c("1444", "1421", "1545", "1582", "1478")

pe_df <- facs_24mo %>% filter(FM %in% pe_patients)
pe_df <- pe_df %>% melt(id = c("timepoint","Study.nro","FM", "HRUH", "pe")) %>% dplyr::select(-Study.nro, -HRUH) %>% dplyr::rename(name = FM) %>% dplyr::select(name,timepoint,variable,value) %>% mutate(type = "facs")
pe_df <- pe_df %>% filter(!is.na(name)) %>% filter(timepoint != "24mo")

## Krskl
p.df <- lapply(unique(pe_df$variable), FUN = function(x){
  message(x)
  y <- pe_df %>% filter(variable == x)
  if(length(unique(y$timepoint)) > 2){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)
fwrite(p.df, "results/manuscript/pe/krskl_pe_only_p_df.txt", sep = "\t", quote = F, row.names = F)


## addition of IFNa
pe_df_ifna <- pe_df %>% filter(timepoint %in% c("3mo", "12mo"))

p.df <- lapply(unique(pe_df_ifna$variable), FUN = function(x){
  message(x)
  y <- pe_df_ifna %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 2){
    median.x = y %>% filter(timepoint == "3mo") %>% pull(value) %>% median(na.rm = T)
    median.y = y %>% filter(timepoint == "12mo") %>% pull(value) %>% median(na.rm = T)
    wilcox.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x, median.3m = median.x, median.12m = median.y, dir = ifelse(log2(median.y/median.x) < 0, "up", "down"))
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) %>% left_join(var_df)
fwrite(p.df, "results/manuscript/pe/ifna_p_df.txt", sep = "\t", quote = F, row.names = F)


## addition of ifna, no PE
no_pe_df <- facs_24mo %>% filter(!FM %in% pe_patients)
no_pe_df <- no_pe_df %>% melt(id = c("timepoint","Study.nro","FM", "HRUH", "pe")) %>% dplyr::select(-Study.nro, -HRUH) %>% dplyr::rename(name = FM) %>% dplyr::select(name,timepoint,variable,value) %>% mutate(type = "facs")
no_pe_df <- no_pe_df %>% filter(!is.na(name)) %>% filter(timepoint != "24mo")
no_pe_df_ifna <- no_pe_df %>% filter(timepoint %in% c("3mo", "12mo"))

p.df <- lapply(unique(no_pe_df_ifna$variable), FUN = function(x){
  message(x)
  y <- no_pe_df_ifna %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 2){
    median.x = y %>% filter(timepoint == "3mo") %>% pull(value) %>% median(na.rm = T)
    median.y = y %>% filter(timepoint == "12mo") %>% pull(value) %>% median(na.rm = T)
    wilcox.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x, median.3m = median.x, median.12m = median.y, dir = ifelse(log2(median.y/median.x) < 0, "up", "down"))
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) %>% left_join(var_df)
fwrite(p.df, "results/manuscript/pe/ifna_p_df.txt", sep = "\t", quote = F, row.names = F)



matrix(c(4, 1, 14, 21),ncol=2,byrow = T) %>% fisher.test(alternative = "greater")
matrix(c(5, 1, 13, 21),ncol=2,byrow = T) %>% fisher.test(alternative = "greater")
matrix(c(6, 2, 12, 20),ncol=2,byrow = T) %>% fisher.test(alternative = "two.sided")
matrix(c(2, 4, 16, 18),ncol=2,byrow = T) %>% fisher.test(alternative = "two.sided")




facs_24mo$pe <- ifelse(facs_24mo$FM %in% pe_patients, "PE/PAH", "no PE/PAH")

facs_24mo %>% filter(timepoint=="0mo") %>% 
  ggplot(aes(pe,CD8.TEMRA)) + geom_boxplot() + ggpubr::stat_compare_means() + geom_jitter() 

facs_24mo %>% filter(timepoint != "24mo") %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("0mo", "3mo", "12mo", "24mo"))) %>% 
  ggplot(aes(timepoint,CD8.TEMRA)) + geom_boxplot() + ggpubr::stat_compare_means() + geom_jitter() + facet_wrap(~pe) + ggsignif::geom_signif(comparisons = list(c("3mo", "12mo")))

facs_24mo %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("0mo", "3mo", "12mo", "24mo"))) %>% 
  ggplot(aes(timepoint,NK)) + geom_boxplot() + ggpubr::stat_compare_means() + geom_jitter() + facet_wrap(~pe)




## 
cells.to.keep   <- cml_seurat@meta.data %>% filter(effusion == "effusion") %>% pull(barcode)
effusion_seurat <- subset(cml_seurat, cells = cells.to.keep) %>% getLatentUMAP()

cells.to.keep   <- cml_seurat@meta.data %>% filter(effusion == "no effusion") %>% pull(barcode)
noeffusion_seurat <- subset(cml_seurat, cells = cells.to.keep) %>% getLatentUMAP()

deg_no_effusion <- lapply(unique(noeffusion_seurat$cluster), getDEGbyCluster, seurat_object = noeffusion_seurat)
deg_effusion    <- lapply(unique(effusion_seurat$cluster), getDEGbyCluster, seurat_object = effusion_seurat)

deg_no_effusion <- deg_no_effusion %>% rbindlist()
deg_effusion    <- deg_effusion %>% rbindlist()

deg_no_effusion <- deg_no_effusion %>% mutate(key = paste(gene,cluster,timepoint,direction))
deg_effusion <- deg_effusion %>% mutate(key = paste(gene,cluster,timepoint,direction))

deg_no_effusion_changes <- deg_no_effusion[!deg_no_effusion$key %in% deg_effusion$key, ]
deg_effusion_changes    <- deg_effusion[!deg_effusion$key %in% deg_no_effusion$key, ]

deg_no_effusion_changes %>% 
  group_by(timepoint,cluster) %>% summarise(n=n()) %>% 
  filter(timepoint == "12mov3mo") %>% 
  ggplot(aes(reorder(cluster,n),n)) + geom_point() + coord_flip() + facet_wrap(~timepoint)

deg_effusion_changes %>% 
  group_by(timepoint,cluster) %>% summarise(n=n()) %>% 
  filter(timepoint == "12mov3mo") %>% 
  ggplot(aes(reorder(cluster,n),n)) + geom_point() + coord_flip() + facet_wrap(~timepoint)

df <- deg_effusion_changes %>% filter(grepl("EMRA", cluster)) %>% 
  filter(timepoint == "12mov3mo") %>% arrange((p_val_adj))

df %>% 
  ggplot(aes(-avg_logFC,-log10(p_val_adj), label=gene)) + geom_point() + ggrepel::geom_text_repel(data = head(subset(df, !gene %in% unwanted_genes), 50), font.face = 3) + xlim(c(-5,5)) + facet_wrap(~cluster)





### Density map

## Density plot
cml_seurat$celltype <- do.call(cml_seurat$cluster %>% extractCoarsePhenotype(), what = "c")
cml_seurat$celltype <- ifelse(cml_seurat$celltype %in% c("B-cell", "CD4", "CD8", "Monocyte", "NK"), cml_seurat$celltype, NA)

median_df <- cml_seurat@reductions$umap@cell.embeddings %>% as.data.frame %>% bind_cols(cml_seurat@meta.data) %>% group_by(celltype) %>% summarise(umap1=median(UMAP_1), umap2=median(UMAP_2))

cml_seurat@reductions$umap@cell.embeddings %>% as.data.frame %>% bind_cols(cml_seurat@meta.data) %>%
  ggplot(aes(UMAP_1,UMAP_2)) + 
  stat_density_2d(aes(fill = ..level..), color="gray50", geom = "polygon") + 

  stat_ellipse(aes(group=celltype,color=celltype),linetype="dashed",type = "norm", segments = 10, size = 0.5) +
  facet_wrap(effusion~timepoint) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") +
  ggrepel::geom_label_repel(data=median_df,aes(umap1,umap2,label=celltype,color=celltype), min.segment.length = 1, nudge_y = 5) + 
  facets_nice + scale_color_manual(values = getPalette(9)) 
ggsave("results/manuscript/pe/umap_dens_timepoint.pdf", width = 9, height = 6)



eff_df <- cml_seurat@meta.data %>% group_by(orig.ident,patient.x,timepoint,effusion) %>% summarise(n=n()) %>% dplyr::select(-n)
  
cml_seurat@meta.data %>% group_by(orig.ident,cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(eff_df) %>% 
  ggplot(aes(effusion,prop,fill=effusion)) + geom_boxplot() + facet_wrap(cluster~timepoint,ncol=3, scales = "free_y") + scale_fill_manual(values = getPalette5(4)) + theme(legend.position = "none")
ggsave("results/manuscript/pe/box_timepoint.png", width = 6, height = 20)




median_df <- cd8_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(cd8_seurat@meta.data) %>% group_by(cluster) %>% summarise(umap1=median(latentumap_1), umap2=median(latentumap_2))

cd8_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(cd8_seurat@meta.data) %>%
  ggplot(aes(latentumap_1,latentumap_2)) + 
  stat_density_2d(aes(fill = ..level..), color="gray50", geom = "polygon") + 
  
  stat_ellipse(aes(group=cluster,color=cluster),linetype="dashed",type = "norm", segments = 10, size = 0.5) +
  facet_wrap(effusion~timepoint) + scale_fill_distiller(direction=1) + theme_bw(base_size = 12) + theme(legend.position = "right") + labs(x = "UMAP 1", y = "UMAP 2") +
  # ggrepel::geom_label_repel(data=median_df,aes(umap1,umap2,label=cluster,color=cluster), min.segment.length = 1, nudge_y = 5) + 
  facets_nice + scale_color_manual(values = getPalette(9)) 
ggsave("results/manuscript/pe/umap_cd8_dens_timepoint.pdf", width = 11, height = 6)


DimPlot(cd8_seurat)
