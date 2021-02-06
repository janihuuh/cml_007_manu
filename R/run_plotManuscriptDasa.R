
dir.create("results/manuscript/dasatinib/")

## Olink
olink_df <- fread("results/manuscript/olink_clean.txt") %>% dplyr::select(-V1) 
olink_dasa_df <- olink_df %>% melt(id = c("name", "timepoint")) %>% filter(timepoint %in% c("0mo", "3mo")) %>%   mutate(fill = ifelse(timepoint == "0mo", "diagnosis", "dasatinib+IFNa")) %>% 
  mutate(fill = ifelse(timepoint %in% c("3mo", "24mo"), "dasatinib", fill)) %>% mutate(type = "olink")

## Compare before and after dasatinib
facs_3mo <- fread("results/manuscript/facs_3mo_df.txt") 

dasa_facs_df <- facs_3mo %>% melt(id = c("timepoint","Study.nro","FM", "HRUH")) %>% 
  mutate(fill = ifelse(timepoint == "0mo", "diagnosis", "dasatinib+IFNa")) %>% 
  mutate(fill = ifelse(timepoint %in% c("3mo", "24mo"), "dasatinib", fill)) %>%
  filter(timepoint %in% c("0mo", "3mo")) %>% dplyr::select(-Study.nro, -HRUH) %>% dplyr::rename(name = FM) %>% dplyr::select(name,timepoint,variable,value,fill) %>% mutate(type = "facs")

dasa_facs_df <- dasa_facs_df %>% filter(!is.na(name))
dasa_df <- rbind(olink_dasa_df, dasa_facs_df)

p <- dasa_facs_df %>% ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + 
  # ggpubr::stat_compare_means(label = "p.signif") + 
  facet_wrap(~variable, scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)[c(2,1)]) + geom_jitter(size = 0.5)
ggsave(plot = p, "results/manuscript/dasatinib/violin_overall.pdf", width = 12, height = 12)


var_df <- dasa_df %>% group_by(variable, type) %>% summarise(n = n()) %>% dplyr::select(-n)
p.df <- lapply(unique(dasa_df$variable), FUN = function(x){
  message(x)
  y <- dasa_df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 2){
    median.x = y %>% filter(timepoint == "0mo") %>% pull(value) %>% median(na.rm = T)
    median.y = y %>% filter(timepoint == "3mo") %>% pull(value) %>% median(na.rm = T)
    wilcox.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x, median.0m = median.x, median.3m = median.y, dir = ifelse(log2(median.y/median.x) > 0, "up", "down"))
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) %>% left_join(var_df)
p.df$variable <- gsub("\\.", "\\_", p.df$variable)
fwrite(p.df, "results/manuscript/dasatinib/dasa_p_df.txt", sep = "\t", quote = F, row.names = F)
xlsx::write.xlsx(x = p.df, file = "results/manuscript/dasatinib/dasa_p_df.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)



p.df %>% 
p.df <- fread("results/manuscript/dasatinib/dasa_p_df.txt")


# p.df <- fread("results/manuscript/dasatinib/dasa_p_df.txt")

p.df %>% filter(p.adj < 0.1) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = type)) + geom_bar(stat = "identity") + coord_flip() + labs(x = "") + facet_wrap(~dir, scales = "free_y")
ggsave("results/manuscript/dasatinib/bar_sigf.pdf", width = 5, height = 6)

p.df %>% filter(type == "facs") %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  mutate(col = ifelse(p.adj < 0.05, dir, "no")) %>% 
  mutate(dir = plyr::revalue(dir, replace = c("down" = "up in diagnosis", "up" = "up in on dasatinib"))) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = col)) + 
  # geom_bar(stat = "identity") + 
  geom_segment(aes(x = reorder(variable,-log10(p.adj)), xend = reorder(variable,-log10(p.adj)), y=0, yend=-log10(p.adj)), color = "gray30") + geom_point(shape = 21, size = 5) +
  coord_flip() + labs(x = "") + facet_wrap(~dir) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + facets_nice + scale_fill_manual(values = c(getPalette(8)[2], "lightgrey", getPalette(8)[1])) + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/loliplot_facs_sigf.pdf", width = 5.5, height = 6)


p.df %>% filter(type == "facs") %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p.adj < 0.2) %>% 
  mutate(col = ifelse(p.adj < 0.05, dir, "no")) %>% 
  mutate(dir = plyr::revalue(dir, replace = c("down" = "up in diagnosis", "up" = "up in on dasatinib"))) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = col)) + 
  # geom_bar(stat = "identity") + 
  geom_segment(aes(x = reorder(variable,-log10(p.adj)), xend = reorder(variable,-log10(p.adj)), y=0, yend=-log10(p.adj)), color = "gray30") + geom_point(shape = 21, size = 5) +
  coord_flip() + labs(x = "") + facet_wrap(~dir) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + facets_nice + scale_fill_manual(values = c(getPalette(8)[2], "lightgrey", getPalette(8)[1])) + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/loliplot_facs_sigf_manu.pdf", width = 7, height = 6)



p.df %>% filter(type == "olink") %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  # filter(p.adj < 0.10) %>% 
  mutate(col = ifelse(p.adj < 0.05, dir, "no")) %>% 
  mutate(dir = plyr::revalue(dir, replace = c("down" = "up in diagnosis", "up" = "up in on dasatinib"))) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = col)) + 
  # geom_bar(stat = "identity") + 
  geom_segment(aes(x = reorder(variable,-log10(p.adj)), xend = reorder(variable,-log10(p.adj)), y=0, yend=-log10(p.adj)), color = "gray30") + geom_point(shape = 21, size = 5) +
  coord_flip() + labs(x = "") + facet_wrap(~dir) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + facets_nice + scale_fill_manual(values = c(getPalette(8)[2], "lightgrey", getPalette(8)[1])) + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/loliplot_olink_sigf.pdf", width = 5.5, height = 9)


p.df %>% filter(type == "olink") %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  filter(p.adj < 0.10) %>% 
  mutate(col = ifelse(p.adj < 0.05, dir, "no")) %>% 
  mutate(dir = plyr::revalue(dir, replace = c("down" = "up in diagnosis", "up" = "up in on dasatinib"))) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = col)) + 
  # geom_bar(stat = "identity") + 
  geom_segment(aes(x = reorder(variable,-log10(p.adj)), xend = reorder(variable,-log10(p.adj)), y=0, yend=-log10(p.adj)), color = "gray30") + geom_point(shape = 21, size = 5) +
  coord_flip() + labs(x = "") + facet_wrap(~dir) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + facets_nice + scale_fill_manual(values = c(getPalette(8)[2], "lightgrey", getPalette(8)[1])) + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/loliplot_olink_sigf_manu.pdf", width = 5.5, height = 4)





sigf_variables <- p.df %>% filter(p.adj < 0.1) %>% pull(variable)

dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  ggplot(aes(timepoint, value, fill = fill)) + geom_boxplot(draw_quantiles = 0.5, outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.signif", label.y.npc = 0.9, color = "black") + 
  theme_classic(base_size = 12) +
  facet_wrap(~variable, scales = "free_y") + 
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)[c(2,1)]) + geom_jitter(size = 0.5, color = "darkgrey") + facets_nice + theme(legend.position = "none") + labs(x = "")
ggsave("results/manuscript/dasatinib/violin_facs_sigf.pdf", width = 9, height = 6)

sigf_variables <- p.df %>% filter(p.adj < 0.1)%>% filter(type == "olink")  %>% pull(variable) 

dasa_df$variable <- gsub("\\.", "\\_", dasa_df$variable)

dasa_df %>% filter(variable %in% sigf_variables) %>% 
  ggplot(aes(timepoint, value)) + geom_boxplot(aes(fill = fill), draw_quantiles = 0.5, outlier.shape = NA) + 
  theme_classic(base_size = 12) +
  ggpubr::stat_compare_means(label = "p.signif", label.y.npc = 0.8, color = "black") + 
  facet_wrap(~variable, scales = "free_y") + 
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)[c(2,1)]) + geom_jitter(size = 0.5, color = "darkgrey") + facets_nice + theme(legend.position = "none") + labs(x = "")
ggsave("results/manuscript/dasatinib/violin_olink_sigf.pdf", width = 9, height = 6)




#### NK FACS
sigf_variables <- p.df %>% filter(p.adj < 0.1) %>% pull(variable)

dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  # filter(variable %in% c("NK.abs", "NK_CD16",	"NK_CD27",	"NK_CD45RA",	"NK_CD56BRIGHT",	"NK_CD57", "NK_CD62L", "NK_GrB")) %>% 
  filter(variable %in% c("NK_CD16",	"NK_CD27",	"NK_CD45RA",	"NK_CD57", "NK_CD62L", "NK_GrB")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + 
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot(draw_quantiles = 0.5, outlier.shape = NA) + 
  # facet_wrap(~variable, scales = "free_y", ncol=3) + 
  facet_wrap(~variable, ncol=3) + 
  theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette(8)) + geom_jitter(size = 0.1, color = "black") + facets_nice +
  ggpubr::stat_compare_means(label = "p.signif", label.y.npc = 0.9) + 
  # ggsignif::geom_signif(comparisons = list(c("0mo", "3mo")), map_signif_level = T) +
  labs(x = "", y = "prop of NK cells (%)") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/box_nk_facs_sigf.pdf", width = 4, height = 3)


dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("NK_abs")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + 
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot(draw_quantiles = 0.5, outlier.shape = NA) + 
  # facet_wrap(~variable, scales = "free_y", ncol=3) + 
  facet_wrap(~variable, ncol=3) + 
  theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette(8)) + geom_jitter(size = 0.1, color = "black") + facets_nice +
  ggpubr::stat_compare_means(label = "p.signif", label.y.npc = 0.9) + 
  # ggsignif::geom_signif(comparisons = list(c("0mo", "3mo")), map_signif_level = T) +
  labs(x = "", y = "abs number (10e9/l)") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/box_nk_abs_facs_sigf.pdf", width = 1.5, height = 2)




dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("CD3_abs", "CD8_abs", "CD8_CD57_abs")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + 
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot(draw_quantiles = 0.5, outlier.shape = NA) + 
  facet_wrap(~variable, scales = "free_y", ncol=3) + theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette(8)) + geom_jitter(size = 0.1) + facets_nice +
  ggpubr::stat_compare_means(label = "p.signif", label.y.npc = 0.9) + 
  # ggsignif::geom_signif(comparisons = list(c("0mo", "3mo")), map_signif_level = T) +
  labs(x = "", y = "abs number (10e9/l)") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/box_cd8_facs_sigf.pdf", width = 5, height = 3)



dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("NKCD16",	"NKCD27",	"NKCD45RA",	"CD56BRIGHT",	"NKCD57", "NKCD62L", "NK_GrB")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + geom_path(aes(group=name)) +
  # facet_wrap(~reorder(variable, p.adj), scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  facet_wrap(~variable, scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5) + facets_nice +
  labs(x = "", y = "prop of NK cells")



dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("NK", "NK.abs")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + 
  facet_wrap(~variable, scales = "free_y") +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5) + facets_nice +
  labs(x = "") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/violin_nk_only.pdf", width = 4, height = 3)



#### CD8 FACS
dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% filter(grepl("CD8", variable) | grepl("CD3.abs", variable)) %>% 
  filter(variable != "CD8.EM") %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 

  ggplot(aes(timepoint, value, fill = fill)) + 
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot(draw_quantiles = 0.5, outlier.shape = NA) + 
  facet_wrap(~variable, scales = "free_y", ncol=4) + theme_classic(base_size = 12) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette(8)) + geom_jitter(size = 0.1) + facets_nice +
  # ggpubr::stat_compare_means(label = "p.signif") + 
  # ggsignif::geom_signif(comparisons = list(c("0mo", "3mo")), map_signif_level = T) +
  labs(x = "", y = "value") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/box_cd8_facs_sigf.pdf", width = 5, height = 3)

dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("NK", "NK.abs")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + 
  facet_wrap(~variable, scales = "free_y") +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5) + facets_nice +
  labs(x = "") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/violin_nk_only.pdf", width = 4, height = 3)



### OLINK results
dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "olink") %>% 
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + 
  # facet_wrap(~reorder(variable, p.adj), scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  facet_wrap(~variable, scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5) + facets_nice
ggsave("results/manuscript/dasatinib/violin_olink_sigf.pdf", width = 7, height = 4)

volcano_df <- dasa_df %>% left_join(p.df) %>% 
  left_join(var_df) %>% filter(type == "olink") %>% mutate(log2fc = log2(median.3m / median.0m))

label_df <- volcano_df %>% dplyr::select(-name,-timepoint,-value) %>% distinct()
label_df <- label_df[!duplicated(label_df$variable), ] %>% mutate(p.adj = p.adjust(p.value, method = "BH"))

label_df %>% 
  ggplot(aes(log2fc,-log10(p.value), color = ifelse(p.adj < 0.05, "sigf", "unsigf"))) + geom_point() + xlim(values = c(-5,5)) + theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = subset(label_df, p.adj < 0.05), aes(label = variable)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + geom_vline(xintercept = -1, linetype = "dotted") + geom_vline(xintercept = 1, linetype = "dotted")





## Focus only on these two time point
cells.to.keep  <- cml_seurat@meta.data %>% filter(timepoint %in% c("dg", "3mo")) %>% pull(barcode)
dasa_seurat    <- subset(cml_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat()
dasa_seurat    <- RunUMAP(dasa_seurat, reduction = "latent", dims = 1:30)

DimPlot(dasa_seurat, reduction = "umap", label = T, repel = T) #, split.by = "timepoint")





### scRNAseq; DEGs
df <- DEG_cluster_df %>% 
  filter(timepoint == "3movdg") %>% 
  mutate(celltype = extractCoarsePhenotype(cluster)) %>% mutate(cluster = cluster %>% reorderClusters()) %>% 
  group_by(celltype,cluster, timepoint) %>% summarise(n=n()) %>% ungroup() %>% 
  mutate(celltype = ifelse(celltype %in% c("cycling", "Cytotoxic", "pDC"), "other", celltype)) %>% 
  filter(timepoint %in% c("3movdg", "12mov3mo")) %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("3movdg", "12mov3mo")))

ggplot(df, aes(reorder(cluster,n),n,fill=celltype)) + geom_segment(aes(x=reorder(cluster,n),xend=reorder(cluster,n),y=0,yend=n)) + geom_point(shape=21,size=5) + 
  scale_fill_manual(values=getPalette(8)) + facets_nice + labs(x = "", y = "DEGs") + coord_flip() #+ ggrepel::geom_text_repel(data = subset(df, timepoint == "3movdg"), aes(label=cluster), nudge_y = 100) #+ theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/lolliplot_deg.pdf", width = 6, height = 5)


nk_deg <- DEG_cluster_df %>% filter(timepoint == "3movdg")  %>% filter(cluster %in% c("0 NK CD56dim", "3 NK CD56dim IFNg cytotoxic", "4 NK CD56dim IFNg exhausted", "17 NK CD56bright")) 
nk_deg %>% filter(gene %in% cytotoxic_markers)

dasa_seurat <- AddModuleScore(dasa_seurat, features = list(c(cytotoxic_markers[-c(1,10)], "NKG7", "FCGR3A")), name = "cytotoxicity")

dasa_seurat@meta.data %>% filter(cluster %in% c("0 NK CD56dim", "3 NK CD56dim IFNg cytotoxic","17 NK CD56bright")) %>% 
  ggplot(aes(timepoint,cytotoxicity1, fill = timepoint)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1) + facet_wrap(~cluster) + scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") +
  ggpubr::stat_compare_means(label = "p.format") + facets_nice
ggsave("results/manuscript/dasatinib/vln_cytotoxicity.pdf", width = 6, height = 5)


dasa_seurat@meta.data %>% filter(cluster %in% c("0 NK CD56dim", "3 NK CD56dim IFNg cytotoxic","17 NK CD56bright")) %>% 
  ggplot(aes(timepoint,cytotoxicity1, fill = timepoint)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1) + scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") +
  ggpubr::stat_compare_means(label = "p.format") + facets_nice
ggsave("results/manuscript/dasatinib/vln_cytotoxicity.pdf", width = 6, height = 5)





# Get enrichment
universe_df <- rownames(cml_seurat)

up_clusters <- lapply(unique(DEG_cluster_df$cluster), function(x){
  p <- DEG_cluster_df %>% filter(timepoint == "3movdg" & cluster == x & direction == "up") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(p)){
    p %>% dplyr::mutate(cluster = x)
  }
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% arrange(p.adjust)

fwrite(up_clusters, "results/manuscript/dasatinib/enrichment_hallmark_up.txt", sep = "\t", quote = F, row.names = F)

down_clusters <- lapply(unique(subset(DEG_cluster_df,timepoint == "3movdg" & direction == "down")$cluster), function(x){
  message(x)
  df <- DEG_cluster_df %>% filter(timepoint == "3movdg" & cluster == x & direction == "down") 
  if(!is.null(df)){
    p <- try(df %>% getHypergeometric(universe_df = universe_df, term_df = hallmark))
  }
  if(!is.null(p)){
    p %>% dplyr::mutate(cluster = x)
  }
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% arrange(p.adjust)

fwrite(down_clusters, "results/manuscript/dasatinib/enrichment_hallmark_dn.txt", sep = "\t", quote = F, row.names = F)



## Focus on NK-cells 
cells.to.keep  <- cml_seurat@meta.data %>% filter(timepoint %in% c("dg", "3mo") & cluster %in% c("0 NK CD56dim", "3 NK CD56dim IFNg cytotoxic", "4 NK CD56dim IFNg exhausted", "17 NK CD56bright")) %>% pull(barcode)
# cells.to.keep  <- cml_seurat@meta.data %>% filter(timepoint %in% c("dg", "3mo") & cluster %in% c("0 NK CD56dim", "3 NK CD56dim IFNg cytotoxic", "17 NK CD56bright")) %>% pull(barcode)
nk_dasa_seurat <- subset(cml_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat()
nk_dasa_seurat <- nk_dasa_seurat %>% getLatentClustering()

q <- NULL; i <- 1
for(clustering_column in clustering_columns){
  q[[i]] <- nk_dasa_seurat@meta.data[,clustering_column] %>% levels %>% length; i <- i + 1
}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/manuscript/dasatinib/scatter_res_nClusters.png", width = 4, height = 3)

DimPlot(nk_dasa_seurat, label = T, repel = T, group.by = "cluster", cols = getPalette3(4)) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/umap_nk.png", width = 5, height = 4)

DimPlot(nk_dasa_seurat, label = T, repel = T, group.by = "RNA_snn_res.0.4", cols = getPalette3(5)) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")
nk_dasa_seurat$reclustered <- nk_dasa_seurat$RNA_snn_res.0.4
Idents(nk_dasa_seurat) <- nk_dasa_seurat$reclustered

nk_dasa_deg <- FindAllMarkers(nk_dasa_seurat, test.use = "t")
nk_dasa_deg <- nk_dasa_deg %>% filter(p_val_adj < 0.05) %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down")) %>% filter(dir == "up")
nk_dasa_deg %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% View


nk_dasa_orig_deg <- FindAllMarkers(nk_dasa_seurat, test.use = "t") %>% filter(p_val_adj < 0.05) %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down")) %>% filter(dir == "up")
nk_dasa_orig_deg %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% View


DotPlot(nk_dasa_seurat, features = rev(c("FCGR3A", "CD27", "PTPRC",  "NCAM1", "B3GAT1", "SELL", "KLRG1")), cols = "RdYlBu") + ggpubr::rotate_x_text(90)

Idents(nk_dasa_seurat) <- nk_dasa_seurat$cluster
DotPlot(nk_dasa_seurat, features = rev(smith_mark), cols = "RdYlBu") + ggpubr::rotate_x_text(90)
ggsave("results/manuscript/dasatinib/dot_nk.pdf", width = 12, height = 4)


getClusterPhenotypes <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 NK CD56dim",
                                                    "1"  = "1 NK CD56dim CD16high",
                                                    "2"  = "2 NK CD56dim",
                                                    "3"  = "3 NK CD56dim",
                                                    "4"  = "4 NK CD56bright"))
                                              
  return(clusters)
}

df        <- nk_dasa_seurat@meta.data %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
median_df <- df %>% group_by(timepoint, cluster) %>% summarise(median = median(prop)) %>% top_n(n = 5, wt = median) %>% arrange(desc(median))

umap_means_df <- nk_dasa_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(nk_dasa_seurat@meta.data) %>% group_by(cluster) %>% summarise(latent_umap_1 = median(latent_umap_1), latent_umap_2 = median(latent_umap_2)) 
umap_means_df <- median_df %>% left_join(umap_means_df)

nk_dasa_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(nk_dasa_seurat@meta.data) %>%
  #  filter(!is.na(orig.ident)) %>% filter(orig.ident != "NA") %>% mutate(orig.ident = droplevels(as.factor(orig.ident))) %>% 
  ggplot(aes(latent_umap_1, latent_umap_2, color = timepoint)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~timepoint, ncol = 7) + 
  scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  ggrepel::geom_text_repel(data = umap_means_df, aes(latent_umap_1,latent_umap_2,label=cluster), color = "black", size = 3.5) +
  scale_color_manual(values = getPalette3(3)) + guides(color=FALSE) + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/umap_nk_dens.png", width = 8, height = 4)



nk_degs <- DEG_cluster_df %>% filter(timepoint == "3movdg") %>% filter(cluster %in% c("0 NK CD56dim", "3 NK CD56dim IFNg cytotoxic", "17 NK CD56bright")) 
nk_degs %>% filter(cluster == "0 NK CD56dim") %>% filter(direction == "up")





## Pseudotime
require(slingshot)
require(SingleCellExperiment)
require(SummarizedExperiment)


nk_dasa_seurat <- RunUMAP(nk_dasa_seurat, reduction = "latent", dims = 1:30, n.neighbors = 5, min.dist = 0.3)
DimPlot(nk_dasa_seurat, reduction = "umap")


## Total
nk_dasa_seurat$orig.clusters <- nk_dasa_seurat$cluster
nk_dasa_sce   <- as.SingleCellExperiment(nk_dasa_seurat)
nk_dasa_sling <- slingshot(data = nk_dasa_sce, clusterLabels = as.character(nk_dasa_sce$orig.clusters), reducedDim = "UMAP", start.clus = "4 NK CD56dim IFNg exhausted")
sling_curve   <- SlingshotDataSet(nk_dasa_sling)

plot(reducedDims(nk_dasa_sling)[["LATENT_UMAP"]], col = cololors, pch = 16, asp = 1)
lines(sling_curve@curves[[1]], lwd = 4, col = curve_cols[1])
lines(sling_curve@curves[[2]], lwd = 4, col = curve_cols[2])

nk_dasa_seurat$pseudotime <- nk_dasa_sling$slingPseudotime_1
nk_dasa_seurat <- AddModuleScore(nk_dasa_seurat, features = list(c("GZMA", "GZMB", "GZMH", "GZMK", "GNLY", "PRF1", "CTSW")), name = "cytotoxicity")
nk_dasa_seurat <- AddModuleScore(nk_dasa_seurat, features = list(c("TIGIT", "KLRD1", "KLRB1", "LAIR1", "KLRC1", "TNFSF10", "FASLG")), name = "inhibitory")
nk_dasa_seurat <- AddModuleScore(nk_dasa_seurat, features = list(c("FCGR3A", "NCR1", "NCR3", "KLRK1", "KLRC2", "KLRF1", "CD226", "CD244")), name = "activation")
nk_dasa_seurat <- AddModuleScore(nk_dasa_seurat, features = list(c("IFNG", "LTB", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2")), name = "cytokines")
nk_dasa_seurat <- AddModuleScore(nk_dasa_seurat, features = list(c("IL7R", "IL2RB", "CXCR3", "CXCR4", "CCR7", "CXCR1", "S1PR5")), name = "cytokine_receptors")
nk_dasa_seurat <- AddModuleScore(nk_dasa_seurat, features = list(c("ITGAL", "CD2", "CD58", "ITGB2", "SELL", "CD96")), name = "activation")



nk_dasa_seurat@meta.data %>% 
  ggplot(aes(cluster,cytotoxicity1,color=cluster)) +geom_violin()

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(cluster,inhibitory1,color=cluster)) +geom_violin()

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(cluster,activation1,color=cluster)) +geom_violin()

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(cluster,cytokines1,color=cluster)) +geom_violin()




nk_dasa_seurat@meta.data %>% 
  ggplot(aes(pseudotime,cytotoxicity1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(pseudotime,inhibitory1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(pseudotime,activation1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(pseudotime,cytokines1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat@meta.data %>% 
  ggplot(aes(pseudotime,activation1,color=cluster)) + geom_point() + geom_smooth(color = "black")






## Baseline
cells.to.keep     <- nk_dasa_seurat@meta.data %>% filter(timepoint == "dg") %>% pull(barcode)
nk_dasa_seurat_bl <- subset(nk_dasa_seurat, cells = cells.to.keep)
nk_dasa_seurat_bl$orig.clusters <- nk_dasa_seurat_bl$cluster
nk_dasa_bl_sce   <- as.SingleCellExperiment(nk_dasa_seurat_bl)
nk_dasa_bl_sling <- slingshot(data = nk_dasa_bl_sce, clusterLabels = as.character(nk_dasa_bl_sce$orig.clusters), reducedDim = "LATENT_UMAP") #, start.clus = "17 NK CD56bright")
sling_bl_curve   <- SlingshotDataSet(nk_dasa_bl_sling)

## Dasa
cells.to.keep     <- nk_dasa_seurat@meta.data %>% filter(timepoint != "dg") %>% pull(barcode)
nk_dasa_seurat_3m <- subset(nk_dasa_seurat, cells = cells.to.keep)
nk_dasa_seurat_3m$orig.clusters <- nk_dasa_seurat_3m$cluster
nk_dasa_3m_sce   <- as.SingleCellExperiment(nk_dasa_seurat_3m)
nk_dasa_3m_sling <- slingshot(data = nk_dasa_3m_sce, clusterLabels = as.character(nk_dasa_3m_sce$orig.clusters), reducedDim = "LATENT_UMAP", start.clus = "17 NK CD56bright")
sling_3m_curve   <- SlingshotDataSet(nk_dasa_3m_sling)


cololors                 <- nk_dasa_sling$orig.clusters %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 0]  <- getPalette3(4)[1]
cololors[cololors == 3]  <- getPalette3(4)[2]
cololors[cololors == 17] <- getPalette3(4)[3]

curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")

pdf("results/manuscript/dasatinib/umap_sling_baseline.pdf", width = 7, height = 4)
plot(reducedDims(nk_dasa_sling)[["LATENT_UMAP"]], col = cololors, pch = 16, asp = 1)
plot(reducedDims(nk_dasa_bl_sling)[["LATENT_UMAP"]], col = cololors, pch = 16, asp = 1)

lines(sling_bl_curve@curves[[1]], lwd = 4, col = curve_cols[1])
# lines(sling_3m_curve@curves[[1]], lwd = 4, col = curve_cols[2])
dev.off()


png("results/public/umap_sling_nk.png", width = 5, height = 5, units = "in", res = 1024)
plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()



nk_dasa_seurat_bl$pseudotime <- nk_dasa_bl_sling$slingPseudotime_1
nk_dasa_seurat_bl <- AddModuleScore(nk_dasa_seurat_bl, features = list(c("GZMA", "GZMB", "GZMH", "GZMK", "GNLY", "PRF1", "CTSW")), name = "cytotoxicity")
nk_dasa_seurat_bl <- AddModuleScore(nk_dasa_seurat_bl, features = list(c("TIGIT", "KLRD1", "KLRB1", "LAIR1", "KLRC1", "TNFSF10", "FASLG")), name = "inhibitory")
nk_dasa_seurat_bl <- AddModuleScore(nk_dasa_seurat_bl, features = list(c("FCGR3A", "NCR1", "NCR3", "KLRK1", "KLRC2", "KLRF1", "CD226", "CD244")), name = "activation")
nk_dasa_seurat_bl <- AddModuleScore(nk_dasa_seurat_bl, features = list(c("IFNG", "LTB", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2")), name = "cytokines")
nk_dasa_seurat_bl <- AddModuleScore(nk_dasa_seurat_bl, features = list(c("IL7R", "IL2RB", "CXCR3", "CXCR4", "CCR7", "CXCR1", "S1PR5")), name = "cytokine_receptors")
nk_dasa_seurat_bl <- AddModuleScore(nk_dasa_seurat_bl, features = list(c("ITGAL", "CD2", "CD58", "ITGB2", "SELL", "CD96")), name = "activation")

nk_dasa_seurat_bl@meta.data %>% 
  ggplot(aes(pseudotime,cytotoxicity1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat_bl@meta.data %>% 
  ggplot(aes(pseudotime,inhibitory1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat_bl@meta.data %>% 
  ggplot(aes(pseudotime,activation1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat_bl@meta.data %>% 
  ggplot(aes(pseudotime,cytokines1,color=cluster)) + geom_point() + geom_smooth(color = "black")

nk_dasa_seurat_bl@meta.data %>% 
  ggplot(aes(pseudotime,activation1,color=cluster)) + geom_point() + geom_smooth(color = "black")




## Focus on B-cells


#### B FACS
dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("B.cells")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + 
  # facet_wrap(~reorder(variable, p.adj), scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  facet_wrap(~variable, scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5) + facets_nice +
  labs(x = "", y = "prop of NK cells") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/violin_b_facs_sigf.pdf", width = 2, height = 2)



cells.to.keep  <- cml_seurat@meta.data %>% filter(timepoint %in% c("dg", "3mo") & cluster %in% c("9 B-cell naive", "11 B-cell memory")) %>% pull(barcode)
b_dasa_seurat <- subset(cml_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat() %>% getLatentClustering()

q <- NULL; i <- 1
for(clustering_column in clustering_columns){
  q[[i]] <- b_dasa_seurat@meta.data[,clustering_column] %>% levels %>% length; i <- i + 1
}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/manuscript/dasatinib/scatter_b_res_nClusters.png", width = 4, height = 3)

DimPlot(b_dasa_seurat, label = T, repel = T, group.by = "cluster", split.by = "timepoint", cols = getPalette(8)) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/umap_b.png", width = 6, height = 3)


down_clusters %>% filter(cluster %in% c("9 B-cell naive", "11 B-cell memory")) %>% 
  mutate(ID = gsub("HALLMARK\\_", "", ID)) %>% 
  ggplot(aes(reorder(ID, -log10(p.adjust)),-log10(p.adjust))) + geom_bar(stat = "identity", fill = "tomato") + coord_flip() + labs(x = "")
ggsave("results/manuscript/dasatinib/bar_b_pathways.pdf", width = 5, height = 3)



b_dasa_seurat@meta.data %>% group_by(cluster, timepoint, patient) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% 
  ggplot(aes(cluster,prop)) + geom_boxplot() + facet_wrap(~timepoint) + ggpubr::stat_compare_means() #+ geom_path(aes(group=patient)) #+ ggpubr::stat_compare_means()

DimPlot(b_dasa_seurat, label = T, repel = T, group.by = "RNA_snn_res.0.7", split.by = "timepoint", cols = getPalette3(8)) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")
Idents(b_dasa_seurat) <- b_dasa_seurat$RNA_snn_res.0.7

b_dasa_deg <- FindAllMarkers(b_dasa_seurat, test.use = "t") %>% filter(p_val_adj < 0.05) %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down")) %>% filter(dir == "up")
b_dasa_deg %>% filter(cluster == "5") %>% View



DotPlot(b_dasa_seurat, features = rev(c("CD19", "CD20", "IGHM", "IGHG", "CCR7", "CD27", "CD1C", "CD24", "NURR77", "CD99", "CD38", "BCL6", "CD14")), cols = "RdYlBu")  +
  ggpubr::rotate_x_text(90) + theme(legend.position = "top")












## Focus on cd8-cells

#### cd8+ FACS
dasa_df %>% filter(variable %in% sigf_variables) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>% 
  left_join(var_df) %>% filter(type == "facs") %>% 
  filter(variable %in% c("CD3.abs", "CD8.abs")) %>% 
  
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + 
  # facet_wrap(~reorder(variable, p.adj), scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  facet_wrap(~variable, scales = "free_y") + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5) + facets_nice +
  labs(x = "", y = "prop of NK cells") + theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/violin_t_facs_sigf.pdf", width = 4, height = 4)



cells.to.keep   <- cml_seurat@meta.data %>% filter(timepoint %in% c("dg", "3mo") & cluster %in% c("1 CD8 effector/exhausted", "2 CD8 EM", "6 CD8 EMRA", "12 CD8 naive/CM")) %>% pull(barcode)
cd8_dasa_seurat <- subset(cml_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat() 
cd8_dasa_seurat <- cd8_dasa_seurat %>% getLatentClustering()

q <- NULL; i <- 1
for(clustering_column in clustering_columns){
  q[[i]] <- cd8_dasa_seurat@meta.data[,clustering_column] %>% levels %>% length; i <- i + 1
}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/manuscript/dasatinib/scatter_b_res_nClusters.png", width = 4, height = 3)


DimPlot(cd8_dasa_seurat, label = T, repel = T, group.by = "cluster", split.by = "timepoint", cols = getPalette3(8)) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")


heatmap_up_df <- up_clusters %>% filter(cluster %in% c("1 CD8 effector/exhausted", "2 CD8 EM", "6 CD8 EMRA", "12 CD8 naive/CM")) %>% 
  mutate(ID = gsub("HALLMARK\\_", "", ID)) %>% group_by(ID, cluster) %>% mutate(p.adj = -log10(p.adjust))

heatmap_dn_df <- down_clusters %>% filter(cluster %in% c("1 CD8 effector/exhausted", "2 CD8 EM", "6 CD8 EMRA", "12 CD8 naive/CM")) %>% 
  mutate(ID = gsub("HALLMARK\\_", "", ID)) %>% group_by(ID, cluster) %>% mutate(p.adj = -log10(p.adjust))

  
heatmap_df <- rbind(heatmap_dn_df, heatmap_up_df)
heatmap_df <- dcast(ID~cluster,data=heatmap_df,value.var="p.adj")
rownames(heatmap_df) <- heatmap_df$ID
heatmap_df <- heatmap_df[,-1]
heatmap_df[is.na(heatmap_df)] <- 0
heatmap_df <- as.matrix(heatmap_df)

p <- pheatmap::pheatmap(heatmap_df)
ggsave(plot = print(p), "results/manuscript/dasatinib/heatmap_cd8_up.pdf", width = 6, height = 5)






DimPlot(cd8_dasa_seurat, label = T, repel = T, group.by = "RNA_snn_res.0.2", split.by = "timepoint", cols = getPalette3(8)) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")
Idents(cd8_dasa_seurat) <- cd8_dasa_seurat$RNA_snn_res.0.2

zhang_cd8_markers
DotPlot(cd8_dasa_seurat, features = (zhang_cd8_genes), cols = "RdYlBu") + ggpubr::rotate_x_text(90)


getClusterPhenotypesDasaCD8 <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 ",
                                                    "1"  = "1 ",
                                                    "2"  = "2 ",
                                                    "3"  = "3 ",
                                                    "4"  = "4 "))
  
  return(clusters)
}

b_dasa_deg <- FindAllMarkers(b_dasa_seurat, test.use = "t") %>% filter(p_val_adj < 0.05) %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down")) %>% filter(dir == "up")
b_dasa_deg %>% filter(cluster == "5") %>% View



DotPlot(b_dasa_seurat, features = rev(c("CD19", "CD20", "IGHM", "IGHG", "CCR7", "CD27", "CD1C", "CD24", "NURR77", "CD99", "CD38", "BCL6", "CD14")), cols = "RdYlBu")  +
  ggpubr::rotate_x_text(90) + theme(legend.position = "top")









## scrnaseq

#######################

df <- deg_all %>% 
  filter(timepoint == "3movdg") %>% 
  filter(cluster == "5 CD8 EMRA ZNF683+") %>% arrange(desc(abs(avg_logFC)))

df %>% 
  ggplot(aes(-avg_logFC,-log10(p_val_adj), label=gene)) + geom_point() + ggrepel::geom_text_repel(data = head(subset(df, !gene %in% unwanted_genes), 30), font.face = 3) + xlim(c(-100,100))



table(deg_all$cluster)
df <- deg_all %>% 
  filter(timepoint == "3movdg") %>% 
  filter(grepl("Monocyte", cluster)) %>% arrange(desc(abs(avg_logFC)))

df %>% 
  ggplot(aes(-avg_logFC,-log10(p_val_adj), label=gene)) + geom_point() + ggrepel::geom_text_repel(data = head(subset(df, !gene %in% unwanted_genes), 30), font.face = 3) + xlim(c(-100,100)) + facet_wrap(~cluster)

table(cml_seurat$cluster)
cells.to.keep <- cml_seurat@meta.data %>% filter(grepl("Monocyte", cluster) | grepl("pDC", cluster)) %>% pull(barcode)
mono_seurat   <- subset(cml_seurat, cells = cells.to.keep)  %>% getLatentUMAP() %>% AddModuleScore(features = list(grep("^HLA\\-D", rownames(cml_seurat), value = T)), name = "hla_ii", nbin = 15)

mono_seurat@meta.data %>% 
  filter(timepoint %in% c("dg", "3mo")) %>% 
  filter(cluster == "6 Monocyte CD16- HAVCR2+") %>% 
  ggplot(aes(timepoint,hla_ii1,fill=timepoint)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(~cluster) + theme_classic(base_size = 12) + facets_nice + scale_fill_manual(values = getPalette(8)) + labs(x = "", y = "HLA class II score") +
  theme(legend.position = "none")
ggsave("results/manuscript/dasatinib/violin_monocytes_hla.pdf", width = 3, height = 4)
