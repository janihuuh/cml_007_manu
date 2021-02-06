
facs_0m <- fread("data/facs/007_0mo.txt") 
colnames(facs_0m) <- make.names(make.unique(colnames(facs_0m)))
facs_0m <- facs_0m %>% mutate(timepoint = "0mo")
olink.cols <- grep("^X", colnames(facs_0m), value = T)
empty.cols <- grep("^V", colnames(facs_0m), value = T)
facs_0m <- facs_0m %>% dplyr::select(-olink.cols) %>% dplyr::select(-empty.cols)

facs_3m <- fread("data/facs/007_3mo.txt") 
colnames(facs_3m) <- make.names(make.unique(colnames(facs_3m)))
facs_3m <- facs_3m %>% mutate(timepoint = "3mo")
olink.cols <- grep("^X", colnames(facs_3m), value = T)
empty.cols <- grep("^V", colnames(facs_3m), value = T)
facs_3m <- facs_3m %>% dplyr::select(-olink.cols) %>% dplyr::select(-empty.cols)

facs_12m <- fread("data/facs/007_12mo.txt") 
colnames(facs_12m) <- make.names(make.unique(colnames(facs_12m)))
facs_12m <- facs_12m %>% mutate(timepoint = "12mo")
olink.cols <- grep("^X", colnames(facs_12m), value = T)
empty.cols <- grep("^V", colnames(facs_12m), value = T)
facs_12m <- facs_12m %>% dplyr::select(-olink.cols) %>% dplyr::select(-empty.cols)

facs_24m <- fread("data/facs/007_24mo.txt") 
colnames(facs_24m) <- make.names(make.unique(colnames(facs_24m)))
facs_24m <- facs_24m %>% mutate(timepoint = "24mo")
olink.cols <- grep("^X", colnames(facs_24m), value = T)
empty.cols <- grep("^V", colnames(facs_24m), value = T)
facs_24m <- facs_24m %>% dplyr::select(-olink.cols) %>% dplyr::select(-empty.cols)

features.to.use.3mo  <- Reduce(intersect,  list(colnames(facs_0m), colnames(facs_3m)))
features.to.use.12mo <- Reduce(intersect,  list(colnames(facs_0m), colnames(facs_3m), colnames(facs_12m)))
features.to.use.24mo <- Reduce(intersect,  list(colnames(facs_0m), colnames(facs_3m), colnames(facs_12m), colnames(facs_24m)))

facs_3mo <- rbind(facs_0m %>% dplyr::select(features.to.use.3mo), facs_3m %>% dplyr::select(features.to.use.3mo)) %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("0mo", "3mo")))
fwrite(facs_3mo, "results/manuscript/facs_3mo_df.txt", sep = "\t", quote = F, row.names = F)

facs_12mo <- rbind(facs_0m %>% dplyr::select(features.to.use.12mo), facs_3m %>% dplyr::select(features.to.use.12mo), facs_12m %>% dplyr::select(features.to.use.12mo)) %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("0mo", "3mo", "12mo")))
fwrite(facs_12mo, "results/manuscript/facs_12mo_df.txt", sep = "\t", quote = F, row.names = F)

facs_24mo <- rbind(facs_12m %>% dplyr::select(features.to.use.12mo), facs_24m %>% dplyr::select(features.to.use.24mo)) %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("12mo", "24mo")))
fwrite(facs_24mo, "results/manuscript/facs_24mo_df.txt", sep = "\t", quote = F, row.names = F)

facs_24mo <- rbind(facs_0m %>% dplyr::select(features.to.use.24mo), facs_3m %>% dplyr::select(features.to.use.24mo), facs_12m %>% dplyr::select(features.to.use.24mo), facs_24m %>% dplyr::select(features.to.use.24mo)) %>% 
  mutate(timepoint = factor(as.character(timepoint), levels = c("0mo", "3mo", "12mo", "24mo")))
fwrite(facs_24mo, "results/manuscript/facs_df.txt", sep = "\t", quote = F, row.names = F)

facs_24mo %>% melt(id = c("timepoint","Study.nro","FM", "HRUH")) %>% 
  mutate(fill = ifelse(timepoint == "0mo", "diagnosis", "dasatinib+IFNa")) %>% 
  mutate(fill = ifelse(timepoint %in% c("3mo", "24mo"), "dasatinib", fill)) %>%
  
  ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + facet_wrap(~variable) + # ggsignif::geom_signif(comparisons = list(c("0mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5)
ggsave("results/manuscript/overall/violin_overall.pdf", width = 12, height = 12)







######## PCA time

facs_full      <- fread("results/manuscript/facs_df.txt")
facs_full_mtx  <- facs_full %>% dplyr::select(-c(timepoint:HRUH))
facs_full_clin <- facs_full %>% dplyr::select(c(timepoint:HRUH))

## Remove rows with over 75% NA
ind1 <- rowSums(is.na(facs_full_mtx)) != ncol(facs_full_mtx) 
ind2 <- apply(facs_full_mtx, 1, function(x)table(is.na(x))[1]) / 36 > 0.75

facs_full_mtx <- facs_full_mtx[intersect(which(ind1),which(ind2)), ]
facs_full_mtx <- facs_full_mtx %>% as.data.frame()

table(is.na(facs_full_mtx))

facs_full_clin <- facs_full_clin[intersect(which(ind1),which(ind2)), ]
facs_full_clin <- facs_full_clin %>% as.data.frame()

## Impute with col means
for(i in 1:ncol(facs_full_mtx)){
  facs_full_mtx[is.na(facs_full_mtx[,i]), i] <- mean(facs_full_mtx[,i], na.rm = TRUE)
}

facs_full_pca <- prcomp((facs_full_mtx))
pca_df <- facs_full_clin %>% bind_cols(as.data.frame(facs_full_pca$x)[,1:10]) %>% bind_cols(facs_full_mtx)

pca_df$timepoint <- factor(as.character(pca_df$timepoint), levels = c("0mo", "3mo", "12mo", "24mo"))
pca_df$timepoint <- plyr::revalue(pca_df$timepoint, replace = c("0mo"  = "0m",
                                                                "3mo"  = "3m (dasa)",
                                                                "12mo" = "12m (dasa+ifna)",
                                                                "24mo" = "24m (dasa)"))
fwrite(pca_df, "results/manuscript/overall/pca_df.txt", sep = "\t", quote = F, row.names = F)

ggplot(pca_df, aes(PC1,PC2, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12)  + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs.pdf", width = 5, height = 5)

ggplot(pca_df, aes(PC1,NK, fill=log10(NK))) + geom_point(shape = 21, size = 3) + theme_bw(base_size = 12) + theme(legend.position = "right") + ggpubr::stat_cor() + scale_fill_viridis_c()
ggsave("results/manuscript/overall/pca_facs_nk.pdf", width = 5, height = 4)

ggplot(pca_df, aes(PC1,CD8.NAIVE, fill=log10(NK))) + geom_point(shape = 21, size = 3) + theme_bw(base_size = 12) #+ stat_ellipse()
ggplot(pca_df, aes(PC1,CD8_GrB, fill=log10(NK))) + geom_point(shape = 21, size = 3) + theme_bw(base_size = 12) #+ stat_ellipse()

ggplot(pca_df, aes(PC3,PC4, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12)+ theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_34.pdf", width = 5, height = 5)

ggplot(pca_df, aes(PC5,PC6, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_56.pdf", width = 5, height = 5)

ggplot(pca_df, aes(PC2,PC3, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_23.pdf", width = 5, height = 5)


ggplot(pca_df, aes(timepoint,PC1, fill=timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + labs(x="") + ggpubr::stat_compare_means()
ggplot(pca_df, aes(timepoint,PC2, fill=timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + labs(x="") + ggpubr::stat_compare_means()
ggplot(pca_df, aes(timepoint,PC3, fill=timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + labs(x="") + ggpubr::stat_compare_means()
ggplot(pca_df, aes(timepoint,PC4, fill=timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + labs(x="") + ggpubr::stat_compare_means()
ggplot(pca_df, aes(timepoint,PC5, fill=timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + labs(x="") + ggpubr::stat_compare_means()
ggplot(pca_df, aes(timepoint,PC6, fill=timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "none") + labs(x="") + ggpubr::stat_compare_means()

  

## Rotation
rot_df <- facs_full_pca$rotation %>% as.data.frame() %>% dplyr::select(PC1:PC6) %>% add_rownames(var = "feature") 

View(rot_df) 
melt(rot_df, id.vars = "feature")



## Cluster based on the top 6 PC:s
set.seed(232)
cluster_df <- pca_df %>% dplyr::select(PC1:PC6) %>% kmeans(centers = 5)
pca_df$cluster <- cluster_df$cluster %>% as.factor()

ggplot(pca_df, aes(PC1,PC2, fill=cluster)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) #+ stat_ellipse()
ggplot(pca_df, aes(PC2,PC3, fill=cluster)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) #+ stat_ellipse()
ggplot(pca_df, aes(PC3,PC4, fill=cluster)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) #+ stat_ellipse()

ggplot(pca_df, aes(cluster,NK,fill=cluster)) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) #+ stat_ellipse()
ggplot(pca_df, aes(cluster,CD8.TEMRA,fill=cluster)) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p")
ggplot(pca_df, aes(timepoint,CD8.TEMRA,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)


## Cluster only the base line
facs_0m_mtx <- facs_full_mtx[which(facs_full_clin$timepoint == "0mo"), ]
facs_0m_pca <- prcomp(facs_0m_mtx)

pca_0m_df <- facs_full_clin[which(facs_full_clin$timepoint == "0mo"), ] %>% bind_cols(as.data.frame(facs_0m_pca$x)[,1:10]) %>% bind_cols(facs_0m_mtx)

set.seed(1)
# cluster_0m_df <- pca_0m_df %>% dplyr::select(PC1:PC6) %>% kmeans(centers = 2)
cluster_0m_df <- facs_0m_mtx %>% kmeans(centers = 3)
facs_0m_mtx %>% dplyr::select(-NK) %>% scale() %>% pheatmap::pheatmap()


facs_0m_mtx %>% dplyr::select(-NK) %>% scale() %>% pheatmap::pheatmap()

table(cluster_0m_df$cluster)
pca_0m_df$cluster <- cluster_0m_df$cluster %>% as.factor()

ggplot(pca_0m_df, aes(PC1,PC2, fill=cluster)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) #+ stat_ellipse()

df <- pca_0m_df %>% dplyr::select(FM, cluster)

facs_clustered <- facs_full_clin %>% left_join(df) %>% bind_cols(facs_full_mtx) %>% filter(!is.na(cluster))
facs_clustered$timepoint <- factor(as.character(facs_clustered$timepoint), levels = c("0mo", "3mo", "12mo", "24mo"))
facs_clustered$timepoint <- plyr::revalue(facs_clustered$timepoint, replace = c("0mo"  = "0m",
                                                                                "3mo"  = "3m (dasa)",
                                                                                "12mo" = "12m (dasa+ifna)",
                                                                                "24mo" = "24m (dasa)"))
ggplot(pca_0m_df, aes(cluster,CD8.TEMRA,fill=cluster)) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p")

ggplot(facs_clustered, aes(timepoint,NK,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)
ggplot(facs_clustered, aes(timepoint,CD8.TEMRA,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)
ggplot(facs_clustered, aes(timepoint,CD8CD57,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)
ggplot(facs_clustered, aes(timepoint,CD8_GrB,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)

facs_clustered_infa <- facs_clustered %>% filter(timepoint %in% c("3m (dasa)", "12m (dasa+ifna)"))
var_df <- facs_clustered_infa %>% dplyr::select(FM, timepoint, cluster:CD8_TNF_IFN) %>% melt(id = c("timepoint", "FM", "cluster"))

p.df <- lapply(unique(var_df$variable), FUN = function(x){
  message(x)
  y <- var_df %>% filter(variable == x)
  df <- NULL
  i <- 1
  for(cluster_temp in unique(y$cluster)){
    # message(cluster_temp)
    z <- y %>% filter(cluster == cluster_temp)  
    if(length(unique(z$timepoint)) == 2){
      median.x = z %>% filter(timepoint == "3m (dasa)") %>% pull(value) %>% median(na.rm = T)
      median.y = z %>% filter(timepoint == "12m (dasa+ifna)") %>% pull(value) %>% median(na.rm = T)
      df[[i]] <- wilcox.test(value~timepoint, data = z) %>% broom::tidy() %>% mutate(variable = x, median.0m = median.x, median.3m = median.y, dir = ifelse(log2(median.y/median.x) > 0, "up", "down"), cluster = cluster_temp)
      i <- i + 1
    }
  }
  
  p <- df %>% rbindlist()
  return(p)
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) 
# fwrite(p.df, "results/manuscript/dasatinib/dasa_p_df.txt", sep = "\t", quote = F, row.names = F)

var_df_wo_0m <- var_df %>% filter(timepoint != "0mo")
krskl.p.df <- lapply(unique(var_df_wo_0m$variable), FUN = function(x){
  message(x)
  y <- var_df_wo_0m %>% filter(variable == x)
  df <- NULL
  i <- 1
  for(cluster_temp in unique(y$cluster)){
    # message(cluster_temp)
    z <- y %>% filter(cluster == cluster_temp)  
    df[[i]] <- kruskal.test(value~timepoint, data = z) %>% broom::tidy() %>% mutate(variable = x, cluster = cluster_temp)
    i <- i + 1
  }
  
  p <- df %>% rbindlist()
  return(p)
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) 
# fwrite(p.df, "results/manuscript/dasatinib/dasa_p_df.txt", sep = "\t", quote = F, row.names = F)


p.df %>% filter(cluster == 1) %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% View
krskl.p.df %>% filter(cluster == 1) %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% View

facs_clustered %>% filter(timepoint %in% c("3m (dasa)", "12m (dasa+ifna)")) %>% 
  ggplot(aes(timepoint,CD8.TEMRA,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)

facs_clustered %>% filter(timepoint %in% c("3m (dasa)", "12m (dasa+ifna)")) %>% 
  ggplot(aes(timepoint,NKCD27,fill=timepoint)) + facet_wrap(~cluster) + geom_boxplot() + scale_fill_manual(values = getPalette3(10)) + theme_bw(base_size = 12) + ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45)













######## PCA time; this time, remove the diagnosis 

facs_full      <- fread("results/manuscript/facs_df.txt")
facs_full_mtx  <- facs_full %>% filter(timepoint != "0mo") %>% dplyr::select(-c(timepoint:HRUH))
facs_full_clin <- facs_full %>% filter(timepoint != "0mo") %>% dplyr::select(c(timepoint:HRUH))

## Remove rows with over 75% NA
ind1 <- rowSums(is.na(facs_full_mtx)) != ncol(facs_full_mtx) 
ind2 <- apply(facs_full_mtx, 1, function(x)table(is.na(x))[1]) / 36 > 0.75

facs_full_mtx <- facs_full_mtx[intersect(which(ind1),which(ind2)), ]
facs_full_mtx <- facs_full_mtx %>% as.data.frame()

facs_full_clin <- facs_full_clin[intersect(which(ind1),which(ind2)), ]
facs_full_clin <- facs_full_clin %>% as.data.frame()

## Impute with col means
for(i in 1:ncol(facs_full_mtx)){
  facs_full_mtx[is.na(facs_full_mtx[,i]), i] <- mean(facs_full_mtx[,i], na.rm = TRUE)
}

facs_full_pca <- prcomp((facs_full_mtx))
pca_df <- facs_full_clin %>% bind_cols(as.data.frame(facs_full_pca$x)[,1:10]) %>% bind_cols(facs_full_mtx)

pca_df$timepoint <- factor(as.character(pca_df$timepoint), levels = c("3mo", "12mo", "24mo"))
pca_df$timepoint <- plyr::revalue(pca_df$timepoint, replace = c("3mo"  = "3m (dasa)",
                                                                "12mo" = "12m (dasa+ifna)",
                                                                "24mo" = "24m (dasa)"))
fwrite(pca_df, "results/manuscript/overall/pca_df_no_dg.txt", sep = "\t", quote = F, row.names = F)

ggplot(pca_df, aes(PC1,PC2, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_no_dg.pdf", width = 5, height = 5)

facs_full_pca$sdev / sum(facs_full_pca$sdev)


ggplot(pca_df, aes(PC3,PC4, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_34_no_dg.pdf",  width = 5, height = 5)

ggplot(pca_df, aes(PC5,PC6, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_56_no_dg.pdf",  width = 5, height = 5)

ggplot(pca_df, aes(PC2,PC3, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "")
ggsave("results/manuscript/overall/pca_facs_23_no_dg.pdf",  width = 5, height = 5)


pca_df$CD4.NAIVE
ggplot(pca_df, aes(PC3,CD4.NAIVE)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "") + ggpubr::stat_cor()
ggplot(pca_df, aes(PC3,CD8.NAIVE)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + theme(legend.position = "top") + labs(fill = "") + ggpubr::stat_cor()

