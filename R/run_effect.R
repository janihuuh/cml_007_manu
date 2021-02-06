
dir.create("results/effectOnClusters/",  showWarnings = F)

## Run DEG on different timepoints by cluster to detemrine which clusters are affected
getDEGbyCluster <- function(seurat_object, cluster){
  
  message(paste0("===== ", cluster, " ====="))
  
  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= 50) return(NULL)
  
  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$timepoint
  
  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()
  
  n1 <- n_df %>% filter(Var1 == "dg") %>% pull(Freq) >= 50
  n2 <- n_df %>% filter(Var1 == "3mo") %>% pull(Freq) >= 50
  n3 <- n_df %>% filter(Var1 == "12mo") %>% pull(Freq) >= 50
  
  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE
  
  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL
  
  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "dg",  ident.2 = "3mo",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "2v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n1) cluster_markers_3v1 <- FindMarkers(object = seurat_cluster, ident.1 = "dg",  ident.2 = "12mo", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n2) cluster_markers_3v2 <- FindMarkers(object = seurat_cluster, ident.1 = "3mo", ident.2 = "12mo", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v2") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  
  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2) 
  
  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)
  
}


## Meta
DEG_cluster    <- lapply(levels(Idents(cml_seurat)), getDEGbyCluster, seurat_object = cml_seurat)
DEG_cluster_df <- DEG_cluster %>% rbindlist()
write.table(DEG_cluster_df, "results/effectOnClusters/deg_cluster.txt", sep = "\t", quote = F)

## Only R
r.cells.to.use   <- cml_seurat@meta.data %>% filter(effusion == "no effusion") %>% pull(barcode)
cml_seurat_r     <- subset(cml_seurat, cells = r.cells.to.use)
DEG_cluster_r    <- lapply(levels(Idents(cml_seurat)), getDEGbyCluster, seurat_object = cml_seurat_r)
DEG_cluster_r_df <- DEG_cluster_r %>% rbindlist()
write.table(DEG_cluster_r_df, "results/effectOnClusters/deg_r_cluster.txt", sep = "\t", quote = F)


## Only N
n.cells.to.use   <- cml_seurat@meta.data %>% filter(effusion == "effusion") %>% pull(barcode)
cml_seurat_n     <- subset(cml_seurat, cells = n.cells.to.use)
DEG_cluster_n    <- lapply(levels(Idents(cml_seurat)), getDEGbyCluster, seurat_object = cml_seurat_n)
DEG_cluster_n_df <- DEG_cluster_n %>% rbindlist()
write.table(DEG_cluster_n_df, "results/effectOnClusters/deg_n_cluster.txt", sep = "\t", quote = F)



#### #### #### #### #### #### #### #### 
#### #### ####  Analyze  #### #### #### 
#### #### #### #### #### #### #### #### 


DEG_cluster_df   <- read.delim("results/effectOnClusters/deg_cluster.txt")    %>% mutate(celltype = extractCoarsePhenotype(as.character(cluster)), cluster = getClusterPhenotypes(as.character(cluster)))
DEG_cluster_r_df <- read.delim("results/effectOnClusters/deg_r_cluster.txt")  %>% mutate(celltype = extractCoarsePhenotype(as.character(cluster)), cluster = getClusterPhenotypes(as.character(cluster)))
DEG_cluster_n_df <- read.delim("results/effectOnClusters/deg_n_cluster.txt")  %>% mutate(celltype = extractCoarsePhenotype(as.character(cluster)), cluster = getClusterPhenotypes(as.character(cluster)))


## Meta picture
DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_nDEG_2v1.pdf", width = 8, height = 6)

DEG_cluster_df %>% group_by(direction, timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1") %>% # %>% filter(cluster %in% 0:15) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(direction~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_nDEG_direction_2v1.pdf", width = 8, height = 6)

DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_nDEG_3v2.pdf", width = 8, height = 6)

DEG_cluster_df %>% group_by(direction, timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2") %>% # %>% filter(cluster %in% 0:15) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(direction~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_nDEG_direction_3v2.pdf", width = 8, height = 6)


## Meta picture responder
DEG_cluster_r_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_r_nDEG_2v1.pdf", width = 8, height = 6)

DEG_cluster_r_df %>% group_by(direction, timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1") %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(direction~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_r_nDEG_direction_2v1.pdf", width = 8, height = 6)

DEG_cluster_r_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_r_nDEG_3v2.pdf", width = 8, height = 6)

DEG_cluster_r_df %>% group_by(direction, timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2") %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(direction~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_r_nDEG_direction_3v2.pdf", width = 8, height = 6)


## Meta picture non-responder
DEG_cluster_n_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_n_nDEG_2v1.pdf", width = 8, height = 6)

DEG_cluster_n_df %>% group_by(direction, timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(direction~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_n_nDEG_direction_2v1.pdf", width = 8, height = 6)


DEG_cluster_n_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_n_nDEG_3v2.pdf", width = 8, height = 6)

DEG_cluster_n_df %>% group_by(direction, timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2")  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(reorder(cluster, n), n, fill = celltype, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(7)) +
  labs(x = "", y = "number of DEG as a function of time") + facet_wrap(direction~timepoint) + theme_bw() + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/effectOnClusters/bar_n_nDEG_direction_3v2.pdf", width = 8, height = 6)








## combine
DEG_cluster_r_summary_2v1 <- DEG_cluster_r_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "2v1")
DEG_cluster_n_summary_2v1 <- DEG_cluster_n_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "2v1")

deg_cluster_cummary_2v1 <- merge(DEG_cluster_r_summary_2v1[,-1], DEG_cluster_n_summary_2v1[,-1], by = "cluster", all = T) 
colnames(deg_cluster_cummary_2v1)[2:3] <- c("R", "N") 

deg_cluster_cummary_2v1 %>% melt(id = "cluster") %>% mutate(cluster = reorder(cluster, value)) %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(cluster,value,fill=variable)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + scale_fill_manual(values = c("salmon", "grey")) +
  theme_bw()
ggsave("results/effectOnClusters/bar_r_v_n_nDEG_2v1.pdf", width = 8, height = 6)

deg_cluster_cummary_2v1 %>% mutate(log2fc = log2(R/N))  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(cluster,log2fc,fill = cluster)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + scale_fill_manual(values = getPalette(21)) + 
  theme_bw() + theme(legend.position = "none") + 
  ylim(values = c(-5,5)) +
  geom_hline(yintercept = 1, linetype = "dotted") + geom_hline(yintercept = -1, linetype = "dotted")
ggsave("results/effectOnClusters/bar_r_v_n_nDEG_log2fc_2v1.pdf", width = 8, height = 6)




DEG_cluster_r_summary_3v2 <- DEG_cluster_r_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "3v2")
DEG_cluster_n_summary_3v2 <- DEG_cluster_n_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "3v2")

deg_cluster_cummary_3v2 <- merge(DEG_cluster_r_summary_3v2[,-1], DEG_cluster_n_summary_3v2[,-1], by = "cluster", all = T) 
colnames(deg_cluster_cummary_3v2)[2:3] <- c("R", "N") 

deg_cluster_cummary_3v2 %>% melt(id = "cluster") %>% mutate(cluster = reorder(cluster, value)) %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(cluster,value,fill=variable)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + scale_fill_manual(values = c("salmon", "grey")) +
  theme_bw()
ggsave("results/effectOnClusters/bar_r_v_n_nDEG_3v2.pdf", width = 8, height = 6)

deg_cluster_cummary_3v2 %>% mutate(log2fc = log2(R/N))  %>% filter(cluster %in% Idents(cml_seurat)) %>% 
  ggplot(aes(cluster,log2fc,fill = cluster)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + scale_fill_manual(values = getPalette(21)) + 
  theme_bw() + theme(legend.position = "none") + 
  ylim(values = c(-5,5)) + 
  geom_hline(yintercept = 1, linetype = "dotted") + geom_hline(yintercept = -1, linetype = "dotted")
ggsave("results/effectOnClusters/bar_r_v_n_nDEG_log2fc_3v2.pdf", width = 8, height = 6)

