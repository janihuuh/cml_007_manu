

initCellphonedb <- function(df, seurat_object, name = "", sample_n = 50, folder = "results/cellphonedb/input_files/", seed = 123){
  
  set.seed(seed)
  
  ## If one wants to use subsampling; remove clusters with less than sample_n cells
  if(!is.null(sample_n)){
    
    clusters.to.rm  <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n < sample_n) %>% pull(cluster)
    message(paste("remove", clusters.to.rm))
    clusters.to.use <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n >= sample_n) %>% pull(cluster)
    df              <- df %>% filter(cluster %in% clusters.to.use) %>% group_by(cluster)
    df              <- df %>% sample_n(size = sample_n) %>% ungroup()
  }
  
  ## Remove empty clusters
  df <- df %>% filter(cluster != "")
  
  counts <- seurat_object@assays$RNA@counts[,df$barcode]
  counts[Matrix::rowSums(counts) != 0, ] %>% as.data.frame %>% add_rownames(var = "Gene") %>%
    dplyr::select(Gene, everything()) %>%
    fwrite(paste0(folder, name, "_counts.txt"), sep = "\t", quote = F, row.names = F)
  
  df %>% dplyr::select(barcode, cluster) %>% dplyr::rename(Cell = barcode, cell_type = cluster) %>%
    fwrite(paste0(folder, name, "_meta.txt"), sep = "\t", quote = F, row.names = F)
  
  return(NULL)
  
}

getPairAmounts <- function(df){
  
  ## Detect columns with pairs
  cols.to.use    <- grep("\\|", colnames(df), value = T)
  # cols.to.use    <- c("interacting_pair", grep("\\|", colnames(df), value = T))
  df             <- df %>% dplyr::select(cols.to.use)
  pairs.complete <- lapply(cols.to.use, function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c")
  pairs.unique   <- unique(pairs.complete)
  # pairs.unique   <- unique(pairs.complete)[-1]
  
  i <- 1
  df_list <- NULL
  
  for(pair in pairs.unique){
    
    # message(pair)
    # pair <- c("interacting_pair", pair)
    df_pair <- df %>% dplyr::select(colnames(df)[pairs.complete %in% pair])
    
    df_list[[i]] <- apply(df_pair[,], 2, function(x) table(!is.na(x))[2]) %>% # bind_cols(data.frame(interacting_pair = df_pair$interacting_pair)) %>%
      as.data.frame() %>% add_rownames(var = "pair") %>% dplyr::rename(value = ".") %>% mutate(cluster1 = lapply(pair, function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c"),
                                                                                               cluster2 = lapply(pair, function(x) strsplit(x, "[|]")[[1]][2]) %>% do.call(what = "c"))
    
    i <- i + 1
    
  }
  
  df1 <- df_list %>% rbindlist()
  
  
}

filterRedundantPairs <- function(df1){
  
  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
  
  
  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))
  df1 <- df1 %>% filter(pair %in% pairs$pairs)
  return(df1)
  
}

getSigfPairs <- function(df){
  
  # df = sigf_means_r_pre
  
  ## Detect columns with pairs
  # cols.to.use    <- grep("\\|", colnames(df), value = T)
  cols.to.use    <- c("interacting_pair", grep("\\|", colnames(df), value = T))
  df             <- df %>% dplyr::select(cols.to.use)
  pairs.complete <- lapply(cols.to.use, function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c")
  # pairs.unique   <- unique(pairs.complete)
  pairs.unique   <- unique(pairs.complete)[-1]
  
  i <- 1
  df_list <- NULL
  
  for(pair in pairs.unique){
    
    message(pair)
    
    pair <- c("interacting_pair", pair)
    df_pair <- df %>% dplyr::select(colnames(df)[pairs.complete %in% pair])
    
    df_list[[i]]  <- melt(df_pair, id = "interacting_pair") %>% filter(!is.na(value)) %>% mutate(cluster1 = lapply(as.character(variable), function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c"),
                                                                                                 cluster2 = lapply(as.character(variable), function(x) strsplit(x, "[|]")[[1]][2]) %>% do.call(what = "c"))
    
    i <- i + 1
    
  }
  
  df1 <- df_list %>% rbindlist()
  
  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
  
  
  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))
  df1 <- df1 %>% filter(variable %in% pairs$pairs)
  return(df1)
  
}

getNewInteractions <- function(df1,df2){
  
  df1 <- df1 %>% mutate(temp = paste(interacting_pair, variable))
  df2 <- df2 %>% mutate(temp = paste(interacting_pair, variable))
  
  new_df <- df2[!df2$temp %in% df1$temp, ]
  
}

plotPairs <- function(df1, df2){
  
  # sigf_r <- sigf_r %>% dplyr::rename(group1 = median_pre, group2 = median_post, p = p.value) %>% mutate(y.position = 30)
  # sigf_r <- sigf_r %>% select(group1, group2, p, y.position)
  
  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)
  
  nClustersTemp <- df$cluster1.x %>% unique() %>% length()
  df <- df %>% dplyr::select(pair, pre, post, cluster1.x, cluster2.x) %>% melt(id = c("pair", "cluster1.x", "cluster2.x"))
  df <- df %>% mutate(phenotype = extractCoarsePhenotype(cluster2.x))
  
  ggplot(df, aes(variable,value,group=pair,fill=phenotype)) +
    geom_point(shape = 21, size = 3) + geom_path() +
    # ggpubr::stat_pvalue_manual(sigf_r) +
    # ggsignif::geom_signif(comparisons = list(c("pre", "post")), test.args = c(paired = T)) +
    # ggrepel::geom_label_repel(data = subset(df, variable == "post"), aes(variable,value,label=cluster2.x,color=cluster2.x), fill = NA) +
    facet_wrap(~cluster1.x, scales = "free") + facets_nice + labs(fill = "", x = "", y = "# of sigf receptor-ligand interactions") +
    # scale_fill_manual(values = getPalette(nClustersTemp))
    scale_fill_manual(values = brewer.pal(9, "Pastel1"))
  
}

plotNewPairs <- function(df1, df2){
  
  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)
  
  nClustersTemp <- df$cluster1.x %>% unique() %>% length()
  df <- df %>% dplyr::select(pair, pre, post, cluster1.x, cluster2.x)
  
  df$post <- c(df$post - df$pre) %>% as.numeric
  df$pre  <- 0
  
  nClusters <- df$cluster2.x %>% unique %>% length
  
  
  df <- df %>% melt(id = c("pair", "cluster1.x", "cluster2.x")) %>% mutate(phenotype = extractCoarsePhenotype(cluster2.x))
  
  df %>% filter(variable == "post") %>%
    ggplot(aes(cluster2.x, value, fill = cluster2.x)) +
    geom_bar(stat = "identity") + facet_wrap(~cluster1.x) + scale_fill_manual(values = getPalette(nClusters)) +
    theme(axis.text.x = element_blank())
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}

getNewPairs <- function(df1, df2){
  
  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)
  
  nClustersTemp <- df$cluster1.x %>% unique() %>% length()
  df <- df %>% dplyr::select(pair, pre, post, cluster1.x, cluster2.x)
  
  df$post <- c(df$post - df$pre) %>% as.numeric
  df$pre  <- 0
  
  nClusters <- df$cluster2.x %>% unique %>% length
  
  df <- df %>% melt(id = c("pair", "cluster1.x", "cluster2.x")) %>% mutate(phenotype = extractCoarsePhenotype(cluster2.x))
  
  
}

testSigfPairs <- function(df1, df2, paired = T){
  
  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)
  
  df_list <- NULL
  i <- 1
  
  for(cluster in unique(df$cluster1.x)){
    
    # cluster = unique(df$cluster1.x)[1]
    df_temp <- df %>% filter(cluster1.x == cluster)
    
    df_list[[i]] <- wilcox.test(df_temp$pre, df_temp$post, paired = paired) %>% broom::tidy() %>%
      mutate(cluster = cluster, median_pre = median(df_temp$pre, na.rm = T), median_post = median(df_temp$post, na.rm = T),
             direction = ifelse(median_post > median_pre, "up", "down"))
    i <- i + 1
    
  }
  
  df_list %>% rbindlist()
  
}


expand.grid.unique <- function(x, y, incl.eq = TRUE){
  g <- function(i){
    z <- setdiff(y, x[seq_len(i - incl.eq)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

f <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}


plotDiffHeatmap <- function(df1, df2){
  
  require(igraph)
  require(ComplexHeatmap)
  
  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
  
  f <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  
  
  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))
  
  df1 <- df1 %>% filter(pair %in% pairs$pairs)
  df2 <- df2 %>% filter(pair %in% pairs$pairs)
  
  df_mat1 <- cast(df1, cluster1 ~ cluster2, sum)
  rownames(df_mat1) <- df_mat1$cluster1
  
  df_mat2 <- cast(df2, cluster1 ~ cluster2, sum)
  rownames(df_mat2) <- df_mat2$cluster1
  
  to.use1 <- rownames(df_mat1) %in% rownames(df_mat2)
  to.use2 <- rownames(df_mat2) %in% rownames(df_mat1)
  
  df_mat2 <- df_mat2 %>% as.data.frame()
  df_mat2[to.use2,c(F,to.use2)]
  
  
  df_diff <- log2(df_mat2[to.use2,c(F,to.use2)] / df_mat1[to.use1,c(F,to.use1)])
  df_diff <- apply(df_diff, 1, function(x){ x[!is.finite(x)] <- NA; return(x)} )
  
  
  anno_df <- data.frame(phenotype = extractCoarsePhenotype(rownames(df_mat1)))
  rownames(anno_df) <- rownames(df_mat1)
  
  nAnno <- anno_df$phenotype %>% unique %>% length()
  col = getPalette(nAnno)
  col = brewer.pal(nAnno, "Pastel1")
  
  pheatmap::pheatmap(f(t(df_diff)), 
                     cluster_rows = F, cluster_cols = F, 
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(20),
                     annotation_row = anno_df, 
                     annotation_col = anno_df, 
                     
                     # annotation_colors = anno_col,
                     na_col = "lightgrey")
  
  
  
}


plotDensity <- function(seurat_object){
  
  nProjects <- length(unique(seurat_object$project))+2
  patient_df <- seurat_object@meta.data %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
  
  df            <- seurat_object@meta.data %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>% filter(!is.na(project))
  median_df     <- df %>% group_by(project, cluster) %>% summarise(median = median(prop)) %>% group_by(cluster) %>% top_n(n = 5, wt = median) %>% arrange(desc(median)) %>% filter(!is.na(project))
  umap_means_df <- seurat_object@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(seurat_object@meta.data) %>% 
    group_by(cluster) %>% summarise(umap1 = median(latentumap_1), umap2 = median(latentumap_2)) %>% 
    left_join(median_df) %>% filter(!is.na(project))
  
  seurat_object@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(seurat_object@meta.data) %>%
    ggplot(aes(latentumap_1, latentumap_2, color = project)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
    facet_wrap(~project, ncol = 4) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
    ggrepel::geom_text_repel(data = umap_means_df, aes(umap1,umap2,label=cluster), color = "black", size = 3.5, direction = "x",angle = 0,vjust = 0, segment.size = 0.2) +
    scale_color_manual(values = getPalette3(nProjects)) + guides(color=FALSE) + theme(legend.position = "none")
  
}


## Calculate slingshot curve

getSlingClonotype <- function(seurat_object, clone){
  
  cells.to.keep       <- seurat_object@meta.data %>% filter(new_clonotypes_id == clone) %>% pull(barcode)
  project_seurat      <- subset(seurat_object, cells = cells.to.keep) %>% as.SingleCellExperiment()
  sce                 <- slingshot(data = project_seurat, clusterLabels = 'cluster', reducedDim = "LATENT_UMAP", start.clus = "8 CD8 naive TCF7+") %>% SlingshotDataSet
  return(sce)
  
}



readGliphFile <- function(filename){
  
  gliph_df   <- fread(filename, fill = T, stringsAsFactors = F)
  
  gliph_df$vb_score <- as.numeric(gliph_df$vb_score)
  gliph_df$Fisher_score <- as.numeric(gliph_df$Fisher_score)
  gliph_df$final_score <- as.numeric(gliph_df$final_score)
  gliph_df$number_unique_cdr3 <- as.numeric(gliph_df$number_unique_cdr3)
  
  gliph_df$number_subject <- as.numeric(gliph_df$number_subject)
  gliph_df$length_score <- as.numeric(gliph_df$length_score)
  gliph_df$cluster_size_score <- as.numeric(gliph_df$cluster_size_score)
  
  ## Remove everything after last numeric index
  max_index  <- max(as.numeric(gliph_df$index), na.rm = T)
  max_indexs <- which(gliph_df$index == max_index)
  gliph_df   <- gliph_df[1:max_indexs[length(max_indexs)], ]
  return(gliph_df)
  
}

vdjToGliph <- function(vdj_df){
  
  ## Write gliph files to the vdj files
  # @ param
  # input: df from vdj
  
  df <- vdj_df %>% dplyr::select(cdr3aa, v, j, name, count) %>%
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = count)
  
  return(df)
  
}

require(clusterProfiler)
require(org.Hs.eg.db)

if(me == "hru"){
  
  hallmark   <- read.gmt("/Users/hru/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")
  
  
}


if(me == "janihuuh"){
  
  hallmark   <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")
  
}




## == No need to modify

plotHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  
  if(table(enrich@result$p.adjust < 0.05) %>% length() > 1){
    heatplot(enrich)
  }
  
  else(NULL)
  
}

getHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  if(!is.null(enrich)){
    
    enrich <- do.call(rbind, enrich@result) %>% t %>% as.data.frame()
    enrich[,c(5:7, 9)] <- sapply(enrich[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
    return(enrich)
    
  }
  
  
}



plotDiffHeatmap <- function(df1,df2){
  
  require(igraph)
  require(ComplexHeatmap)
  
  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
  
  f <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  
  
  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))
  
  df1 <- df1 %>% filter(pair %in% pairs$pairs)
  df2 <- df2 %>% filter(pair %in% pairs$pairs)
  
  df_mat1 <- cast(df1, cluster1 ~ cluster2, sum)
  rownames(df_mat1) <- df_mat1$cluster1
  
  df_mat2 <- cast(df2, cluster1 ~ cluster2, sum)
  rownames(df_mat2) <- df_mat2$cluster1
  
  to.use1 <- rownames(df_mat1) %in% rownames(df_mat2)
  to.use2 <- rownames(df_mat2) %in% rownames(df_mat1)
  
  df_mat2 <- df_mat2 %>% as.data.frame()
  df_mat2[to.use2,c(F,to.use2)]
  
  df_diff <- log2(df_mat2[to.use2,c(F,to.use2)] / df_mat1[to.use1,c(F,to.use1)])
  df_diff <- apply(df_diff, 1, function(x){ x[!is.finite(x)] <- NA; return(x)} )
  
  pheatmap::pheatmap(f(t(df_diff)),
                     cluster_rows = T,
                     cluster_cols = T,
                     
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, 0.1))),
                     breaks = seq(-2, 2, 0.1),
                     # color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(-3, 3, 0.1))),
                     # breaks = seq(-3, 3, 0.1),
                     
                     display_numbers = T,
                     number_format ="%.1f",
                     
                     # annotation_row = anno_df,
                     # annotation_col = anno_df,
                     # annotation_colors = anno_col,
                     na_col = "lightgrey")
  
}


hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
blueprint <- SingleR::BlueprintEncodeData()

getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){
  
  # hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
  # blueprint <- SingleR::BlueprintEncodeData()
  
  ## @ params
  ## cluster = possible cluster vec, if not provided, tries to find in meta.data$cluster
  ## method = if "cluster", then performs preds based on clusters, not cells
  ## sample = to subsample or not 
  
  if(!is.null(sample)){
    
    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])
    
  }
  
  sce       <- as.SingleCellExperiment(seurat_object)
  
  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)
    
    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.hca$first.labels
      return(seurat_object)  
    }
    
    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }
    
  }
  
  
  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}



getClusterPhenotypes <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "0"  = "0 NK CD56dim",
    "1"  = "1 CD8 EM GZMK+ CXCR3+" ,
    "2"  = "2 CD8 cytotoxic GZMH+" ,
    "3"  = "3 NK CD56dim HAVCR2+" ,
    "4"  = "4 CD4 CM/naive" ,
    "5"  = "5 CD8 EMRA ZNF683+" ,
    "6"  = "6 Monocyte CD16- HAVCR2+" ,
    "7"  = "7 CD4 Th1-like CTLA4+ CXCR3+" ,
    "8"  = "8 CD8 naive TCF7+" ,
    "9"  = "9 CD8 EMRA IFNG+" ,
    "10" = "10 B-cell immature" ,
    "11" = "11 B-cell memory" ,
    "12" = "12 CD8 MAIT-like CXCR3+" ,
    "13" = "13 NK CD56bright" ,
    "14" = "14 CD8/B-cell doublet" ,
    "15" = "15 Monocyte CD16+" ,
    "16" = "16 B-cell immature" ,
    "17" = "17 cycling" ,
    "18" = "18 T-cell HAVCR2+" ,
    "19" = "19 pDC" ))
  
  return(clusters)
  
}

getClusterPhenotypesXXX <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "0"  = "0 ",
    "1"  = "1 " ,
    "2"  = "2 " ,
    "3"  = "3 " ,
    "4"  = "4 " ,
    "5"  = "5 " ,
    "6"  = "6 " ,
    "7"  = "7 " ,
    "8"  = "8 " ,
    "9"  = "9 " ,
    "10" = "10 " ,
    "11" = "11 " ,
    "12" = "12 " ,
    "13" = "13 " ,
    "14" = "14 " ,
    "15" = "15 " ,
    "16" = "16 " ,
    "17" = "17 " ,
    "18" = "18 " ,
    "19" = "19 " ,
    "20" = "20 " ,
    "21" = "21 " ,
    "22" = "22 ",
    "23" = "23 ",
    "24" = "24 ",
    "25" = "25 "))
  
  return(clusters)
  
}

plotClustering <- function(seurat_object){
  
  res <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]
  
  q <- NULL; i <- 1
  
  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }
  
  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
  
}

putLatentsSeurat <- function(seurat_object, latent){
  
  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)
  
  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)
  
  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)
  
  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))
  
  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}

getScviInput <- function(seurat_object, folder){
  
  dir.create(folder, showWarnings = F)
  genes_to_keep <- list()
  i <- 1
  Idents(seurat_object) <- seurat_object$orig.ident
  
  for(patient in unique(seurat_object$orig.ident)){
    
    message(patient)
    seurat_temp <- subset(seurat_object, idents = patient)
    counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame
    
    genes_to_keep[[i]] <- rownames(counts_temp)
    i <- i + 1
    
    counts_temp <- seurat_object@assays$RNA@counts[ ,seurat_object$orig.ident == patient] %>% as.data.frame
    counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
    data.table::fwrite(counts_temp, paste0(folder, patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)
    
  }
}



getQC <- function(seurat_object){
  
  ###################
  
  min_mito     <- 0
  max_mito     <- 10
  
  min_ribo     <- 10
  max_ribo     <- 50
  
  min_features <- 600
  max_features <- 3000
  
  min_counts   <- 1.75e3
  max_counts   <- 10e3
  
  min_pct50    <- 25
  max_pct50    <- 60
  
  cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                    "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                    "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                    "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                    "TUBB", "TYMS", "UBE2C")
  
  ###################
  
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^RP", col.name = "percent.ribo")
  # seurat_object  <- PercentageFeatureSet(seurat_object, features = cycle.genes, col.name = "percent.cycle")
  seurat_object@meta.data$barcode <- colnames(seurat_object)
  
  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.mt") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.ribo") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nFeature_RNA") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nCount_RNA") + geom_hline(yintercept = 30e3) + geom_hline(yintercept = 5e2)
  
  
  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()
  
  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)
  
  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))
  
  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))
  
  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)
  
}


getHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  enrich <- do.call(rbind, enrich@result) %>% t %>% as.data.frame()
  enrich[,c(5:7, 9)] <- sapply(enrich[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(enrich)
  
}

getDEGbyCluster <- function(seurat_object, cluster){
  
  
  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= 50) return(NULL)
  message(paste0("===== ", cluster, " ====="))
  
  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$timepoint
  
  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()
  
  n1 <- n_df %>% filter(Var1 == "dg") %>% pull(Freq) >= 10
  n2 <- n_df %>% filter(Var1 == "3mo") %>% pull(Freq) >= 10
  n3 <- n_df %>% filter(Var1 == "12mo") %>% pull(Freq) >= 10
  
  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE
  
  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL
  
  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "dg", ident.2 = "3mo", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3movdg") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n1) cluster_markers_3v1 <- FindMarkers(object = seurat_cluster, ident.1 = "dg", ident.2 = "12mo", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "12movdg") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n2) cluster_markers_3v2 <- FindMarkers(object = seurat_cluster, ident.1 = "3mo", ident.2 = "12mo", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "12mov3mo") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  
  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2)
  
  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)
  
}


reorderClusters <- function(cluster_vec){
  
  ## Get clusters in order
  clusters <- cluster_vec %>% unique()
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)
  
}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){

  str1 <- substr(str1, 1, nchar(str1) - 27)
  extractFileName(str1)
}

extractTimepoint <- function(strs){
  
  strs2 <-NULL
  i <- 1
  for(str1 in strs){
    strs2[[i]] <- strsplit(str1, "[_]")[[1]][2]
    i <- i + 1
  }
  
  return(strs2)

}


plotQcViolin <- function(viz_df, var_to_plot, grouping, min, max){

  ## Plot univariate violin plots with filter thresholds

  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable

  viz_df_temp <- viz_df %>% dplyr::select(var_to_plot)

  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table

  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +

    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +

    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +

    labs(x = "", title = var_to_plot) + theme(legend.position = "none")

}




getClusterPhenotypesOld <- function(clusters){

  # input : vector of clusters 

  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 CD8+ effector",

                                                    "1"  = "1 NK CD56dim",
                                                    "2"  = "2 NK CD56dim",
                                                    "3"  = "3 CD4+ T-cells",
                                                    "4"  = "4 Monocytes CD16-",
                                                    "5"  = "5 CD4+ Tem",
                                                    "6"  = "6 B-cell naive",
                                                    "7"  = "7 CD8+ effector/exhausted",
                                                    "8"  = "8 Poor quality T cells",
                                                    "9"  = "9 Monocytes CD16+ small",
                                                    "10" = "10 B-cell class-switched memory",

                                                    "11" = "11 CD8+ effector memory",
                                                    "12" = "12 NK cells",
                                                    "13" = "13 B-cell",
                                                    "14" = "14 NK CD56bright",
                                                    "15" = "15 NK cells",
                                                    "16" = "16 Monocytes CD16+",
                                                    "17" = "17 Monocytes CD16+",
                                                    "18" = "18 Cycling cells",
                                                    "19" = "19 Erythrocytes",
                                                    "20" = "20 Monocyte pDC"
                                                    
                                                    
                                                    ))


  return(clusters)

}


getClusterPhenotypesOldButNotThatOld <- function(clusters){
  
  # input : vector of clusters 
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 NK CD56dim",
                                                    
                                                    "1"  = "1 CD8 effector/exhausted",
                                                    "2"  = "2 CD8 EM",
                                                    "3"  = "3 NK CD56dim IFNg cytotoxic",
                                                    "4"  = "4 NK CD56dim IFNg exhausted",
                                                    "5"  = "5 CD4 naive",
                                                    "6"  = "6 CD8 EMRA",
                                                    "7"  = "7 CD4 Th1-like",
                                                    "8"  = "8 Monocytes CD16-",
                                                    "9"  = "9 B-cell naive",
                                                    "10" = "10 Cytotoxic lymphocyte",
                                                    
                                                    "11" = "11 B-cell memory",
                                                    "12" = "12 CD8 naive/CM",
                                                    "13" = "13 Cytotoxic lymphocyte",
                                                    "14" = "14 CD4 EMRA",
                                                    "15" = "15 CD4 CM",
                                                    "16" = "16 T cell low quality",
                                                    "17" = "17 NK CD56bright",
                                                    "18" = "18 Monocytes CD16+",
                                                    "19" = "19 cycling",
                                                    "20" = "20 Monocytes low quality",
                                                    "21" = "21 Monocytes low quality",
                                                    "22" = "22 pDC"))
  
  
  return(clusters)
  
}





fixSeurat <- function(seurat_object){
  
  ## Fix meta data if it brokes
  
  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)
  
  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]
  
  rownames(meta.data) <- meta.data$barcode
  
  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)
  
  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg
  
  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)
  
}





extractCoarsePhenotype <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }
  
  return(p)
  
}
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


getLatentUMAP <- function(seurat_object){
  
  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df
  
  return(seurat_object)
  
}






calculateFoldchange <- function(viz_df){
  
  df_temp <- viz_df %>%
    tidyr::complete(timepoint, cluster, fill = list(z = 0)) %>%
    
    group_by(timepoint, cluster) %>%
    summarise(n = n()) %>% mutate(prop = n / sum(n))
  
  df_temp1 <- df_temp %>% filter(timepoint == "dg")
  df_temp2 <- df_temp %>% filter(timepoint == "3mo")
  df_temp3 <- df_temp %>% filter(timepoint == "12mo")
  
  df_temp <- data.frame(cluster = df_temp1$cluster,
                        "0" = 0,
                        "1" = log2(df_temp2$prop / df_temp1$prop),
                        "3" = log2(df_temp3$prop / df_temp2$prop)) %>% melt(id = "cluster")
  
  return(df_temp)
  
}



plotBarFoldchange <- function(df){
  
  df <- df %>%
    mutate(direction = "stable") %>%
    mutate(direction = ifelse(value > 1, "increase", direction)) %>%
    mutate(direction = ifelse(value < -1, "decrease", direction)) %>%
    mutate(direction = factor(direction, levels = c("increase", "stable", "decrease")))
  
  limits <- max(abs(df$value))
  
  ggplot(df, aes(x = cluster, value, fill = direction)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_hline(yintercept = -1, linetype = "dotted") + theme_bw() +
    labs(y = "log2(fc)", x = "") +
    ylim(-limits, limits)
  
  
}


doLatentClustering <- function(seurat_object){
  
  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = 1:30)
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)
  
}




getSingler <- function(seurat_object, folder){
  
  ### Run singler on public
  ## Try SingleR
  dir.create(folder, showWarnings = F)
  
  for(patient in unique(seurat_object$orig.ident)){
    
    require(SingleR)
    message(patient)
    
    ## Cells to select
    cells.to.keep <- seurat_object@meta.data[seurat_object$orig.ident == patient, ] %>% pull(barcode)
    public_temp      = subset(seurat_object, cells = cells.to.keep)
    
    annot = data.frame(public_temp@meta.data)
    rownames(annot) <- colnames(public_temp)
    
    singler = SingleR::CreateSinglerObject(public_temp@assays$RNA@counts,
                                           project.name = patient,
                                           min.genes = 0,
                                           technology = "10X",
                                           species = "Human",
                                           do.signatures = F,
                                           fine.tune = F,
                                           clusters = Idents(public_temp))
    
    saveRDS(singler, file = paste0(folder, patient, '.rds'))
    
  }
  
  
  ## Combine the data sets
  singler.objects.file <- list.files(folder,  full.names = T)
  singler.objects.file <- grep(".rds", singler.objects.file, value = T)
  singler.objects      <- lapply(singler.objects.file, FUN = function(x) {message(x); readRDS(x)})
  
  singler = SingleR.Combine(singler.objects,
                            order    = colnames(seurat_object),
                            clusters = Idents(seurat_object),
                            xy       = seurat_object@reductions$umap@cell.embeddings)
  
  metadata <- data.frame(barcode        = names(singler$meta.data$orig.ident),
                         hpca_pred      = singler$singler[[1]]$SingleR.single$labels,
                         blueprint_pred = singler$singler[[2]]$SingleR.single$labels) %>% add_rownames(var = "barcode")
  
  public_metadata <- merge(seurat_object@meta.data, metadata, by = "barcode")
  
  
  ## Add into Seurat
  public_metadata                      <- public_metadata[match(colnames(seurat_object), public_metadata$barcode), ]
  seurat_object$singler_hpca_pred      <- public_metadata$hpca_pred
  seurat_object$singler_blueprint_pred <- public_metadata$blueprint_pred
  
  return(seurat_object)
}


## Run seurat
preprocessSeurat <- function(orig_object, cells.to.use){
  
  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)
  
  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta
  
  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)
  
  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]
  
  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))
  
  ## Scale data
  object <- ScaleData(object, features = hvg)
  
  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))
  
  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)
  
  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df
  
  return(object)
  
}


plotSlingshot <- function(slingshot_object, reducedDim){
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  
  cololors                <- slingshot_object$orig.clusters %>% extractClusterNumber() %>% as.numeric()
  cololors[cololors == 0] <- getPalette5(8)[1]
  cololors[cololors == 1] <- getPalette5(8)[2]
  cololors[cololors == 2] <- getPalette5(8)[3]
  cololors[cololors == 3] <- getPalette5(8)[4]
  cololors[cololors == 4] <- getPalette5(8)[5]
  cololors[cololors == 5] <- getPalette5(8)[6]
  cololors[cololors == 6] <- getPalette5(8)[7]
  cololors[cololors == 7] <- getPalette5(8)[8]
  
  
  curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")
  
  colors2 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(9)
  
  sling_curve <- SlingshotDataSet(slingshot_object)
  
  par(mfrow = c(1,2))
  plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
  for(i in 1:2){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }
  
  plot(reducedDims(slingshot_object)[[reducedDim]], col = colors[cut(slingshot_object$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
  # lines(SlingshotDataSet(slingshot_object), lwd = 2, col = "black")
  for(i in 1:2){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }
  
  par(mfrow = c(1,1))
  
  
}


runSlingshot <- function(sce_object, cluster_with_earliset_timepoint, reducedDim){
  
  require(slingshot); require(SummarizedExperiment)
  sce <- slingshot(data = sce_object, clusterLabels = 'orig.clusters', reducedDim = reducedDim, start.clus = cluster_with_earliset_timepoint)
  
}

getNewClusters <- function(clusters){
  clusters %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml()
}

extractCoarsePhenotype <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }
  
  return(p)
  
}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

reorderClusters <- function(cluster_vec){
  
  ## Get clusters in order
  clusters <- cluster_vec %>% unique()
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)
  
}


extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}







## For GSEA
forGSEA <- function(de_df){
  
  ## Create .rnk files for GSEA
  rnk  <- de_df %>% arrange(desc(avg_logFC)) %>% select(gene, avg_logFC)
  colnames(rnk) <- c('Name', 'metric')
  return(rnk)
  
}



extractName = function(str1){
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  sub(".*\\/", "", str1)
}



getLatentUMAP <- function(seurat_object){
  
  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df
  
  return(seurat_object)
  
}

fixSeurat <- function(seurat_object){
  
  ## Fix meta data if it brokes
  
  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)
  
  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]
  
  rownames(meta.data) <- meta.data$barcode
  
  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)
  
  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg
  
  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)
  
}


getLatentClustering <- function(seurat_object){
  
  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)
  
}
