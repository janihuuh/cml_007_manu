
runSeurat    <- function(seurat_object, cells.to.use, doPCA = F){
  
  clonotype_seurat_meta           <- filter(seurat_object@meta.data, barcode %in% cells.to.use)
  rownames(clonotype_seurat_meta) <- clonotype_seurat_meta$barcode
  
  clonotype_seurat                <- subset(seurat_object, cells = clonotype_seurat_meta$barcode)
  clonotype_seurat                <- clonotype_seurat %>% NormalizeData(verbose = F) %>% ScaleData(verbose = F)
  clonotype_seurat@meta.data      <- clonotype_seurat_meta[colnames(clonotype_seurat), ]
  clonotype_seurat@meta.data$orig.clusters <- Idents(clonotype_seurat)
  
  ## For variable genes, we use the genes that are already calculated in the meta seurat function
  message(paste("analysing", length(cells.to.use), "cells"))
  
  ## Compute up to 50 PCs (if under 50 cells then obviosuly under 50)
  
  if(doPCA == T){
    
    message(paste("using", length(VariableFeatures(seurat_object)),  "HVGs"))
    pcs.compute = 50
    
    if(ncol(clonotype_seurat) <= 50){
      pcs.compute = ncol(clonotype_seurat@data) - 1
    }
    
    clonotype_seurat <- RunPCA(object = clonotype_seurat, pc.genes = VariableFeatures(seurat_object), pcs.compute = pcs.compute, verbose = F)
    nPCs             <- sum(clonotype_seurat[["pca"]]@stdev > 1.5)
    
    # Cluster the cells with graph-based clustering after PCA
    # clonotype_seurat    <- FindNeighbors(object = clonotype_seurat, reduction.type = "pca", dims.use = 1:nPCs, verbose = F)
    # clonotype_seurat    <- FindClusters(object = clonotype_seurat, resolution = 1, verbose = F)
    
  }
  
  return(clonotype_seurat)
  
}

runSlingshot <- function(sce_object, cluster_with_earliset_timepoint, reducedDim){
  
  require(slingshot); require(SummarizedExperiment)
  sce <- slingshot(data = sce_object, clusterLabels = 'orig.clusters', reducedDim = reducedDim, start.clus = cluster_with_earliset_timepoint)
  
}

plotSlingshot <- function(slingshot_object, reducedDim){
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  
  cololors <- slingshot_object$timepoint
  cololors[cololors == 1] <- "salmon"
  cololors[cololors == 2] <- "darkolivegreen"
  cololors[cololors == 3] <- "dodgerblue"
  
  curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")
  
  colors2 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(9)
  
  sling_curve <- SlingshotDataSet(slingshot_object)
  
  par(mfrow = c(1,2))
  plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
  for(i in 1:length(sling_curve@curves)){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }
  
  plot(reducedDims(slingshot_object)[[reducedDim]], col = colors[cut(slingshot_object$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
  # lines(SlingshotDataSet(slingshot_object), lwd = 2, col = "black")
  for(i in 1:length(sling_curve@curves)){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }
  
  par(mfrow = c(1,1))
  
  
}

calculateSlingshot <- function(slingshot_object, pseudotime){
  
  require(gam)
  require(SingleCellExperiment)
  
  # slingshot_object = clonotype_sling
  # for time, only look at the 1,000 most variable genes
  
  # t     <- slingshot_object$slingPseudotime_1 # + 1 %>% log
  t     <- pseudotime
  
  
  Y     <- log1p(assays(slingshot_object)$count)
  var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
  Y     <- Y[var1K,]
  
  # fit a GAM with a loess term for pseudotime
  gam.pval <- apply(Y,1,function(z){
    
    d   <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    p   <- summary(tmp)[4][[1]][1,5]
    
    return(p)
    
  })
  
  gam.pval <- data.frame(gam.pval) %>%
    add_rownames("gene")
  
  return(gam.pval)
  
}

plotPseudotimeGenes <- function(clonotype_sce, genes_of_choice){
  
  df <- data.frame(pseudotime = clonotype_sce$Pseudotime,
                   timepoint  = clonotype_sce$timepoint,
                   t(as.matrix(assays(clonotype_sce)$logcounts[genes_of_choice, ])))
  
  df %>% melt(id = c("pseudotime", "timepoint")) %>%
    ggplot(aes(log(pseudotime),value,color=timepoint)) +
    geom_point(shape = 21, fill = "grey", color = 'black') +
    stat_smooth(method = "loess", color = "darkred", fill = "lightgrey", alpha = 0.3) +
    facet_wrap(~variable, ncol = 5)
  
}

analyseClonotypes <- function(seurat_object, clonotypename, output_dir){
  
  
  ########
  ## This is the workhorse function
  ########
  
  ## Takes one clonotype, makes a Seurat object of it and plots it
  ## Then counts pseudotime with slingshot and the most DE-genes along pseudotime and visualizes them
  ## Does not output anything, but saves files into output_dir
  
  # @ param:
  # input :
  # seurat_object = seurat_object
  # clonotypename = char, that is in seurat_object@new_clonotypes_id
  # output_dir    = path, where to plot the files
  
  
  dir.create(output_dir, showWarnings = F)
  # clonotypename <- substr(clonotypename, 10, nchar(clonotypename))
  
  message(paste("=============", clonotypename, "============="))
  cells.to.use <- seurat_object@meta.data %>% filter(new_clonotypes_id == clonotypename) %>% select(barcode)
  cells.to.use <- as.character(cells.to.use$barcode)
  
  ## Visualise
  require(gridExtra)
  
  # i) seurat
  clonotype_seurat <- runSeurat(seurat_object = seurat_object, cells.to.use = cells.to.use, doPCA = F)
  clonotype_sce    <- as.SingleCellExperiment(clonotype_seurat, to = "sce")
  message("Seurat done")
  
  orig.clusters <- unique(Idents(clonotype_seurat))[order(unique(Idents(clonotype_seurat)))]
  orig.clusters <- substr(orig.clusters, 1,2) %>% as.numeric()
  
  ## PCA
  pca_df <- cbind(clonotype_seurat@meta.data,
                  # clonotype_seurat[["pca"]]@cell.embeddings[,1:6],
                  # clonotype_seurat[["umap"]]@cell.embeddings,
                  clonotype_seurat[["latent_umap"]]@cell.embeddings,
                  cluster = Idents(clonotype_seurat)) %>%
    dplyr::rename(latent_umap_1 = latent_1, latent_umap_2 = latent_2)
  
  
  # a      <- ggplot(pca_df, aes(PC_1, PC_2, color = timepoint)) + geom_point() + scale_color_brewer(palette = "Set3")
  # b      <- ggplot(pca_df, aes(PC_1, PC_2, color = cluster)) + geom_point() + scale_color_manual(values = getPalette(nClusters)[orig.clusters+1])
  
  # c      <- ggplot(pca_df, aes(latent_umap_1, latent_umap_2, color = timepoint)) + geom_point() + scale_color_brewer(palette = "Paired")
  # d      <- ggplot(pca_df, aes(latent_umap_1, latent_umap_2, color = cluster)) + geom_point() + scale_color_manual(values = getPalette(nClusters)[orig.clusters+1])
  #
  # pdf(paste0(output_dir, clonotypename, "_pca.pdf"), width = 10, height = 4)
  # grid.arrange(c,d,ncol=2)
  # dev.off()
  #
  
  
  ## Plot clonotype evolution
  clonotype_seurat@meta.data %>% tidyr::complete(timepoint, orig.clusters, fill = list(z = 0)) %>%
    group_by(timepoint, orig.clusters) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) %>%
    
    ggplot(aes(timepoint,prop, fill = orig.clusters, group=orig.clusters)) + geom_area() + scale_fill_manual(values = getPalette(nClusters)[orig.clusters+1]) + theme_bw()
  ggsave(paste0(output_dir, clonotypename, "_evolution1.png"), width = 6, height = 4)
  
  clonotype_seurat@meta.data %>% tidyr::complete(timepoint, orig.clusters, fill = list(z = 0)) %>%
    group_by(timepoint, orig.clusters) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) %>%
    
    ggplot(aes(timepoint,n, fill = orig.clusters, group = orig.clusters)) + geom_area() + scale_fill_manual(values = getPalette(nClusters)[orig.clusters+1]) + theme_bw()
  ggsave(paste0(output_dir, clonotypename, "_evolution2.png"), width = 6, height = 4)
  
  message("Visualising seurat done")
  
  
  ## Select the cluster with the most cells at the earliest timepoint to be the start to find meaningful trajectory
  top_cl <- clonotype_seurat@meta.data %>%
    group_by(orig.clusters, timepoint)
  
  earliest <- unique(top_cl$timepoint)[1]
  top_cl <- top_cl %>% dplyr::filter(timepoint == earliest) %>%
    summarise(n = n()) %>%
    arrange(desc(n))
  
  cluster_with_earliset_timepoint <- top_cl[1,1] %>% as.numeric()
  
  
  # ii) pseudotime
  ## Calcualate pseudotime with slingshot. It can result in multiple different trajectories; let's take them all in consieration
  clonotype_sling     <- runSlingshot(clonotype_sce, cluster_with_earliset_timepoint = cluster_with_earliset_timepoint, reducedDim = "LATENT_UMAP")
  
  
  ## If more than one trajectory appears, plot them all
  pseudotimes <- grep("slingPseudotime*", colnames(colData(clonotype_sling)), value = T)
  
  for(pseudotime in pseudotimes){
    
    message(paste("Calculating", pseudotime, "..."))
    
    clonotype_sce$Pseudotime = clonotype_sling[[pseudotime]]
    
    slingshot_genes     <- calculateSlingshot(slingshot_object = clonotype_sling, pseudotime = clonotype_sce$Pseudotime )
    slingshot_genes_top <- slingshot_genes %>% dplyr::filter(gam.pval < 0.05) %>% top_n(20, wt = -gam.pval)
    
    
    if(nrow(slingshot_genes_top) != 0){
      
      p <- plotPseudotimeGenes(clonotype_sce = clonotype_sce, genes_of_choice = slingshot_genes_top$gene)
      ggsave(plot = p, paste0(output_dir, clonotypename, "_", pseudotime, "_sling_genes.png"), width = 12, height = 9)
      
      pdf(paste0(output_dir, clonotypename, "_", pseudotime, "_pseudo_sling.pdf"), width = 12, height = 6)
      plotSlingshot(clonotype_sling, reducedDim = "LATENT_UMAP")
      dev.off()
      
      p <- FeaturePlot(clonotype_seurat, features = slingshot_genes_top$gene, cols = c("lightgrey", "salmon"), reduction =  "latent_umap", ncol = 5)
      ggsave(plot = p, paste0(output_dir, clonotypename, "_", pseudotime, "_sling_features.png"), width = 14, height = 10)
      
    }
    
    
    write.table(arrange(slingshot_genes, gam.pval), paste0(output_dir, clonotypename, "_", pseudotime, "_sling_genes.txt"), sep = "\t", row.names = F)
    
    
  }
  
  
}

getDEGbyClonotype <- function(seurat_object, clonotypename){
  
  ## If under 30 cells to begin with
  if(is.na(clonotypename)) return(NULL)
  if(table(seurat_object$new_clonotypes_id) %>% as.data.frame() %>% filter(Var1 == clonotypename) %>% pull(Freq) <= 30) return(NULL)
  
  ## Subet to only cluster
  message(paste("=============", clonotypename, "============="))
  cells.to.use <- seurat_object@meta.data %>% filter(new_clonotypes_id == clonotypename) %>% pull(barcode)
  
  clonotype_seurat_meta           <- filter(seurat_object@meta.data, barcode %in% cells.to.use)
  rownames(clonotype_seurat_meta) <- clonotype_seurat_meta$barcode
  clonotype_seurat                <- subset(seurat_object, cells = clonotype_seurat_meta$barcode)
  clonotype_seurat                <- clonotype_seurat %>% NormalizeData(verbose = F) %>% ScaleData(verbose = F)
  
  Idents(clonotype_seurat) <- clonotype_seurat$timepoint
  
  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(clonotype_seurat)) %>% as.data.frame()
  
  n1 <- n_df %>% filter(Var1 == 1) %>% pull(Freq) >= 5
  n2 <- n_df %>% filter(Var1 == 2) %>% pull(Freq) >= 5
  n3 <- n_df %>% filter(Var1 == 3) %>% pull(Freq) >= 5
  
  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE
  
  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL
  
  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = clonotype_seurat, ident.1 = "1", ident.2 = "2", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "2v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n1) cluster_markers_3v1 <- FindMarkers(object = clonotype_seurat, ident.1 = "1", ident.2 = "3", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n2) cluster_markers_3v2 <- FindMarkers(object = clonotype_seurat, ident.1 = "2", ident.2 = "3", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v2") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  
  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2)
  
  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(new_clonotypes_id = clonotypename, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)
  
}

plotLatentUmapClonotype <- function(seurat_object, clonotype, output_dir, name = ""){
  
  message(clonotype)
  
  if(nchar(name) == 0) name = clonotype
  
  ## Init visualisation object
  viz_df <- cbind(seurat_object@meta.data,
                  # seurat_object[["pca"]]@cell.embeddings[,1:6],
                  # seurat_object[["umap"]]@cell.embeddings,
                  seurat_object[["latent_umap"]]@cell.embeddings,
                  cluster = Idents(seurat_object)) %>% 
    mutate(cluster = factor(cluster, levels = levels(cluster)[order(levels(cluster))])) %>% 
    dplyr::rename(umap_1 = latent_1 , umap_2 = latent_2)
  
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
