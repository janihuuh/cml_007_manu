
getHVGscran <- function(seurat_object){

  ## @ params
  # object = seurat object

  # Using scran to pull the list of hvgs
  scran_object <- as.SingleCellExperiment(seurat_object)
  fit          <- scran::trendVar(scran_object, parametric = TRUE, use.spikes = FALSE)
  hvgs         <- scran::decomposeVar(scran_object, fit)
  return(hvgs)

}


## Plot HVGs
plotSeuratHVG <- function(object, top_n = 20){

  ## @ params
  # object = seurat object
  # top_n  = numeric, how many HVG:s to name

  var.features <- VariableFeatures(object)
  top_genes    <- var.features %>% head(n = top_n)

  hvf.info <- HVFInfo(object = object)  %>% add_rownames(var = "gene") %>% as.data.frame()
  hvf.info <- hvf.info %>% dplyr::mutate(var.status = ifelse(test = hvf.info$gene %in% var.features,  yes = "yes", no = "no"))

  ggplot() +
    geom_point(data = subset(hvf.info, var.status == "no"), aes(mean, variance.standardized), alpha = 1, color = "grey", size = 0.8) +
    geom_point(data = subset(hvf.info, var.status == "yes"), aes(mean, variance.standardized), alpha = 0.8, color = "dodgerblue", size = 0.8) +
    ggrepel::geom_text_repel(data = subset(hvf.info, gene %in% top_genes), aes(mean, variance.standardized, label = gene), fontface = "italic") +

    scale_x_log10() + scale_y_log10() +
    theme_bw() + labs(x = "Average expression", y = "Standardized variance")

}


plotScranHVG <- function(scran_hvg, top_n = 20){

  ## @ params
  # scran_hvg = df, scran_hvg
  # top_n  = numeric, how many HVG:s to name
  
  top_genes    <- scran_hvg %>% filter(FDR < 0.05) %>% top_n(wt = total, n = top_n) %>% pull(gene)

  ggplot(scran_hvg, aes(mean, total, fill = ifelse(FDR < 0.05, "yes", "no"))) + geom_point(pch = 21, size = 2) +
    # geom_smooth(fill = "darkred") +
    ggrepel::geom_text_repel(data = subset(scran_hvg, gene %in% top_genes), aes(mean, total, label = gene), fontface = "italic") +
    scale_fill_manual(values = c("grey", "darkred")) +

    labs(fill = "p.adj < 0.05", y = "total variance", x = "mean expression")


}
