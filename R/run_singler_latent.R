
## Try SingleR with latent 
require(SingleR)
dir.create("results/singler/", showWarnings = F)

for(patient in unique(cml_seurat$orig.ident)[2:12]){
  
  message(patient)
  
  ## Cells to select
  cells.to.keep <- cml_seurat@meta.data[cml_seurat$orig.ident == patient, ] %>% pull(barcode)
  # cells.to.keep = cml_seurat@meta.data %>% filter(orig.ident == patient) %>% pull(barcode)
  cml_temp      = subset(cml_seurat, cells = cells.to.keep)
  
  annot = data.frame(cml_temp@meta.data)
  rownames(annot) <- colnames(cml_temp)
  
  singler = SingleR::CreateSinglerObject(cml_temp@assays$RNA@counts, 
                                         # annot = annot,
                                         project.name = patient, 
                                         min.genes = 0,
                                         technology = "10X", 
                                         species = "Human", 
                                         # citation = "Huuhtanen et al., 2020",
                                         do.signatures = F, 
                                         fine.tune = F,  
                                         clusters = Idents(cml_temp))
  
  singler$meta.data$orig.ident = cml_temp@meta.data$orig.ident                   # the original identities, if not supplied in 'annot'
  singler$meta.data$xy         = cml_temp@reductions$latent_umap@cell.embeddings # the latent umap coordinates
  singler$meta.data$clusters   = Idents(cml_temp)                                # the Seurat clusters (if 'clusters' not provided)
  
  median_df <- data.frame(singler$meta.data$xy, label = singler$singler[[1]]$SingleR.single$labels) %>% group_by(label) %>% summarise(UMAP_1 = median(latentumap_1), UMAP_2 = median(latentumap_2))
  
  p <- data.frame(singler$meta.data$xy, label = singler$singler[[1]]$SingleR.single$labels) %>% 
    ggplot(aes(latentumap_1, latentumap_2, color = label)) + geom_point(size = 0.3) + add_guide + scale_color_manual(values = getPalette(45)) +
    ggrepel::geom_label_repel(data = median_df, aes(UMAP_1, UMAP_2, label = label)) # + theme(legend.position = "none")
  ggsave(plot = p, paste0("results/singler/latent_umap_predictions_", patient, "_hpca.png"), width = 16, height = 8)
  
  median_df <- data.frame(singler$meta.data$xy, label = singler$singler[[2]]$SingleR.single$labels) %>% group_by(label) %>% summarise(UMAP_1 = median(latentumap_1), UMAP_2 = median(latentumap_2))
  
  p <- data.frame(singler$meta.data$xy, label = singler$singler[[2]]$SingleR.single$labels) %>% 
    ggplot(aes(latentumap_1, latentumap_2, color = label)) + geom_point(size = 0.3) + add_guide + scale_color_manual(values = getPalette(45)) +
    ggrepel::geom_label_repel(data = median_df, aes(UMAP_1, UMAP_2, label = label)) # + theme(legend.position = "none")
  ggsave(plot = p, paste0("results/singler/latent_umap_predictions_", patient, "_blueprint.png"), width = 16, height = 8)
  
  saveRDS(singler, file = paste0('results/singler/cml_temp_latent_', patient, '.rds'))
  
}




## Combine the data sets
singler.objects.file <- list.files('results/singler/', pattern = 'cml_temp', full.names = T)
singler.objects.file <- grep(".rds", singler.objects.file, value = T)
singler.objects      <- lapply(singler.objects.file, FUN = function(x) {message(x); readRDS(x)})

singler = SingleR.Combine(singler.objects, 
                          order    = colnames(cml_seurat), 
                          clusters = Idents(cml_seurat),
                          xy       = cml_seurat@reductions$latent_umap@cell.embeddings)
saveRDS(singler, 'results/singler/cml_singler.rds')
singler <- readRDS('results/singler/cml_singler.rds')


## HPCA based predictions
latent_umap_df <- cml_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% add_rownames(var = "barcode")
label_df       <- data.frame(label = singler$singler[[1]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

p <- latent_umap_df %>% left_join(label_df) %>% 
  ggplot(aes(latent_umap_1, latent_umap_2, color = label)) + geom_point(size = 0.3) + add_guide + scale_color_manual(values = getPalette(63)) + facet_wrap(~label, ncol = 6)
ggsave(plot = p, "results/singler/latent_umap_predictions_hpca.png", width = 30, height = 40)

## Blueprint-encode based predictions
latent_umap_df <- cml_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% add_rownames(var = "barcode")
cluster_df     <- Idents(cml_seurat) %>% as.data.frame()  %>% add_rownames(var = "barcode")
colnames(cluster_df) <- c("barcode", "cluster")
label_df       <- data.frame(label = singler$singler[[2]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

p <- latent_umap_df %>% left_join(label_df) %>% 
  ggplot(aes(latent_umap_1, latent_umap_2, color = label)) + geom_point(size = 0.3) + add_guide + scale_color_manual(values = getPalette(30)) + facet_wrap(~label, ncol = 6)
ggsave(plot = p, "results/singler/latent_umap_predictions_blueprint.png", width = 30, height = 20)


## Based on cluster-level predictions
singler$singler[[1]]$SingleR.clusters$labels
singler$singler[[2]]$SingleR.clusters$labels






## HPCA
latent_umap_df <- cml_seurat@reductions$umap@cell.embeddings %>% as.data.frame %>% add_rownames(var = "barcode")
cluster_df     <- Idents(cml_seurat) %>% as.data.frame()  %>% add_rownames(var = "barcode")
colnames(cluster_df) <- c("barcode", "cluster")
label_df       <- data.frame(label = singler$singler[[1]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

df <- latent_umap_df %>% left_join(label_df)  %>% left_join(cluster_df)
tot_df <- table(df$label, df$cluster) %>% melt %>% group_by(Var2) %>% summarise(n = sum(value))

table(df$label, df$cluster) %>% melt %>% group_by(Var2) %>% 
  top_n(n = 5, wt = value) %>% left_join(tot_df) %>% mutate(prop = round(value/n, 3)) %>% filter(prop > 0.25) %>% 
  
  ggplot(aes(Var2,prop,fill=Var1,label=Var1)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = getPalette(18)) +
  theme_bw() + coord_flip() + geom_hline(yintercept = 0.5, linetype = "dotted") + 
  ggrepel::geom_text_repel() + 
  facet_wrap(~Var2)
ggsave("results/singler/bar_cluster_predictions_hpca.png", width = 16, height = 12)




## Blueprint
latent_umap_df <- cml_seurat@reductions$umap@cell.embeddings %>% as.data.frame %>% add_rownames(var = "barcode")
cluster_df     <- Idents(cml_seurat) %>% as.data.frame()  %>% add_rownames(var = "barcode")
colnames(cluster_df) <- c("barcode", "cluster")
label_df       <- data.frame(label = singler$singler[[2]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

df <- latent_umap_df %>% left_join(label_df)  %>% left_join(cluster_df)
tot_df <- table(df$label, df$cluster) %>% melt %>% group_by(Var2) %>% summarise(n = sum(value))

table(df$label, df$cluster) %>% melt %>% group_by(Var2) %>% 
  top_n(n = 5, wt = value) %>% left_join(tot_df) %>% mutate(prop = round(value/n, 3)) %>% filter(prop > 0.25) %>% 
  
  ggplot(aes(Var2,prop,fill=Var1,label=Var1)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = getPalette(18)) +
  theme_bw() + coord_flip() + geom_hline(yintercept = 0.5, linetype = "dotted") + 
  ggrepel::geom_text_repel() + 
  facet_wrap(~Var2)
ggsave("results/singler/bar_cluster_predictions_blueprint.png", width = 16, height = 12)



cml_seurat$hpca      <- Idents(cml_seurat) %>% getClusterPhenotypesHPCA()
cml_seurat$blueprint <- Idents(cml_seurat) %>% getClusterPhenotypesBlueprint()

p <- DimPlot(cml_seurat, reduction = "umap", group.by = "hpca", label = T, repel = T, cols = getPalette(15)) + theme_void()
ggsave(plot = p, "results/singler/latent_umap_predictions_cluster_hpca.png", width = 12, height = 8)

p <- DimPlot(cml_seurat, reduction = "umap", group.by = "blueprint", label = T, repel = T, cols = getPalette(21)) + theme_void()
ggsave(plot = p, "results/singler/latent_umap_predictions_cluster_blueprint.png", width = 12, height = 8)

