
dir.create("results/scvi/", showWarnings = F)
dir.create("results/scvi/input_files", showWarnings = F)
dir.create("results/scvi/results", showWarnings = F)

clonality_genes <- getClonalityGenes(cml_seurat)

## cml-samples to scVI; write .csv files
for(patient in unique(cml_seurat$orig.ident)){

  message(patient)
  # seurat_temp <- subset(cml_seurat, idents = patient)
  counts_temp <- cml_seurat@assays$RNA@counts[ ,cml_seurat$orig.ident == patient] %>% as.data.frame
  counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]

  data.table::fwrite(counts_temp, paste0("results/scvi/input_files/", patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)

}




## Get latent representation; UMAP it
cml_latent      <- fread("results/scvi/results/cml_oneshot_latent.csv")
cml_batch       <- fread("results/scvi/results/ccml_oneshot_indices.csv")
cml_latent_umap <- uwot::umap(cml_latent) # %>% select(UMAP1:UMAP2)
colnames(cml_latent_umap) <- c("UMAP1", "UMAP2")

cml_umap_df <- data.frame(orig.ident = cml_batch$V1, cml_latent_umap)
nClusters   <- cml_batch$V1 %>% unique() %>% length()

p <- ggplot() + 
  geom_point(data = cml_umap_df, aes(UMAP1, UMAP2, color = as.factor(orig.ident)), size = 0.3) + 
  add_guide + scale_color_manual(values = getPalette(nClusters)) + labs(color = "orig.ident")
ggsave(plot = p, "results/dimred/latent_umap.png", width = 10, height = 8)

## Put embeddings into Seurat object
cml_latent      <- as.matrix(cml_latent)
cml_latent_umap <- as.matrix(cml_latent_umap)

rownames(cml_latent)      <- colnames(cml_seurat)
rownames(cml_latent_umap) <- colnames(cml_seurat)

latent_dim_red              <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = cml_latent))
latent_umap_dim_red         <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = cml_latent_umap))

cml_seurat[['latent']]      <- latent_dim_red
cml_seurat[['latent_umap']] <- latent_umap_dim_red

# saveRDS(cml_seurat, "results/cml_seurat.rds")

DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T, group.by = "orig.ident") + theme(legend.position = "none")
DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T) + theme(legend.position = "none")
DimPlot(cml_seurat, reduction = "latent_umap", group.by = "cluster", label = T, repel = T) + theme(legend.position = "none")

cml_seurat <- RunUMAP(cml_seurat, reduction = "latent", dims = 1:30, n.neighbors = 10)
DimPlot(cml_seurat, reduction = "umap", label = T, repel = T, group.by = "orig.ident") + theme(legend.position = "none")
DimPlot(cml_seurat, reduction = "umap", label = T, repel = T) + theme(legend.position = "none")



## Do clustering
res                <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
cml_seurat         <- FindNeighbors(cml_seurat, reduction = "latent", dims = 1:30)
cml_seurat         <- FindClusters(object = cml_seurat, resolution = res)
clustering_columns <- grep("res", colnames(cml_seurat@meta.data), value = T)
saveRDS(cml_seurat, "results/cml_seurat_scvi.rds")


## Plot clustering results
p <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  nColors <- cml_seurat@meta.data[,clustering_column] %>% levels %>% length
  
  p[[i]] <- DimPlot(cml_seurat, reduction = "latent_umap", group.by = clustering_column, cols = getPalette(nColors), label = T) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = clustering_column)
  i <- i + 1
  
}

png("results/hvg/umap_per_res.png", width = 1024, height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()


q <- NULL; i <- 1
for(clustering_column in clustering_columns){
  q[[i]] <- cml_seurat@meta.data[,clustering_column] %>% levels %>% length; i <- i + 1
}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/hvg/scatter_res_nClusters.png", width = 4, height = 3)

## Choose 1
DimPlot(cml_seurat, reduction = "latent_umap", group.by = "RNA_snn_res.1", cols = getPalette(27), label = T)
Idents(cml_seurat) <- factor(cml_seurat$RNA_snn_res.1, levels = 0:23)# %>% reorderClusters()

DimPlot(cml_seurat, reduction = "latent_umap", cols = getPalette(27), label = T) + theme(legend.position = "none")
ggsave("results/hvg/latent_umap_clusters.png", width = 6, height = 5)


## Remove clusters 3 and 10, as they cluster only based on mitochondrial genes
idents.to.keep = Idents(cml_seurat)[!Idents(cml_seurat) %in% c("4 NK", "12 ?", "20 ?")] %>% unique()

cml_seurat <- subset(cml_seurat, idents = idents.to.keep)
cml_seurat <- getLatentUMAP(cml_seurat)
cml_seurat <- fixSeurat(cml_seurat)

saveRDS(cml_seurat, "results/cml_seurat_with_mito_clusters.rds")
# cml_seurat_mito <- readRDS("results/cml_seurat_with_mito_clusters.rds")


