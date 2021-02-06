
## Init output options
cml_seurat <- cml_seurat

seurat_file_path <- paste0("results/cml_seurat_qc.rds")
output_dir       <- paste0("results/hvg/")
dir.create(output_dir, showWarnings = F)



## No need to modify
#### =====================



####################

scale_factor = 10e3     ## Scaling factor
n_hvgs       = 3000     ## How many HVGs to calculate (minus the clonality and unwanted variation genes)
n_pc_stdev   = 2        ## The minimun stdev value for PC to be accepted

####################






## Normalize
cml_seurat <- NormalizeData(cml_seurat, normalization.method = "LogNormalize", scale.factor = scale_factor)

## ==== Find HVGs with Seurat
cml_seurat <- FindVariableFeatures(cml_seurat, selection.method = "vst", nfeatures = n_hvgs)
seurat_hvg    <- VariableFeatures(cml_seurat)

## Remove clonality and unwanted genes
clonality_genes <- getClonalityGenes(cml_seurat)
unwanted_genes  <- getUnwantedGenes(cml_seurat)

seurat_hvg      <- seurat_hvg[!seurat_hvg %in% clonality_genes]
seurat_hvg      <- seurat_hvg[!seurat_hvg %in% unwanted_genes]
write.table(seurat_hvg, paste0(output_dir, "seurat_hvg.txt"), sep = "\t", quote = F, row.names = F)

## Visualize the HVG function
VariableFeatures(cml_seurat) <- seurat_hvg

plotSeuratHVG(cml_seurat, top_n = 20)
ggsave(paste0(output_dir, "seurat_hvg.pdf"), width = 8, height = 6)



## ==== Calculate HVGs with scran 
scran_hvg <- getHVGscran(cml_seurat)
scran_hvg <- scran_hvg %>% as.data.frame() %>% add_rownames(var = "gene")
scran_hvg <- scran_hvg[!scran_hvg$gene %in% clonality_genes, ]
scran_hvg <- scran_hvg[!scran_hvg$gene %in% unwanted_genes, ]

plotScranHVG(scran_hvg, top_n = 20)
ggsave(paste0(output_dir, "scran_hvg.pdf"), width = 8, height = 6)


scran_hvg_sigf <- scran_hvg %>% arrange(FDR) %>% filter(total > 0.5 & FDR < 0.05)

write.table(scran_hvg, paste0(output_dir, "scran_hvg.txt"), sep = "\t", quote = F, row.names = F)
write.table(scran_hvg_sigf, paste0(output_dir, "scran_hvg_sigf.txt"), sep = "\t", quote = F, row.names = F)



############ Choose which hvgs to use ############
# intersect(seurat_hvg, scran_hvg_sigf$gene)
# table(scran_hvg_sigf$gene %in% seurat_hvg)

# hvg = seurat_hvg
hvg = scran_hvg_sigf$gene
VariableFeatures(cml_seurat) <- hvg

###############################################



## Scale data based on the HVGs
cml_seurat <- ScaleData(cml_seurat, features = hvg)

## Perform  PCA
cml_seurat <- RunPCA(cml_seurat, features = hvg, npcs = 50)

## Choose nPC based on stdev
nPCs <- sum(cml_seurat[["pca"]]@stdev > n_pc_stdev)
message(paste("nPCs:", nPCs, "\n"))
write.table(nPCs, paste0(output_dir, "nPCs.txt"), sep = "\t", quote = F, row.names = F)

## RunUMAP-function does not work
cml_seurat <- RunUMAP(cml_seurat, dims = 1:nPCs, learning.rate = 1)

# Meanwhile try something hacky-ish
# umap_df <- cml_seurat[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
# umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
# cml_seurat[['umap']] <- umap_df

DimPlot(cml_seurat, split.by = "orig.ident")
ggsave("results/hvg/umap_split.png", width = 24, height = 4)

DimPlot(cml_seurat, reduction = "umap", repel = T, label = T) + theme(legend.position = "none")
ggsave("results/hvg/umap.png", width = 6, height = 4)

DimPlot(cml_seurat, reduction = "umap", group.by = "orig.ident", repel = T, label = T) + theme(legend.position = "none")




## Clustering
res                <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
cml_seurat         <- FindNeighbors(cml_seurat, reduction = "pca", dims = 1:nPCs)
cml_seurat         <- FindClusters(object = cml_seurat, resolution = res)
clustering_columns <- grep("res", colnames(cml_seurat@meta.data), value = T)
# saveRDS(cml_seurat, "results/cml_seurat_scvi.rds")


## Plot clustering results
p <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  nColors <- cml_seurat@meta.data[,clustering_column] %>% levels %>% length
  
  p[[i]] <- DimPlot(cml_seurat, reduction = "umap", group.by = clustering_column, cols = getPalette(nColors), label = T) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = clustering_column)
  i <- i + 1
  
}

png("results/hvg/umap_per_res.png", width = 1024, height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()


q <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  q[[i]] <- cml_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
  
}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/hvg/scatter_res_nClusters.png", width = 4, height = 3)

saveRDS(cml_seurat, file =  "results/cml_seurat_qc.rds")


## Choose 0.6
DimPlot(cml_seurat, reduction = "umap", group.by = "RNA_snn_res.0.6", cols = getPalette(19), label = T)
Idents(cml_seurat) <- factor(cml_seurat$RNA_snn_res.0.6, levels = 0:18)

DimPlot(cml_seurat, label = T) + theme(legend.position = "none")
ggsave("results/hvg/umap_clusters.png", width = 6, height = 5)


cml_seurat$clusters <- getClusterPhenotypes(Idents(cml_seurat))
DimPlot(cml_seurat, label = T, group.by = "clusters", cols = getPalette(19)) + theme(legend.position = "none")
ggsave("results/hvg/umap_cluster_names.png", width = 9, height = 7)


