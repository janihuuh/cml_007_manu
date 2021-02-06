
cml_seurat  <- readRDS("results/cml_seurat_scvi.rds")
# saveRDS(cml_seurat, "results/cml_seurat_scvi.rds")

# tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt")

# cml_seurat <- fixSeurat(cml_seurat)
cml_seurat$patient   <- extractName(cml_seurat$orig.ident)
cml_seurat$timepoint <- extractTimepoint(cml_seurat$orig.ident) %>% factor(levels = c("dg", "3mo", "12mo"))
cml_seurat$effusion  <- ifelse(cml_seurat$patient %in% c(716,720), "effusion", "no effusion")
Idents(cml_seurat)   <- factor(cml_seurat$RNA_snn_res.0.5, levels = 0:20) %>% getClusterPhenotypes()
nClusters            <- length(levels(Idents(cml_seurat)))

DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T, cols = getPalette(nClusters)) + theme_void() + theme(legend.position = "none") 

cml_seurat2  <- readRDS("results/cml_seurat_with_mito_clusters.rds")

DimPlot(cml_seurat, reduction = "latent_umap", label = T, repel = T, cols = getPalette(23)) + theme_void() + theme(legend.position = "none") 
DimPlot(cml_seurat2, reduction = "latent_umap", label = T, repel = T, cols = getPalette(23)) + theme_void() + theme(legend.position = "none") 



Idents(cml_seurat) <- factor(cml_seurat$RNA_snn_res.0.5, levels = 0:20) %>% getClusterPhenotypes()
length(unique(cml_seurat$RNA_snn_res.0.7))
Idents(cml_seurat) <- factor(cml_seurat$RNA_snn_res.0.7, levels = 0:22) %>% getClusterPhenotypes()

DimPlot(cml_seurat, reduction = "umap", label = T, repel = T, cols = getPalette(23)) + theme_void() + theme(legend.position = "none") 


