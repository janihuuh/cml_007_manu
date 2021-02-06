
## Run cellphone db. Make subselection of 1000 cells from responders at time 0 and time 1
dir.create("results/cellphonedb/", showWarnings = F)
dir.create("results/cellphonedb/input_files", showWarnings = F)
dir.create("results/cellphonedb/analysis", showWarnings = F)

# sample_n = 100
# sample_n = 15

################################ Init cellphone input files

clusters.to.rm <- c("18 T-cell HAVCR2+", "17 cycling", "14 CD8/B-cell doublet")

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(timepoint == "dg") %>% group_by(cluster) %>% summarise(n = n()) %>% arrange((n))
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(timepoint == "3mo") %>% group_by(cluster) %>% summarise(n = n()) %>% arrange((n))
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(timepoint == "12mo") %>% group_by(cluster) %>% summarise(n = n()) %>% arrange((n))
  
### Total
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(timepoint == "dg") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "tot_0m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(timepoint == "3mo") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "tot_1m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(timepoint == "12mo") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "tot_3m", sample_n = sample_n)



### Responders
sample_n = 100

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(effusion == "no effusion", timepoint == "dg") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "r_0m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(effusion == "no effusion", timepoint == "3mo") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "r_1m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(effusion == "no effusion", timepoint == "12mo") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "r_3m", sample_n = sample_n)




### Non-responders
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(effusion == "effusion", timepoint == "dg") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "n_0m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(effusion == "effusion", timepoint == "3mo") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "n_1m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(effusion == "effusion", timepoint == "12mo") %>% filter(!cluster %in% clusters.to.rm) %>% 
  initCellphonedb(seurat_object = cml_seurat, name = "n_3m", sample_n = sample_n)





################################ Analyze

sigf_means_tot_0m <- fread("results/cellphonedb/out/tot_0m/significant_means.txt")
sigf_means_tot_1m <- fread("results/cellphonedb/out/tot_1m/significant_means.txt")
sigf_means_tot_3m <- fread("results/cellphonedb/out/tot_3m/significant_means.txt")

sigf_means_r_0m <- fread("results/cellphonedb/out/r_0m/significant_means.txt")
sigf_means_r_1m <- fread("results/cellphonedb/out/r_1m/significant_means.txt")
sigf_means_r_3m <- fread("results/cellphonedb/out/r_3m/significant_means.txt")

sigf_means_n_0m <- fread("results/cellphonedb/out/n_0m/significant_means.txt")
sigf_means_n_1m <- fread("results/cellphonedb/out/n_1m/significant_means.txt")
sigf_means_n_3m <- fread("results/cellphonedb/out/n_3m/significant_means.txt")


## How many interactions are there?
sigf_means_tot_0m_pairs <- getPairAmounts(sigf_means_tot_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()
sigf_means_tot_1m_pairs <- getPairAmounts(sigf_means_tot_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()
sigf_means_tot_3m_pairs <- getPairAmounts(sigf_means_tot_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()

sigf_means_r_0m_pairs <- getPairAmounts(sigf_means_r_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()
sigf_means_r_1m_pairs <- getPairAmounts(sigf_means_r_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()
sigf_means_r_3m_pairs <- getPairAmounts(sigf_means_r_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()

sigf_means_n_0m_pairs <- getPairAmounts(sigf_means_n_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()
sigf_means_n_1m_pairs <- getPairAmounts(sigf_means_n_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()
sigf_means_n_3m_pairs <- getPairAmounts(sigf_means_n_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2)) %>% filterRedundantPairs()



## Which phenotypes gets/loses sigf amount of interactions?
sigf_r    <- testSigfPairs(sigf_means_r_0m_pairs, sigf_means_r_1m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_r, "results/cellphonedb/analysis/de_clusters_r_2v1.txt", sep = "\t", quote = F, row.names = F)

sigf_r    <- testSigfPairs(sigf_means_r_1m_pairs, sigf_means_r_3m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_r, "results/cellphonedb/analysis/de_clusters_r_3v2.txt", sep = "\t", quote = F, row.names = F)

sigf_n    <- testSigfPairs(sigf_means_n_0m_pairs, sigf_means_n_1m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_n, "results/cellphonedb/analysis/de_clusters_n_2v1.txt", sep = "\t", quote = F, row.names = F)

sigf_n    <- testSigfPairs(sigf_means_n_1m_pairs, sigf_means_n_3m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_n, "results/cellphonedb/analysis/de_clusters_n_3v2.txt", sep = "\t", quote = F, row.names = F)


sigf_diff_0m <- testSigfPairs(sigf_means_r_0m_pairs, sigf_means_n_0m_pairs, paired = F) %>% filter(p.value < 0.05)
write.table(sigf_diff_0m, "results/cellphonedb/analysis/de_clusters_r_vs_n_0m.txt", sep = "\t", quote = F, row.names = F)

sigf_diff_1m <- testSigfPairs(sigf_means_r_1m_pairs, sigf_means_n_1m_pairs, paired = F) %>% filter(p.value < 0.05)
write.table(sigf_diff_1m, "results/cellphonedb/analysis/de_clusters_r_vs_n_1m.txt", sep = "\t", quote = F, row.names = F)

sigf_diff_3m <- testSigfPairs(sigf_means_r_3m_pairs, sigf_means_n_3m_pairs, paired = F) %>% filter(p.value < 0.05)
write.table(sigf_diff_3m, "results/cellphonedb/analysis/de_clusters_r_vs_n_3m.txt", sep = "\t", quote = F, row.names = F)


## Visualize Which phenotypes gets/loses sigf amount of interactions?
p <- plotPairs(sigf_means_tot_0m_pairs, sigf_means_tot_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_tot_2v1.png", width = 12, height = 10)

p <- plotPairs(sigf_means_r_0m_pairs, sigf_means_r_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_r_2v1.png", width = 12, height = 10)

p <- plotPairs(sigf_means_n_0m_pairs, sigf_means_n_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_n_2v1.png", width = 12, height = 10)



p <- plotNewPairs(sigf_means_tot_0m_pairs, sigf_means_tot_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_tot_2v1.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_r_0m_pairs, sigf_means_r_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_r_2v1.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_n_0m_pairs, sigf_means_n_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_n_2v1.png", width = 12, height = 10)

r_new_pairs_2v1 <- getNewPairs(df1 = sigf_means_r_0m_pairs, df2 = sigf_means_r_1m_pairs) %>% mutate(overall = "R")
n_new_pairs_2v1 <- getNewPairs(df1 = sigf_means_n_0m_pairs, df2 = sigf_means_n_1m_pairs) %>% mutate(overall = "N")


p <- rbind(r_new_pairs_2v1, n_new_pairs_2v1) %>% 
  ggplot(aes(overall,value, fill = overall)) + geom_boxplot(outlier.shape = 21) + 
  facet_wrap(~cluster2.x, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("N", "R"))) + scale_fill_manual(values = c("lightgrey", "salmon"))
ggsave(plot = p, "results/cellphonedb/analysis/boxplot_sigf_r_v_n_2v1.png", width = 15, height = 13)




p <- plotPairs(sigf_means_tot_1m_pairs, sigf_means_tot_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_tot_3v2.png", width = 12, height = 10)

p <- plotPairs(sigf_means_r_1m_pairs, sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_r_3v2.png", width = 12, height = 10)

p <- plotPairs(sigf_means_n_1m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_n_3v2.png", width = 12, height = 10)



p <- plotNewPairs(sigf_means_tot_1m_pairs, sigf_means_tot_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_tot_3v2.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_r_1m_pairs, sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_r_3v2.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_n_1m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_n_3v2.png", width = 12, height = 10)

r_new_pairs_3v2 <- getNewPairs(df1 = sigf_means_r_1m_pairs, df2 = sigf_means_r_3m_pairs) %>% mutate(overall = "R")
n_new_pairs_3v2 <- getNewPairs(df1 = sigf_means_n_1m_pairs, df2 = sigf_means_n_3m_pairs) %>% mutate(overall = "N")

p <- rbind(r_new_pairs_3v2, n_new_pairs_3v2) %>% 
  ggplot(aes(overall,value, fill = overall)) + geom_boxplot(outlier.shape = 21) + 
  facet_wrap(~cluster2.x, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("N", "R"))) + scale_fill_manual(values = c("lightgrey", "salmon"))
ggsave(plot = p, "results/cellphonedb/analysis/boxplot_sigf_r_v_n_3v2.png", width = 15, height = 13)








#### Heatmaps
p <- cast(sigf_means_tot_0m_pairs, cluster1 ~ cluster2, sum)
p <- p[,-1]
rownames(p) <- colnames(p)
p <- p %>% pheatmap::pheatmap(display_numbers = T, number_format ="%.0f", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0,40,1))), breaks = seq(0, 40, 1))
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_tot_0m.png", width = 8, height = 6)

p <- cast(sigf_means_tot_1m_pairs, cluster1 ~ cluster2, sum)
p <- p[,-1]
rownames(p) <- colnames(p)
p <- p %>% pheatmap::pheatmap(display_numbers = T, number_format ="%.0f", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0,40,1))), breaks = seq(0, 40, 1))
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_tot_1m.png", width = 8, height = 6)

p <- cast(sigf_means_tot_3m_pairs, cluster1 ~ cluster2, sum)
p <- p[,-1]
rownames(p) <- colnames(p)
p <- p %>% pheatmap::pheatmap(display_numbers = T, number_format ="%.0f", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0,40,1))), breaks = seq(0, 40, 1))
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_tot_3m.png", width = 8, height = 6)





p <- cast(sigf_means_r_0m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_0m.png", width = 8, height = 6)

p <- cast(sigf_means_r_1m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_1m.png", width = 8, height = 6)

p <- cast(sigf_means_r_3m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_3m.png", width = 8, height = 6)




p <- cast(sigf_means_n_0m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_0m.png", width = 8, height = 6)

p <- cast(sigf_means_n_1m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_1m.png", width = 8, height = 6)

p <- cast(sigf_means_n_3m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_3m.png", width = 8, height = 6)







p <- plotDiffHeatmap(df1 = sigf_means_tot_0m_pairs, df2 = sigf_means_tot_1m_pairs)
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_tot_2v1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_tot_1m_pairs, df2 = sigf_means_tot_3m_pairs)
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_tot_3v2.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_r_0m_pairs, df2 = sigf_means_r_1m_pairs)
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_r_2v1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_r_1m_pairs, df2 = sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cmanuscript/overall/heatmap_sigf_r_3v2.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_n_0m_pairs, sigf_means_n_1m_pairs)
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_n_2v1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_n_1m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/manuscript/overall/heatmap_sigf_n_3v2.png", width = 8, height = 6)








## What are the sigf interactions?
sigf_means_r_0m_sigf_pairs <- getSigfPairs(sigf_means_r_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_r_1m_sigf_pairs <- getSigfPairs(sigf_means_r_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_r_3m_sigf_pairs <- getSigfPairs(sigf_means_r_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))

sigf_means_n_0m_sigf_pairs <- getSigfPairs(sigf_means_n_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_n_1m_sigf_pairs <- getSigfPairs(sigf_means_n_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_n_3m_sigf_pairs <- getSigfPairs(sigf_means_n_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))


## What are the new interactions?
sigf_means_r_sigf_pairs_new <- getNewInteractions(sigf_means_r_0m_sigf_pairs, sigf_means_r_1m_sigf_pairs)
write.table(sigf_means_r_sigf_pairs_new, "results/cellphonedb/analysis/new_pairs_r.txt", sep = "\t", quote = F, row.names = F)

sigf_means_n_sigf_pairs_new <- getNewInteractions(sigf_means_n_0m_sigf_pairs, sigf_means_n_1m_sigf_pairs)
write.table(sigf_means_n_sigf_pairs_new, "results/cellphonedb/analysis/new_pairs_n.txt", sep = "\t", quote = F, row.names = F)

sigf_means_r_sigf_pairs_new %>% filter(cluster1 == "10 CD8 effector/exhausted") %>% write.table("results/cellphonedb/analysis/new_pairs_r_for_exhausted.txt", sep = "\t", quote = F, row.names = F)

