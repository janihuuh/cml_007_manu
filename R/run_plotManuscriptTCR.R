
## Follow antigen-specific clonotypes, ie clonotypes that are found more than once
clonotypes.to.keep <- as.data.frame(table(cml_seurat$new_clonotypes_id)) %>% filter(Freq>1) %>% pull(Var1)

clusters.to.rm <- c("14 CD8/B-cell doublet", "12 CD8 MAIT-like CXCR3+")
clusters.to.keep 

cells.to.keep    <- cml_seurat@meta.data %>% filter(new_clonotypes_id %in% clonotypes.to.keep) %>% 
  # filter(grepl("CD8", cluster) | grepl("CD4", cluster) | grepl("T-cell", cluster) & !grepl("B-cell", cluster)) %>% 
  filter(grepl("CD8", cluster)) %>% 
  filter(!cluster %in% clusters.to.rm) %>% pull(barcode)
clonotype_seurat <- subset(cml_seurat, cells = cells.to.keep)
clonotype_seurat <- clonotype_seurat %>% getLatentUMAP()



## 
cl_metadata <- clonotype_seurat@meta.data %>% group_by(orig.ident) %>% summarise(total_clonotype_reads = n())
meta <- clonotype_seurat@meta.data %>% left_join(cl_metadata)
rownames(meta) <- meta$barcode
clonotype_seurat@meta.data <- meta

cl_metadata <- clonotype_seurat@meta.data %>% group_by(orig.ident, timepoint, patient.x) %>% summarise(total_clonotype_reads = n())

clonotype_seurat@meta.data %>% 
  filter(temp == "y") %>% 
  group_by(orig.ident, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(cl_metadata) %>% mutate(prop=n/total_clonotype_reads) %>% 
  ggplot(aes(timepoint,prop,fill=cluster)) + geom_bar(stat="identity")



##
# top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% head(n=5) %>% pull(new_clonotypes_id)
top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 5, wt = n)  %>% pull(new_clonotypes_id)
# top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 30, wt = n)  %>% pull(new_clonotypes_id)
# top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 20, wt = n)  %>% pull(new_clonotypes_id)
# top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 10, wt = n)  %>% pull(new_clonotypes_id)

clonotype_seurat@meta.data %>% 
  filter(new_clonotypes_id %in% top_clonotypes) %>% 
  filter(grepl("CD8", cluster)) %>% 
  group_by(new_clonotypes_id, orig.ident, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(cl_metadata) %>% mutate(prop=n/total_clonotype_reads) %>% 
  mutate(timepoint = plyr::revalue(timepoint, replace = c("dg" = 0, "3mo" = 3, "12mo" = 12))) %>% mutate(timepoint = as.numeric(timepoint)) %>% 
  
  ggplot(aes(timepoint,prop,fill=cluster)) + 
  # geom_bar(stat="identity") + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  facet_wrap(~new_clonotypes_id,ncol=5, scales = "free_y") + facets_nice + theme_classic(base_size = 12)
ggsave("results/manuscript/tcr/area_tcr.pdf", width = 11, height = 8)


df <- clonotype_seurat@meta.data %>% 
  filter(new_clonotypes_id %in% top_clonotypes) %>% 
  filter(grepl("CD8", cluster)) %>% 
  group_by(new_clonotypes_id, orig.ident, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(cl_metadata) %>% mutate(prop=n/total_clonotype_reads) %>% 
  mutate(timepoint = plyr::revalue(timepoint, replace = c("dg" = 0, "3mo" = 3, "12mo" = 12))) %>% mutate(timepoint = as.numeric(timepoint)) %>% ungroup() %>% as.data.frame()

df %>% 
  ggplot(aes(as.factor(timepoint),prop)) + geom_boxplot() + facet_wrap(~cluster, scales = "free") + ggpubr::stat_compare_means()


dcast(df, new_clonotypes_id~cluster+timepoint, value.var = "prop")

clone_df <- dcast(df, paste(new_clonotypes_id,timepoint)~cluster, value.var = "prop")
clone_df[is.na(clone_df)] <- 0
clone_df_mtx <- clone_df[,-1]

clone_df_pca <- prcomp((clone_df_mtx))
clone_df_pca_df <- clone_df_pca$x %>% as.data.frame() %>% mutate(clone = clone_df$`paste(new_clonotypes_id, timepoint)`) #%>% mutate(clone = strsplit(clone, split = " ")[1])

ggplot(clone_df_pca_df, aes(PC1,PC2)) + geom_point()
 


## Pseudobulk
cells.to.keep <- clonotype_seurat@meta.data %>% filter(new_clonotypes_id %in% top_clonotypes) %>% filter(grepl("CD8", cluster)) %>% pull(barcode)
top20_clonotypes_seurat <- subset(clonotype_seurat, cells = cells.to.keep) %>% preprocessSeurat(cells.to.use = cells.to.keep)

## Only HVGs
hvg     <- VariableFeatures(top20_clonotypes_seurat)
counts  <- as.data.frame(top20_clonotypes_seurat@assays$RNA@counts[rownames(top20_clonotypes_seurat) %in% hvg, ])
counts2 <- t(counts) %>% as.data.frame()
counts2$key <- paste(top20_clonotypes_seurat$new_clonotypes_id, top20_clonotypes_seurat$timepoint)
# counts2$patient   <- top20_clonotypes_seurat$patient.x
# counts2$timepoint <- top20_clonotypes_seurat$timepoint
counts_sum <- lapply(unique(counts2$key), FUN = function(x) counts2 %>% filter(key == x) %>% dplyr::select(-key) %>% colSums() %>% as.data.frame) # %>% rbindlist()
counts_sum <- counts_sum %>% Laurae::cbindlist() 
counts_sum <- counts_sum %>% as.data.frame()

rownames(counts_sum) <- colnames(counts2)[-1]
colnames(counts_sum) <- unique(counts2$key)
counts_sum_norm <- apply(counts_sum, 2, function(x) as.data.frame(x/colSums(counts_sum))) %>% Laurae::cbindlist()
rownames(counts_sum_norm) <- colnames(counts2)[-1]

counts_sum_norm_pca <- prcomp(t(counts_sum_norm))
counts_sum_norm_pca_df <- data.frame(counts_sum_norm_pca$x)

ggplot(counts_sum_norm_pca_df, aes(PC1,PC2)) + geom_point()




#### Look at which clonotypes expand 
top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 50, wt = n)  %>% pull(new_clonotypes_id)
df <- clonotype_seurat@meta.data %>% 
  filter(new_clonotypes_id %in% top_clonotypes) %>% 
  filter(grepl("CD8", cluster)) %>% 
  group_by(new_clonotypes_id, orig.ident, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(cl_metadata) %>% mutate(prop=n/total_clonotype_reads) %>% 
  mutate(timepoint = plyr::revalue(timepoint, replace = c("dg" = 0, "3mo" = 3, "12mo" = 12))) %>% mutate(timepoint = as.numeric(timepoint)) %>% ungroup() %>% as.data.frame()

df_0   <- df %>% filter(timepoint == 1) %>% group_by(new_clonotypes_id) %>% summarise(prop = sum(prop))
df_3   <- df %>% filter(timepoint == 2) %>% group_by(new_clonotypes_id) %>% summarise(prop = sum(prop))
df_12  <- df %>% filter(timepoint == 3) %>% group_by(new_clonotypes_id) %>% summarise(prop = sum(prop))
df_exp <- df_0 %>% left_join(df_3, by = "new_clonotypes_id") %>% left_join(df_12, by = "new_clonotypes_id") %>% mutate(log2fc_3v0 = log2(prop.y/prop.x), log2fc_12v3 = log2(prop/prop.y))

## 3v0
df_expanded <- df_exp %>% filter(log2fc_3v0 > 1) %>% pull(new_clonotypes_id)
df_shrunk   <- df_exp %>% filter(log2fc_3v0 < -1) %>% pull(new_clonotypes_id)
df_norm     <- df_exp %>% filter(log2fc_3v0 > -1 & log2fc_3v0 < 1) %>% pull(new_clonotypes_id)

df_0_exp    <- df %>% filter(new_clonotypes_id %in% df_expanded) %>% filter(timepoint == 1) %>% mutate(type = "expanded") %>% group_by(new_clonotypes_id) %>% mutate(prop = prop/sum(prop))
df_0_shrunk <- df %>% filter(new_clonotypes_id %in% df_shrunk) %>% filter(timepoint == 1) %>% mutate(type = "shrunk") %>% group_by(new_clonotypes_id) %>% mutate(prop = prop/sum(prop))
df_0_norm   <- df %>% filter(new_clonotypes_id %in% df_norm) %>% filter(timepoint == 1) %>% mutate(type = "normal") %>% group_by(new_clonotypes_id) %>% mutate(prop = prop/sum(prop))

rbind(df_0_exp, df_0_norm, df_0_shrunk) %>% 
  ggplot(aes(type,prop,fill=type)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, ncol = 5) + ggpubr::stat_compare_means(label = "p.format") + 
  theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) + theme(legend.position = "none") + labs(x = "", y = "prop of clonotype at baseline") + facets_nice + geom_jitter(size = 0.5)
ggsave("results/manuscript/tcr/box_3v0_expansion.pdf", width = 10, height = 4)

rbind(df_0_exp, df_0_norm, df_0_shrunk) %>% 
  filter(cluster %in% c("1 CD8 EM GZMK+ CXCR3+", "2 CD8 cytotoxic GZMH+", "5 CD8 EMRA ZNF683+")) %>% 
  ggplot(aes(type,prop,fill=type)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, ncol = 5) + ggpubr::stat_compare_means(label = "p.format") + 
  theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) + theme(legend.position = "none") + labs(x = "", y = "prop of clonotype at baseline") + facets_nice + geom_jitter(size = 0.5) +
  ylim(values = c(0,1.19))
ggsave("results/manuscript/tcr/box_3v0_expansion_ms.pdf", width = 6.5, height = 3.5)
  

## 12v3
df_expanded <- df_exp %>% filter(log2fc_12v3 > 1) %>% pull(new_clonotypes_id)
df_shrunk   <- df_exp %>% filter(log2fc_12v3 < -1) %>% pull(new_clonotypes_id)
df_norm     <- df_exp %>% filter(log2fc_12v3 > -1 & log2fc_12v3 < 1) %>% pull(new_clonotypes_id)

df_0_exp    <- df %>% filter(new_clonotypes_id %in% df_expanded) %>% filter(timepoint == 1) %>% mutate(type = "expanded") %>% group_by(new_clonotypes_id) %>% mutate(prop = prop/sum(prop))
df_0_shrunk <- df %>% filter(new_clonotypes_id %in% df_shrunk) %>% filter(timepoint == 1) %>% mutate(type = "shrunk") %>% group_by(new_clonotypes_id) %>% mutate(prop = prop/sum(prop))
df_0_norm   <- df %>% filter(new_clonotypes_id %in% df_norm) %>% filter(timepoint == 1) %>% mutate(type = "normal") %>% group_by(new_clonotypes_id) %>% mutate(prop = prop/sum(prop))

rbind(df_0_exp, df_0_norm, df_0_shrunk) %>% 
  ggplot(aes(type,prop,fill=type)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, ncol = 5) + ggpubr::stat_compare_means(label = "p.format") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) + theme(legend.position = "none") + labs(x = "", y = "prop of clonotype at 3 months") + facets_nice + geom_jitter(size = 0.5)
ggsave("results/manuscript/tcr/box_12v3_expansion.pdf", width = 10, height = 4)

rbind(df_0_exp, df_0_norm, df_0_shrunk) %>% 
  filter(cluster %in% c("5 CD8 EMRA ZNF683+", "9 CD8 EMRA IFNG+")) %>% 
  ggplot(aes(type,prop,fill=type)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, ncol = 5) + ggpubr::stat_compare_means(label = "p.format") + 
  theme_classic(base_size = 15) +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) + theme(legend.position = "none") + labs(x = "", y = "prop of clonotype at 3mo") + facets_nice + geom_jitter(size = 0.5) +
  ylim(values = c(0,1.19))
ggsave("results/manuscript/tcr/box_12v3_expansion_ms.pdf", width = 5, height = 3.75)




### TCRGP
tcrgp_df <- clonotype_seurat@meta.data %>% dplyr::select(new_clonotypes_id, trb_cdr3s_aa, trb_cdr3s_nt, v_trb, d_trb, j_trb)
fwrite(tcrgp_df, "results/manuscript/tcr/tcrgp_df.txt", sep = "\t", quote = F, row.names = F)
# tcrgp_df <- fread("results/manuscript/tcr/tcrgp_df.txt")
metadata <- clonotype_seurat@meta.data

## Find antigen-specific TCRs, select 0.95 as a threshold
ag_specific           <- lapply(list.files("results/manuscript/tcr/tcrgp/", full.names = T), function(x) fread(x)) %>% bind_cols()
colnames(ag_specific) <- do.call(lapply(list.files("results/manuscript/tcr/tcrgp/", full.names = T), extractFileName) , what = "c") %>% gsub(pattern = ".csv", replacement =  "") %>% gsub(pattern = "007_clonotype_", replacement =  "")
ag_specific[ag_specific < 0.9] <- 0
ag_specific_df   <- metadata %>% bind_cols(ag_specific)
ag_specific_filt <- ag_specific_df # %>% filter(count > 1)

preds_df <- ag_specific_filt %>% dplyr::select(ELAGIGILTV_cdr3b_comb:YVLDHLIVV_cdr3b)
target <- apply(preds_df, 1, function(x) ifelse(max(x) == 0, 16, which.max(x)))
ag_specific_filt$target <- target
ag_specific_filt$target <- c(colnames(preds_df), "none")[ag_specific_filt$target]
ag_specific_filt$target[rowSums(preds_df) > 1] <- "multi"
ag_specific_filt$target[is.na(ag_specific_filt$target)] <- "none"

data.frame(ag_specific_filt$pr1_cdr3b, ag_specific_filt$trb_cdr3s_aa) %>% View
clonotype_seurat$target <- ag_specific_filt$target

tcrgp_res_df <- ag_specific_filt %>% dplyr::select(v_trb, trb_cdr3s_aa, GILGFVFTL_cdr3b:IPSINVHHY_cdr3b, NLVPMVATV_cdr3b, RAKFKQLL_cdr3b:target)
fwrite(tcrgp_res_df, "results/manuscript/tcr/tcrgp_res_df.txt", sep = "\t", quote = F, row.names = F)

table(clonotype_seurat$target)

target_df <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id, target) %>% summarise(n=n()) %>% dplyr::select(-n)

target_df %>% group_by(target) %>% summarise(n=n())

clonotype_seurat@meta.data %>% 
  filter(grepl("CD8", cluster)) %>% 
  group_by(new_clonotypes_id, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(target_df) %>% 

  ggplot(aes(target,prop)) + geom_boxplot() + facet_wrap(~cluster) + ggpubr::rotate_x_text(angle = 45)


clonotype_seurat@meta.data %>% 
  filter(grepl("CD8", cluster)) %>% 
  group_by(new_clonotypes_id, timepoint, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(target_df) %>% 
  
  ggplot(aes(cluster,prop)) + geom_boxplot() + facet_wrap(~target) + ggpubr::rotate_x_text(angle = 45)


df <- clonotype_seurat@meta.data %>% group_by(orig.ident, target, .drop = FALSE) %>% summarise(n=n()) %>% 
  filter(target == "NLVPMVATV_cdr3b") 

cl_metadata %>% left_join(df) %>% mutate(n = ifelse(is.na(n),0,n)) %>%
  mutate(prop=n/total_clonotype_reads) %>% 
  ggplot(aes(timepoint,prop,fill=timepoint)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.7) + scale_fill_manual(values=getPalette3(4)) + theme_classic(base_size = 12) + theme(legend.position = "none") +
  labs(x = "", y = "prop of \nanti-CMV clonotypes in repertoire")
ggsave("results/manuscript/tcr/box_nlv.pdf", width = 3, height = 4)


clonotype_seurat@meta.data %>% 
  filter(target == "NLVPMVATV_cdr3b") %>% 
  ggplot(aes(timepoint,cyto1)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label="p") + facet_wrap(~patient.x)



clonotype_seurat$target_nlv <- ifelse(clonotype_seurat$target == "NLVPMVATV_cdr3b", "yes", "no")
cells.to.highlight <- clonotype_seurat@meta.data %>% filter(target == "NLVPMVATV_cdr3b") %>% pull(barcode)
  DimPlot(clonotype_seurat, cells.highlight = cells.to.highlight, split.by = "timepoint")

clonotype_seurat@meta.data %>% 
  filter(target == "NLVPMVATV_cdr3b") %>% 
  group_by(new_clonotypes_id) %>% summarise(n=n()) 


clonotype_seurat@meta.data %>% 
  group_by(new_clonotypes_id, target) %>% summarise(n=n()) 



patient_df <- clonotype_seurat@meta.data %>% group_by(orig.ident,timepoint,patient.x) %>% summarise(n=n()) %>% dplyr::select(-n)

df_temp <- clonotype_seurat@meta.data %>% 
  filter(target == "NLVPMVATV_cdr3b") %>% 
  group_by(orig.ident,cluster,.drop = FALSE) %>% summarise(n=n(),.drop = FALSE) %>% mutate(prop=n/sum(n)) %>% 
  left_join(cl_metadata %>% group_by(orig.ident,.drop = FALSE) %>% summarise(tot=sum(total_clonotype_reads),.drop = FALSE), by = "orig.ident") %>% mutate(prop_reads = n/tot) %>% left_join(patient_df)

df_temp %>% 
  filter(grepl("CD8", cluster)) %>% 
  filter(!cluster %in% clusters.to.rm) %>% 
  # mutate(pattern = factor(as.character(pattern), levels = c("FGNT", "QTGT", "GSQP"))) %>% 
  # mutate(timepoint = plyr::revalue(timepoint, replace = c("dg" = 1, "3mo" = 2, "12mo" = 3))) %>% 
  mutate(patient.x = factor(as.character(patient.x), levels = c("720", "730", "706"))) %>% 
  ggplot(aes(timepoint,prop,fill=cluster)) + geom_area(aes(group=cluster, fill=cluster),colour="gray30") + facet_wrap(~patient.x, ncol=3) + scale_fill_manual(values = getPalette(7)) +
  theme_classic(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "top") + facets_nice + guides(fill=guide_legend(nrow=2,byrow=TRUE,label.hjust = 0)) + labs(fill = "") + labs(x = "", y = "prop of anti-CMV repertoire")
ggsave("results/manuscript/tcr/area_nlv_tcr.pdf", width = 8, height = 4)



clonotype_seurat@meta.data %>% group_by(orig.ident, target, .drop = FALSE) %>% summarise(n=n()) %>% 
  filter(target == "NLVPMVATV_cdr3b") 

clonotype_seurat$cmv_pos <- ifelse(clonotype_seurat$patient.x == "730", "CMVunknown", "CMVpos")














## GLIPH2
gliph_df <- clonotype_seurat@meta.data %>% dplyr::rename(v=v_trb, d=d_trb, j=j_trb,cdr3nt=trb_cdr3s_nt,cdr3aa=trb_cdr3s_aa,name=patient.x) %>% mutate(count=10,freq=NA,VEnd=NA, DStart=NA,DEnd=NA,JStart=NA) %>% dplyr::select(count,freq,cdr3nt,cdr3aa,v,d,j,VEnd,DStart,DEnd,JStart,name) %>% vdjToGliph()
gliph_df <- gliph_df[!duplicated(gliph_df$CDR3b), ]
gliph_df <- gliph_df[gliph_df$CDR3b != "", ]
fwrite(gliph_df, "results/manuscript/tcr/gliph_df.txt", sep = "\t", quote = F, row.names = F)

gliph_df <- readGliphFile("results/manuscript/tcr/P1032_X5AWQ50KS8.csv") %>% filter(number_unique_cdr3 > 2 & vb_score < 0.05)
gliph_df <- readGliphFile("results/manuscript/tcr/P1032_X5AWQ50KS8.csv") %>% filter(vb_score < 0.05)
gliph_df <- readGliphFile("results/manuscript/tcr/P1032_X5AWQ50KS8.csv") %>% filter(number_unique_cdr3 > 2)

gliph_df <- gliph_df[!is.na(gliph_df$index), ] 
gliph_to_use <- gliph_df %>% dplyr::select(pattern,TcRb) %>% dplyr::rename(trb_cdr3s_aa = TcRb)




## Logoplot
require(ggseqlogo)

plotLogo <- function(seq){
  
  ## Determine the commonest sequence
  max_n <- table(nchar(seq)) %>% which.max %>% names %>% as.numeric()
  seq_cut <- seq[nchar(seq) == max_n]
  seq_cut <- substr(seq_cut, 4, nchar(seq_cut) - 3)
  
  # ggplot() + geom_logo(seq_cut,  method = 'prob') + theme_logo()
  ggplot() + geom_logo(seq_cut,  method = 'bits') + theme_logo()
  
}

plotLogoList <- function(seq_list){
  
  ## Determine the commonest sequence
  getSeq <- function(seq){
    
    max_n <- table(nchar(seq)) %>% which.max %>% names %>% as.numeric()
    seq_cut <- seq[nchar(seq) == max_n]
    seq_cut <- substr(seq_cut, 4, nchar(seq_cut) - 3)
    
  }
  
  seq_list_cut <- lapply(seq_list, getSeq)
  
  ggseqlogo(seq_list_cut,  method = 'bits', ncol=5) + theme_logo()
  
}


logolist <- split(gliph_df, f = gliph_df$index)[1:4]
names(logolist) <- lapply(logolist, FUN = function(x) unique(x$pattern))
seq_list <- lapply(logolist, FUN = function(x) x %>% pull(TcRb))
plotLogoList(seq_list)
ggsave("results/manuscript/tcr/logoplots_top3.pdf", width = 12, height = 8)

logolist <- split(pr1_gliph, f = pr1_gliph$index)[1:10]
names(logolist) <- lapply(logolist, FUN = function(x) unique(x$pattern))
seq_list <- lapply(logolist, FUN = function(x) x %>% pull(TcRb))
plotLogoList(seq_list)
ggsave("results/tcr/pr1/logoplots_top10.pdf", width = 12, height = 4)




meta <- clonotype_seurat@meta.data %>% left_join(gliph_to_use)
rownames(meta) <- meta$barcode
clonotype_seurat@meta.data <- meta

df_temp <- clonotype_seurat@meta.data %>% 
  filter(grepl("CD8", cluster)) %>% 
  filter(!is.na(pattern)) %>% 
  group_by(new_clonotypes_id,timepoint,cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% 
  left_join(meta, by = "new_clonotypes_id") %>% 
  dplyr::select(new_clonotypes_id,timepoint.x,pattern,prop,cluster.x)
  
df_temp <- df_temp[!duplicated(df_temp), ]

df_temp %>% ggplot(aes(pattern,prop)) + geom_boxplot(outlier.shape = NA) + facet_wrap(timepoint.x~cluster.x, ncol=5) + ggpubr::rotate_x_text(angle = 45) + ggpubr::stat_compare_means(label="p.format") + geom_jitter(size = 0.5)

df_temp %>% ggplot(aes(cluster.x,prop)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~pattern, ncol=5) + ggpubr::rotate_x_text(angle = 45) + ggpubr::stat_compare_means(label="p.format") + geom_jitter(size = 0.5)

df_temp %>% ggplot(aes(cluster.x,prop)) + geom_boxplot(outlier.shape = NA) + facet_wrap(timepoint.x~pattern, ncol=3) + ggpubr::rotate_x_text(angle = 45) + ggpubr::stat_compare_means(label="p.format") + geom_jitter(size = 0.5)


df_temp <- clonotype_seurat@meta.data %>% 
  filter(grepl("CD8", cluster)) %>% 
  filter(!is.na(pattern)) %>% 
  group_by(pattern, timepoint,cluster,.drop = FALSE) %>% summarise(n=n(),.drop = FALSE) %>% mutate(prop=n/sum(n)) %>% 
  left_join(cl_metadata %>% group_by(timepoint,.drop = FALSE) %>% summarise(tot=sum(total_clonotype_reads),.drop = FALSE), by = "timepoint") %>% mutate(prop_reads = n/tot)
# cl_df <- cl_metadata %>% group_by(timepoint) %>% summarise(tot=sum(total_clonotype_reads))

df_temp %>% ggplot(aes(timepoint,prop,fill=cluster)) + geom_bar(stat="identity") + facet_wrap(~pattern, ncol=3) + ggpubr::rotate_x_text(angle = 45) 


df_temp %>% 
  filter(grepl("CD8", cluster)) %>% 
  filter(!cluster %in% clusters.to.rm) %>% 
  mutate(pattern = factor(as.character(pattern), levels = c("FGNT", "QTGT", "GSQP"))) %>% 
  # mutate(timepoint = plyr::revalue(timepoint, replace = c("dg" = 1, "3mo" = 2, "12mo" = 3))) %>% 
  ggplot(aes(timepoint,prop,fill=cluster)) + geom_area(aes(group=cluster, fill=cluster),colour="gray30") + facet_wrap(~pattern, ncol=3) + scale_fill_manual(values = getPalette(7)) +
  theme_classic(base_size = 12) + ggpubr::rotate_x_text(angle = 45) + theme(legend.position = "top") + facets_nice + guides(fill=guide_legend(nrow=2,byrow=TRUE,label.hjust = 0)) + labs(fill = "")
ggsave("results/manuscript/tcr/area_gliph.pdf", width = 8, height = 4)


df_temp %>% ggplot(aes(timepoint,prop_reads,fill=cluster)) + geom_bar(stat="identity") + facet_wrap(~pattern, ncol=3) + ggpubr::rotate_x_text(angle = 45) 

df_temp %>% ggplot(aes(timepoint,n,fill=cluster)) + geom_bar(stat="identity") + facet_wrap(~pattern, ncol=3) + ggpubr::rotate_x_text(angle = 45) 


## Do the patterns expand or suppress
clonotype_seurat@meta.data %>% 
  filter(grepl("CD8", cluster)) %>% 
  filter(!is.na(pattern)) %>% 
  group_by(pattern,orig.ident) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% left_join(cl_metadata) %>% mutate(prop=n/total_clonotype_reads) %>% 
  # group_by(pattern,timepoint) %>% summarise(prop=sum(prop)) %>% 
  
  ggplot(aes(timepoint,prop)) + geom_point() + facet_wrap(~pattern) + geom_path(aes(group=patient.x))






####### pseudotime

DimPlot(clonotype_seurat, label = T, repel = T, cols = getPalette5(8)) + theme(legend.position = "none")
ggsave("results/manuscript/tcr/umap_clonotype.png", width = 5, height = 5)



## Calculate pseudotimes
runSlingshot <- function(seurat_object, cells.to.keep, reducedDim = "LATENT_UMAP"){
  
  require(slingshot); require(SummarizedExperiment)
  seurat_temp <- subset(seurat_object, cells = cells.to.keep) 
  timepoint   <- seurat_temp$timepoint %>% as.character() %>% unique()
  timepoint   <- timepoint[1]
  sce_temp    <- as.SingleCellExperiment(seurat_temp)
  sce_temp    <- slingshot(data = sce_temp, clusterLabels = 'cluster', reducedDim = reducedDim, start.clus = "8 CD8 naive TCF7+")
  sling_curve <- SlingshotDataSet(sce_temp)
  return(sling_curve)
  
}

### Get pseudotimes for each clonotype
getClonotypeCurves <- function(seurat_object, clonotype){
  
  message(clonotype)
  cells.to.keep <- seurat_object@meta.data %>% filter(new_clonotypes_id == clonotype) %>% pull(barcode)
  curve         <- runSlingshot(seurat_object, cells.to.keep = cells.to.keep)
  
}



cololors                <- clonotype_seurat$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 5, wt = n)  %>% pull(new_clonotypes_id)
curves     <- lapply(top_clonotypes, getClonotypeCurves, seurat_object = clonotype_seurat)
curve_cols <- getPalette(25)

png("results/manuscript/tcr/sling_clonotype.png", width = 1024, height = 1024)

par(mfrow=c(4,5))
for(i in 1:length(curves)){
  
  message(i)
  plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
  # text(0,3,"8 CD8 naive TCF7+")
  # text(-3,0,"5 CD8 EMRA ZNF683+")
  # text(3,0,"2 CD8 cytotoxic GZMH+")
  # text(0,-3,"1 CD8 EM GZMK+ CXCR3+")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()





top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(patient = substr(new_clonotypes_id, 1,3)) %>% group_by(patient) %>% top_n(n = 20, wt = n)  %>% pull(new_clonotypes_id)
top_clonotypes[58] <- top_clonotypes[59]
top_clonotypes <- top_clonotypes[order(top_clonotypes)]

curves     <- lapply(top_clonotypes, getClonotypeCurves, seurat_object = clonotype_seurat)
curve_cols <- getPalette(length(top_clonotypes))

png("results/manuscript/tcr/sling_clonotype_full.png", width = 1024*5, height = 1024)

par(mfrow=c(4,20))
for(i in 1:length(curves)){
  
  message(i)
  plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
  # text(0,3,"8 CD8 naive TCF7+")
  # text(-3,0,"5 CD8 EMRA ZNF683+")
  # text(3,0,"2 CD8 cytotoxic GZMH+")
  # text(0,-3,"1 CD8 EM GZMK+ CXCR3+")
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()

## 8-5-2 was the most frequenct pattern, focus on that
clones_852 <- c("706_clonotype11", "706_clonotype15", "706_clonotype9",
                "716_clonotype10", "716_clonotype13", "716_clonotype5", 
                "720_clonotype14", "720_clonotype15", "720_clonotype19", 
                "730_clonotype5", "730_clonotype6")

curves     <- lapply(clones_852, getClonotypeCurves, seurat_object = clonotype_seurat)
curve_cols <- getPalette(length(clones_852))

curve_cols = c("black", "blue", "red")
png("results/manuscript/tcr/sling_clonotype_852_helper.png", width = 1024, height = 1024)

par(mfrow=c(4,4))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = clones_852[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[j])
  }
  
}
dev.off()


## Remove "wrong" curves
curves[[1]]@curves[[1]] <- NULL
curves[[2]]@curves[[2]] <- NULL
curves[[3]]@curves[[2]] <- NULL
curves[[4]]@curves[[3]] <- NULL
curves[[4]]@curves[[2]] <- NULL
curves[[5]]@curves[[2]] <- NULL
curves[[6]]@curves[[2]] <- NULL
curves[[10]]@curves[[1]] <- NULL
curves[[11]]@curves[[2]] <- NULL

curve_cols <- getPalette(length(clones_852))
png("results/manuscript/tcr/sling_clonotype_852.png", width = 1024*1/2, height = 1024*1/2)

plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors) #, main = top_clonotypes[i])
text(0,3,"8")
text(-3,0,"5")
text(3,0,"2")
text(0,-3,"1")
for(i in 1:length(curves)){
  message(i)
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
}

dev.off()









## 8-5-1
clones_851 <- c("706_clonotype1", "706_clonotype11", "706_clonotype19", "706_clonotype4",
                "716_clonotype15", "716_clonotype7", 
                "720_clonotype6",
                "730_clonotype8")

curves     <- lapply(clones_851, getClonotypeCurves, seurat_object = clonotype_seurat)

curve_cols <- c("black", "blue", "red")
png("results/manuscript/tcr/sling_clonotype_851_helper.png", width = 1024, height = 1024)

par(mfrow=c(4,4))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = clones_851[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[j])
  }
  
}
dev.off()


## Remove "wrong" curves
curves[[1]]@curves[[1]] <- NULL
curves[[2]]@curves[[2]] <- NULL
curves[[3]]@curves[[1]] <- NULL
curves[[4]]@curves[[3]] <- NULL
curves[[4]]@curves[[2]] <- NULL
curves[[5]]@curves[[2]] <- NULL
curves[[6]]@curves[[1]] <- NULL
curves[[7]]@curves[[3]] <- NULL
curves[[7]]@curves[[2]] <- NULL
curves[[8]]@curves[[2]] <- NULL

curve_cols <- getPalette(length(clones_851))
png("results/manuscript/tcr/sling_clonotype_851.png", width = 1024*1/2, height = 1024*1/2)

plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors) #, main = top_clonotypes[i])
text(0,3,"8")
text(-3,0,"5")
text(3,0,"2")
text(0,-3,"1")
for(i in 1:length(curves)){
  message(i)
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
}

dev.off()














####### at 0months
clonotype_seurat_0m  <- subset(clonotype_seurat, timepoint == "dg")
clonotype_seurat_3m  <- subset(clonotype_seurat, timepoint == "3mo")
clonotype_seurat_12m <- subset(clonotype_seurat, timepoint == "12mo")

top_clonotypes_0m <- clonotype_seurat_dg@meta.data %>% filter(timepoint == "dg") %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% filter(n>30) %>% pull(new_clonotypes_id)
top_clonotypes_0m <- top_clonotypes_0m[order(top_clonotypes_0m)]

curves     <- lapply(top_clonotypes_0m, getClonotypeCurves, seurat_object = clonotype_seurat_0m)
curve_cols <- getPalette(length(top_clonotypes_0m))

cololors                <- clonotype_seurat_0m$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

png("results/manuscript/tcr/sling_clonotype_0m.png", width = 1024*5, height = 1024*1/2)

length(curves)
par(mfrow=c(2,21))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat_0m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()









####### at 3months
clonotype_seurat_3m <- subset(clonotype_seurat, timepoint == "3mo")
top_clonotypes_3m   <- clonotype_seurat@meta.data %>% filter(timepoint == "3mo") %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% filter(n>30) %>% pull(new_clonotypes_id)
top_clonotypes_3m   <- top_clonotypes_3m[order(top_clonotypes_3m)]

curves     <- lapply(top_clonotypes_3m, getClonotypeCurves, seurat_object = clonotype_seurat_3m)
curve_cols <- getPalette(length(top_clonotypes_3m))

cololors                <- clonotype_seurat_3m$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

png("results/manuscript/tcr/sling_clonotype_3m.png", width = 1024*5, height = 1024*1/2)

length(curves)
par(mfrow=c(2,21))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat_3m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()





####### at 12months
clonotype_seurat_12m <- subset(clonotype_seurat, timepoint == "12mo")
top_clonotypes_12m   <- clonotype_seurat_dg@meta.data %>% filter(timepoint == "12mo") %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% filter(n>30) %>% pull(new_clonotypes_id)
top_clonotypes_12m   <- top_clonotypes_12m[order(top_clonotypes_12m)]

curves     <- lapply(top_clonotypes_12m, getClonotypeCurves, seurat_object = clonotype_seurat_12m)
curve_cols <- getPalette(length(top_clonotypes_12m))

cololors                <- clonotype_seurat_12m$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

png("results/manuscript/tcr/sling_clonotype_12m.png", width = 1024*5, height = 1024*1/2)

length(curves)
par(mfrow=c(2,21))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat_12m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()


## baseline common trajes
clonotypes_852 <- c("720_clonotype1",
                    "720_clonotype3",
                    "720_clonotype12",
                    
                    "730_clonotype10",
                    "730_clonotype10",
                    
                    "716_clonotype1",
                    "716_clonotype7",
                    "716_clonotype2",
                    "716_clonotype23")

curves     <- lapply(clonotypes_852, getClonotypeCurves, seurat_object = clonotype_seurat_0m)
curve_cols <- getPalette(length(clonotypes_852))

cololors                <- clonotype_seurat_0m$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

png("results/manuscript/tcr/sling_clonotype_0m_common_detailed.png", width = 1024*1/2, height = 1024*1/4)

length(curves)
par(mfrow=c(2,5))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat_0m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = clonotypes_852[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()


png("results/manuscript/tcr/sling_clonotype_0m_common.png", width = 1024*1/2, height = 1024*1/2)

par(mfrow=c(1,1))

plot(clonotype_seurat_0m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
text(0,3,"8")
text(-3,0,"5")
text(3,0,"2")
text(0,-3,"1")

for(i in 1:length(curves)){
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
}
dev.off()






## 3mo common trajes
clonotypes_82 <- c("730_clonotype2",
                   "730_clonotype3",
                   
                   "730_clonotype7",
                   "720_clonotype9",
                   "716_clonotype1",
                   "716_clonotype4",
                   
                   "716_clonotype1",
                   "730_clonotype10",
                   "720_clonotype12",
                   "706_clonotype3",
                   
                   "730_clonotype12",
                   "716_clonotype16",
                   "716_clonotype3",
                   
                   "716_clonotype9",
                   "716_clonotype12",
                   "716_clonotype8",
                   "716_clonotype23")

curves     <- lapply(clonotypes_82, getClonotypeCurves, seurat_object = clonotype_seurat_3m)
curve_cols <- getPalette(length(clonotypes_82))

cololors                <- clonotype_seurat_3m$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

png("results/manuscript/tcr/sling_clonotype_3m_common_detailed.png", width = 1024*1/2, height = 1024*1/2)

length(curves)
par(mfrow=c(4,5))
for(i in 1:length(curves)){
  
  plot(clonotype_seurat_3m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = clonotypes_82[i])
  text(0,3,"8")
  text(-3,0,"5")
  text(3,0,"2")
  text(0,-3,"1")
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
  
}
dev.off()


png("results/manuscript/tcr/sling_clonotype_3m_common.png", width = 1024*1/2, height = 1024*1/2)

par(mfrow=c(1,1))

plot(clonotype_seurat_3m@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
text(0,3,"8")
text(-3,0,"5")
text(3,0,"2")
text(0,-3,"1")

for(i in 1:length(curves)){
  for(j in 1:length(curves[[i]]@curves)){
    lines(curves[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
  }
}
dev.off()


# clonotypes_85 <- c("730_clonotype1",
#                    "720_clonotype1",
#                    "730_clonotype3",
#                    "720_clonotype4",
#                    "730_clonotype4")



clonotypes_821<- c("720_clonotype2",
                   "720_clonotype6",
                   "706_clonotype1",
                   "720_clonotype12",
                   "716_clonotype6",
                   "730_clonotype6",
                   "720_clonotype12",
                   "716_clonotype6",
                   "730_clonotype12",
                   "716_clonotype16",
                   "716_clonotype3",
                   "716_clonotype11",
                   "716_clonotype8",
                   "716_clonotype23")











#### Follow clonotypes in different time points; how do they evolve?
top_clonotypes <- clonotype_seurat@meta.data  %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% head(10) %>% pull(new_clonotypes_id)
# top_clonotypes <- top_clonotypes[order(top_clonotypes)]

curves_0m  <- lapply(top_clonotypes, getClonotypeCurves, seurat_object = clonotype_seurat_0m)
curves_3m  <- lapply(top_clonotypes, getClonotypeCurves, seurat_object = clonotype_seurat_3m)
curves_12m <- lapply(top_clonotypes, getClonotypeCurves, seurat_object = clonotype_seurat_12m)
curve_cols <- getPalette(length(top_clonotypes))

cololors                <- clonotype_seurat$cluster %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 1] <- getPalette5(8)[1]
cololors[cololors == 2] <- getPalette5(8)[2]
cololors[cololors == 5] <- getPalette5(8)[3]
cololors[cololors == 8] <- getPalette5(8)[4]
cololors[cololors == 9] <- getPalette5(8)[5]

n_clones=10
png("results/manuscript/tcr/sling_clonotype_tot.png", width = 1024, height = 1024*3)

par(mfrow=c(n_clones,3))

for(i in 1:n_clones){
  
  for(k in 1:3){
    
    if(k==1){
      plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
      for(j in 1:length(curves_0m[[i]]@curves)){
        lines(curves_0m[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
      }
    }
    
    if(k==2){
      plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
      for(j in 1:length(curves_3m[[i]]@curves)){
        lines(curves_3m[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
      }
    }
    
    if(k==3){
      plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
      for(j in 1:length(curves_12m[[i]]@curves)){
        lines(curves_12m[[i]]@curves[[j]], lwd = 4, col = curve_cols[i])
      }
    }

    
  }
  
}
dev.off()




curve_cols = c("black", "darkred", "navyblue")
png("results/manuscript/tcr/sling_clonotype_tot.png", width = 1024, height = 1024)

par(mfrow=c(4,3))

for(i in 1:n_clones){
  
  message(i)
  plot(clonotype_seurat@reductions$latent_umap@cell.embeddings, pch = 16, asp = 1, col = cololors, main = top_clonotypes[i])
  for(j in 1:length(curves_0m[[i]]@curves)){
    lines(curves_0m[[i]]@curves[[j]], lwd = 4, col = curve_cols[1])
  }
  for(j in 1:length(curves_3m[[i]]@curves)){
    lines(curves_3m[[i]]@curves[[j]], lwd = 4, col = curve_cols[2])
  }
  for(j in 1:length(curves_12m[[i]]@curves)){
    lines(curves_12m[[i]]@curves[[j]], lwd = 4, col = curve_cols[3])
  }

}
dev.off()






### Plot only time points
top_clonotypes <- clonotype_seurat@meta.data %>% group_by(new_clonotypes_id) %>% summarise(n=n()) %>% arrange(desc(n)) %>% head(49) %>% pull(new_clonotypes_id)
unique(top_clonotypes)

plotClonotype <- function(x){
  
  temp_seurat <- subset(clonotype_seurat, new_clonotypes_id == x)
  p <- DimPlot(temp_seurat, group.by = "timepoint", cols = getPalette3(4), pt.size = 1) + labs(x = "", y = "") + theme_bw() + theme(legend.position = "none")
#  print(p)
  
}

q <- lapply(top_clonotypes, plotClonotype)
p <- cowplot::plot_grid(plotlist = q, ncol = 7, nrow = 7)
ggsave(plot = print(p), "results/manuscript/tcr/umap_top_clonotypes.png", width = 12, height = 11)



transition_clonotypes <- c(
  "720_clonotype1",
  "720_clonotype2",
  "730_clonotype3",
  "730_clonotype5",
  "720_clonotype3",
  
  "706_clonotype1",
  "716_clonotype2",
  "730_clonotype6",
  "716_clonotype1",
  "720_clonotype4",
  
  "706_clonotype2",
  "730_clonotype7",
  "720_clonotype5",
  "720_clonotype6",
  "716_clonotype4"
  
  
)


lefttoright_clonotypes <- c(
  "720_clonotype1",
  "720_clonotype2",
  "730_clonotype3",
  "730_clonotype5",
  "720_clonotype3",
  
  "730_clonotype6",
  "720_clonotype4",
  "720_clonotype5",
  "720_clonotype6",
  "730_clonotype8",
  
  "730_clonotype9",
  "730_clonotype10",
  "720_clonotype7",
  "720_clonotype9")

q <- lapply(lefttoright_clonotypes, plotClonotype) 
p <- cowplot::plot_grid(plotlist = q, ncol = 4)
ggsave(plot = print(p), "results/manuscript/tcr/umap_top_clonotypes_left2right.png", width = 8, height = 8)




lefttoright_clonotypes2 <- c(
  "720_clonotype1",
  "730_clonotype5",
  "720_clonotype3",
  
  "720_clonotype4",
  "720_clonotype5",
  "720_clonotype6",

  "730_clonotype9",
  "730_clonotype10",
  "720_clonotype7")

q <- lapply(lefttoright_clonotypes2, plotClonotype) 
p <- cowplot::plot_grid(plotlist = q, ncol = 3)
ggsave(plot = print(p), "results/manuscript/tcr/umap_top_clonotypes_left2right2.png", width = 7, height = 6.5)



DimPlot(temp_seurat, group.by = "timepoint", cols = getPalette3(4), pt.size = 1) + labs(x = "", y = "") + theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave( "results/manuscript/tcr/umap_timepoint.png", width = 7, height = 6.5)


