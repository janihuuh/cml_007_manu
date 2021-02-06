
## Read in the olink data
olink <- fread("data/facs/Olink/007_olink_data.txt")
names <- olink$V1

olink <- t(olink)
colnames(olink) <- make.names(olink[1,])
olink <- olink[-1,]
olink <- gsub("\\,", "\\.", olink)
olink <- apply(olink, 1, FUN = function(x) x %>% as.character %>% as.numeric)

colnames(olink) <- substr(colnames(olink), 5, nchar(colnames(olink)))
colnames(olink) <- make.names(colnames(olink))
rownames(olink) <- make.names(names)

## Remove non-variance cols
olink[is.nan(olink)] <- NA
colnames(olink)[colSums(is.na(olink)) > 0]
olink <- olink[ ,apply(olink, 2, function(x)!all(is.na(x)))]

## Find cols with over 50% of NAs
olink <- olink[,-which(apply(olink, 2, function(x) table(is.na(x))[2] > table(is.na(x))[1]) %>% as.vector())]

## Replace NA with col means
for(i in 1:ncol(olink)){
  olink[is.na(olink[,i]), i] <- mean(olink[,i], na.rm = TRUE)
}

olink <- as.data.frame(olink)
olink$name <- substr(rownames(olink), 2, 5)
olink$timepoint <- substr(rownames(olink), 7, nchar(rownames(olink)))
olink$timepoint[olink$timepoint == "dg"] <- "0mo"
fwrite(olink, "results/manuscript/olink_clean.txt", sep = "\t", quote = F, row.names = T)



#### PCA time


extractOlinkName <- function(str1){
  # strsplit(str1, "[_]")[[1]][2]
  sub("\\..*", "", str1)
}
extractOlinkTimepoint <- function(str1){
  strsplit(str1, "[.]")[[1]][2]
  # sub("*.\\.", "", str1)
}

###### ## PCA time
olink_full      <- fread("results/manuscript/olink_clean.txt")
olink_full_mtx  <- olink_full %>% dplyr::select(-c(V1,name,timepoint))
olink_full_clin <- olink_full %>% dplyr::select(c(V1,name,timepoint))

## Remove rows with over 75% NA
ind1 <- rowSums(is.na(olink_full_mtx)) != ncol(olink_full_mtx) 
ind2 <- apply(olink_full_mtx, 1, function(x)table(is.na(x))[1]) / 36 > 0.75

olink_full_mtx <- olink_full_mtx[intersect(which(ind1),which(ind2)), ]
olink_full_mtx <- olink_full_mtx %>% as.data.frame()

olink_full_clin <- olink_full_clin[intersect(which(ind1),which(ind2)), ]
olink_full_clin <- olink_full_clin %>% as.data.frame()

olink_full_pca <- prcomp((olink_full_mtx), scale. = T)
pca_olink_df <- olink_full_clin %>% bind_cols(as.data.frame(olink_full_pca$x)[,1:10]) %>% bind_cols(olink_full_mtx)

pca_olink_df$timepoint <- factor(as.character(pca_olink_df$timepoint), levels = c("0mo", "3mo", "12mo", "24mo"))
pca_olink_df$timepoint <- plyr::revalue(pca_olink_df$timepoint, replace = c("0mo"  = "0m",
                                                                            "3mo"  = "3m (dasa)",
                                                                            "12mo" = "12m (dasa+ifna)",
                                                                            "24mo" = "24m (dasa)"))
fwrite(pca_olink_df, "results/manuscript/overall/pca_olink_df.txt", sep = "\t", quote = F, row.names = F)
# pca_olink_df <- fread("results/manuscript/overall/pca_olink_df.txt")

ggplot(pca_olink_df, aes(PC1,PC2, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink.pdf", width = 5, height = 5)

ggplot(pca_olink_df, aes(PC3,PC4, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_34.pdf", width = 5, height = 5)

ggplot(pca_olink_df, aes(PC5,PC6, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_56.pdf", width = 5, height = 5)

ggplot(pca_olink_df, aes(PC2,PC3, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_23.pdf", width = 5, height = 5)









###### ## PCA time; no dg
olink_full      <- fread("results/manuscript/olink_clean.txt")
olink_full_mtx  <- olink_full %>% filter(timepoint != "0mo") %>% dplyr::select(-c(V1,name,timepoint))
olink_full_clin <- olink_full %>% filter(timepoint != "0mo") %>% dplyr::select(c(V1,name,timepoint))

## Remove rows with over 75% NA
ind1 <- rowSums(is.na(olink_full_mtx)) != ncol(olink_full_mtx) 
ind2 <- apply(olink_full_mtx, 1, function(x)table(is.na(x))[1]) / 36 > 0.75

olink_full_mtx <- olink_full_mtx[intersect(which(ind1),which(ind2)), ]
olink_full_mtx <- olink_full_mtx %>% as.data.frame()

olink_full_clin <- olink_full_clin[intersect(which(ind1),which(ind2)), ]
olink_full_clin <- olink_full_clin %>% as.data.frame()

olink_full_pca <- prcomp((olink_full_mtx), scale. = T)
pca_olink_df <- olink_full_clin %>% bind_cols(as.data.frame(olink_full_pca$x)[,1:10]) %>% bind_cols(olink_full_mtx)

pca_olink_df$timepoint <- factor(as.character(pca_olink_df$timepoint), levels = c("3mo", "12mo", "24mo"))
pca_olink_df$timepoint <- plyr::revalue(pca_olink_df$timepoint, replace = c("3mo"  = "3m (dasa)",
                                                                            "12mo" = "12m (dasa+ifna)",
                                                                            "24mo" = "24m (dasa)"))
fwrite(pca_olink_df, "results/manuscript/overall/pca_olink_df_no_dg.txt", sep = "\t", quote = F, row.names = F)
# pca_olink_df <- fread("results/manuscript/overall/pca_olink_df.txt")

ggplot(pca_olink_df, aes(PC1,PC2, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_no_dg.pdf", width = 5, height = 5)

ggplot(pca_olink_df, aes(PC3,PC4, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_34_no_dg.pdf", width = 5, height = 5)

ggplot(pca_olink_df, aes(PC5,PC6, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_56_no_dg.pdf", width = 5, height = 5)

ggplot(pca_olink_df, aes(PC2,PC3, fill=timepoint)) + geom_point(shape = 21, size = 3) + scale_fill_manual(values = getPalette3(4)[-1]) + theme_bw(base_size = 12) + labs(fill="") + theme(legend.position = "top")
ggsave("results/manuscript/overall/pca_olink_23_no_dg.pdf", width = 5, height = 5)

