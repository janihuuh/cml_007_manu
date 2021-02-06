
clin_df <- fread("data/facs/007_0mo.txt") 
colnames(clin_df) <- make.unique(make.names(colnames(clin_df)))
clin_df <- clin_df %>% select(FM, Patient.nro, Sokal, Sokal.1)

mmr_df  <- fread("data/clinical/mmr.csv")
mmr_df <- sapply(mmr_df, as.character)
mmr_df <- gsub("\\,", "\\.", mmr_df)
mmr_df <- apply(mmr_df, 2, as.numeric)
colnames(mmr_df) <- make.unique(make.names(colnames(mmr_df)))
mmr_df <- as.data.frame(mmr_df)

clin_df <- clin_df %>% left_join(mmr_df, by = c("Patient.nro" = "Pat.number"))
clin_df$mr4 <- ifelse(clin_df$M18 < 0.01, "MR4", "other")
clin_df$mr5 <- ifelse(clin_df$M18 < 0.0032, "MR5", "other")
