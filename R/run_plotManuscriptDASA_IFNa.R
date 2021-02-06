
dir.create("results/manuscript/ifna_dasa/", showWarnings = F)

## Compare before and after ifna

## FACS
facs_24mo_only <- fread("results/manuscript/facs_24mo_df.txt") 
ifna_dasa_facs_df <- facs_24mo_only %>% melt(id = c("timepoint","Study.nro","FM", "HRUH")) %>% 
  mutate(fill = ifelse(timepoint == "12mo", "dasa+IFNa", "dasa")) %>% 
  filter(timepoint %in% c("12mo", "24mo")) %>% dplyr::select(-Study.nro, -HRUH) %>% dplyr::rename(name = FM) %>% dplyr::select(name,timepoint,variable,value,fill) %>% mutate(type = "facs")
ifna_dasa_facs_df <- ifna_dasa_facs_df %>% filter(!is.na(name))

## Olink
olink_df <- fread("results/manuscript/olink_clean.txt") %>% dplyr::select(-V1) 
olink_ifna_dasa_df <- olink_df %>% melt(id = c("name", "timepoint")) %>% filter(timepoint %in% c("12mo", "24mo")) %>%   mutate(fill = ifelse(timepoint == "12mo", "dasa+IFNa", "dasa")) %>% 
  mutate(fill = ifelse(timepoint == "12mo", "dasa+IFNa", "dasa")) %>% mutate(type = "olink")

ifna_dasa_df <- rbind(olink_ifna_dasa_df, ifna_dasa_facs_df)

p <- ifna_dasa_df %>% ggplot(aes(timepoint, value, fill = fill)) + geom_violin(draw_quantiles = 0.5) + ggpubr::stat_compare_means(label = "p.signif") + facet_wrap(~variable) + # ggsignif::geom_signif(comparisons = list(c("3mo", "12mo"))) +
  geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)) + geom_jitter(size = 0.5)
# ggsave(plot = p, "results/manuscript/ifna_dasa/violin_overall.pdf", width = 12, height = 12)


var_df <- ifna_df %>% group_by(variable, type) %>% summarise(n = n()) %>% dplyr::select(-n)

p.df <- lapply(unique(ifna_dasa_df$variable), FUN = function(x){
  message(x)
  y <- ifna_dasa_df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 2){
    median.x = y %>% filter(timepoint == "12mo") %>% pull(value) %>% median(na.rm = T)
    median.y = y %>% filter(timepoint == "24mo") %>% pull(value) %>% median(na.rm = T)
    wilcox.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x, median.12m = median.x, median.24m = median.y, dir = ifelse(log2(median.y/median.x) < 0, "up", "down"))
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) %>% left_join(var_df)
fwrite(p.df, "results/manuscript/ifna_dasa/ifna_dasa_p_df.txt", sep = "\t", quote = F, row.names = F)
xlsx::write.xlsx(x = p.df, file = "results/manuscript/dasatinib/ifna_dasa_p_df.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)




p.df %>% filter(p.adj < 0.1) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = type)) + geom_bar(stat = "identity") + coord_flip() + labs(x = "") + facet_wrap(~dir, scales = "free_y")
ggsave("results/manuscript/ifna_dasa/bar_sigf.pdf", width = 5, height = 6)


p.df %>% filter(type == "facs") %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
  mutate(col = ifelse(p.adj < 0.05, dir, "no")) %>% 
  mutate(dir = plyr::revalue(dir, replace = c("down" = "up in dasa+ifna", "up" = "up in on dasa"))) %>% 
  ggplot(aes(reorder(variable,-log10(p.adj)), -log10(p.adj), fill = col)) + 
  # geom_bar(stat = "identity") + 
  geom_segment(aes(x = reorder(variable,-log10(p.adj)), xend = reorder(variable,-log10(p.adj)), y=0, yend=-log10(p.adj)), color = "gray30") + geom_point(shape = 21, size = 5) +
  coord_flip() + labs(x = "") + facet_wrap(~dir) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + facets_nice + scale_fill_manual(values = c("lightgrey", getPalette(8)[2], getPalette(8)[1])) + theme(legend.position = "none")
ggsave("results/manuscript/ifna_dasa/loliplot_facs_sigf.pdf", width = 5.5, height = 6)

