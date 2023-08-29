


s_8 <- get_classes(rlsa, k = 8)
colnames(s_8)[1] <- "class_alter"
s_8$sample_id <- rownames(s_8)
s_8$class_alter <- as.factor(s_8$class_alter)
levels(s_8$class_alter) <- c(5, 3, 6, 7, 4, 1, 8, 2)
cbg8 <- left_join(groupings, s_8[, c("sample_id", "class_alter")])
cl8 <- cbg8 %>% group_by(class_alter, subgroups) %>% summarise(value = n()) %>% as.data.frame()

cl8$class_alter <- paste0("G", as.character(cl8$class_alter))

nodes8 <- data.frame(name = c(cl8$class_alter, cl8$subgroups) %>% unique())
edges8 <- cl8
edges8$Source <- match(cl8$class_alter, nodes8$name) - 1
edges8$Target <- match(cl8$subgroups, nodes8$name) - 1

sankeyNetwork(Links = edges8, Nodes = nodes8, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)


grouping_anno <- s_8
grouping_anno$sample_id <- rownames(grouping_anno)
grouping_anno <- grouping_anno[, c("sample_id", "class_alter")]
colnames(grouping_anno) <- c("sample_id", "groups")
grouping_anno$groups <- paste0("G", grouping_anno$groups)
grouping_anno$groups <- factor(grouping_anno$groups, levels = paste0("G", 1:8))


px <- list()
for (i in 1:length(interest_genes)){
  px[[interest_genes[i]]] <- gp_violin(exp_s, interest_genes[i], grouping_anno)
}

line1 <- ggarrange(px[[1]], px[[2]], px[[3]], ncol = 3, nrow = 1)
line1 <- annotate_figure(line1, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
line2 <- ggarrange(px[[4]], px[[5]], ncol = 3, nrow = 1)
line2 <- annotate_figure(line2, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
merged <- ggarrange(line1, line2, ncol = 1, nrow = 2)
ggsave("/cluster/home/yjliu_jh/projects/leu_j/output/new_violin3.pdf", merged, dpi=300,
       width = 12, height = 5, units = "in")




px2 <- list()
for (i in 1:length(interest_genes)){
  px2[[interest_genes[i]]] <- gp_violin2(exp_s, interest_genes[i], grouping_anno, z = 0.3)
}

line1 <- ggarrange(px2[[1]], px2[[2]], px2[[3]], ncol = 3, nrow = 1)
line1 <- annotate_figure(line1, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
line2 <- ggarrange(px2[[4]], px2[[5]], ncol = 3, nrow = 1)
line2 <- annotate_figure(line2, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
merged <- ggarrange(line1, line2, ncol = 1, nrow = 2)
ggsave("/cluster/home/yjliu_jh/projects/leu_j/output/new_violin3.pdf", merged, dpi=300,
       width = 12, height = 8, units = "in")






# get_stats(rlsa, 8)
#s_8 <- get_classes(rlsa["SD:skmeans"], k = 8)



