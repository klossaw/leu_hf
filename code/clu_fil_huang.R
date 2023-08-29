
# check clustering_small_temp.R


#count_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human_counts.csv"))
exp_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human.csv"))
s_fil <- exp_s[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
s_fil <- s_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
s_fil <- round(as.matrix(s_fil))

# use heatgenes from Prof. Huang
remaining_samples <- sinfo2[sinfo2$disease_type %in% "AML", ]$sample_id
s_filh <- s_fil[, colnames(s_fil) %in% remaining_samples]
s_filh1 <- s_filh[rownames(s_fil) %in% heatgenes, ]


# use most variable genes 
osd <- transform(s_filh, SD = apply(s_filh, 1, sd, na.rm = TRUE))
s_fil_15 <- s_filh[order(osd$SD, decreasing = T), ][1:round(0.15 * nrow(osd)), ]
s_fil_10 <- s_filh[order(osd$SD, decreasing = T), ][1:round(0.1 * nrow(osd)), ]
s_fil_5 <- s_filh[order(osd$SD, decreasing = T), ][1:round(0.05 * nrow(osd)), ]

#length(setdiff(heatgenes, rownames(s_fil_10)))
s_filh2 <- s_fil_10


rls_n1 = run_all_consensus_partition_methods(s_filh1, max_k = 8,
                                             top_value_method = c("SD", "MAD", "CV", "ATC"),
                                             partition_method = c("hclust", "kmeans", "skmeans", "mclust", "pam"),
                                             cores = 20)
rls_n2 = run_all_consensus_partition_methods(s_filh2, max_k = 8,
                                             top_value_method = c("SD", "MAD", "CV", "ATC"),
                                             partition_method = c("hclust", "kmeans", "skmeans", "mclust", "pam"),
                                             cores = 20)


readr::write_rds(rls_n1, "/cluster/home/yjliu_jh/projects/leu_j/data/clustering_8_heatgenes.rds")
readr::write_rds(rls_n2, "/cluster/home/yjliu_jh/projects/leu_j/data/clustering_8_10perc_genes.rds")



s8_n1 <- get_classes(rls_n1, k = 8)
s8_n2 <- get_classes(rls_n2, k = 8)
colnames(s8_n1)[1] <- "class_alter1"
colnames(s8_n2)[1] <- "class_alter2"
s8_n1$sample_id <- rownames(s8_n1)
s8_n2$sample_id <- rownames(s8_n2)
s8_n1$class_alter1 <- paste0("1_G", s8_n1$class_alter1)
s8_n2$class_alter2 <- paste0("2_G", s8_n2$class_alter2)




groupings_fil <- groupings[groupings$sample_id %in% remaining_samples, ]
cbgx <- left_join(groupings_fil, s8_n1[, c("sample_id", "class_alter1")])
cbgx <- left_join(cbgx, s8_n2[, c("sample_id", "class_alter2")])


cl8x <- cbgx %>% group_by(subgroups, class_alter1) %>% summarise(value = n()) %>% as.data.frame()
cl8x2 <- cbgx %>% group_by(class_alter1, class_alter2) %>% summarise(value = n()) %>% as.data.frame()



nodesx <- data.frame(name = c(cbgx$subgroups, cbgx$class_alter1, cbgx$class_alter2) %>% unique())
edgesx <- rbind(cl8x, setNames(cl8x2, names(cl8x)))
edgesx$Source <- match(edgesx$class_alter1, nodesx$name) - 1
edgesx$Target <- match(edgesx$subgroups, nodesx$name) - 1

sankeyNetwork(Links = edgesx, Nodes = nodesx, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)








anno_s <- as.data.frame(rarsamples[, c(1,5)])
anno_s <- unique(anno_s[anno_s$sample_id %in% colnames(small_exp),])
fusions <- anno_s$fusions

ha_s <- HeatmapAnnotation(
  df = as.data.frame(fusions),
  gp = gpar(border = "gray", col = "white", lwd = 0.2),
  height = unit(0.4, "cm"), simple_anno_size_adjust = TRUE,
  col = list(
    fusions = fusions_col
  )
)


hm_s <- heatmap.2(small_expm, hclust = function(x) hclust(x, method = "ward.D"), 
                  distfun = function(x) as.dist((1 - cor(t(x)))/2), scale = "row", 
                  key = TRUE, symkey = FALSE, density.info = "none", trace = "none", 
                  cexRow = 0.5)
# use hclust(ward.D) to cluster the heatmap

mat_scaled_s <- t(apply(small_expm, 1, scale))
mat_scaled_s <- t(hm_s$carpet)[row.names(small_expm),colnames(small_expm)]

ht_s2 <- Heatmap(mat_scaled_s,
                 name = "expression",
                 #col = GenerateHeatmapColor("bluered",5),
                 col = circlize::colorRamp2(c(-1.8, 0, 1.8), c("blue", "white", "red")),
                 cluster_columns = F,
                 cluster_rows = rev(hm_s$rowDendrogram),
                 show_row_names = T,
                 show_column_names = F,
                 top_annotation = ha_s,
                 column_title = NULL
)


ht_s <- draw(ht_s, merge_legend = TRUE)
ht_s2 <- draw(ht_s2, merge_legend = TRUE)



ggsave("/cluster/home/yjliu_jh/projects/rarg/output/figs/eheatmap.pdf",
       ht_s, dpi=300, width = 5, height = 5, units = "in")

