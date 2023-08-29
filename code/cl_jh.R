

# cr: flora
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

of_fil <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_221216.rds")

rm_genes <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/target_pattern_extracted_genes.xlsx", sheet = 2, col_names = F)
rm_genes <- rm_genes[[1]]


black_path <- fs::path_package("jhdata", "extdata/genesets/heatmap_blackgeneset.xlsx")
black1 <- readxl::read_excel(black_path, sheet = 1)
black2 <- readxl::read_excel(black_path, sheet = 2)
black3 <- readxl::read_excel(black_path, sheet = 3)
black1 <- black1[[1]]
black2 <- black2[[1]]
black3 <- black3[[1]]
blackgenes <- unique(c(black1, black2, black3))

MT_genes <- rownames(ordered_mat)[grepl("MT-", rownames(ordered_mat))]
HLA_genes <- rownames(ordered_mat)[grepl("HLA-", rownames(ordered_mat))]
remove_genes <- setdiff(c(MT_genes, HLA_genes), c("MT-ND6", "HLA-DQA2", "HLA-DQB2"))
blackgenes <- c(blackgenes, remove_genes)


var <- matrixStats::rowVars(of_fil)
fpkm_var <- cbind(of_fil, var)
fpkm_ordered <- as.data.frame(fpkm_var) %>% dplyr::arrange(desc(var))
fpkm_ordered <- fpkm_ordered[rownames(fpkm_ordered) %notin% rm_genes, ]
fpkm_mat <- fpkm_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 


#jhuanglabGO::varied_genes_htmap


allcolour = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", 
              "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", 
              "#FFFF00", "#808000", "#FF00FF", "#FA8072", "#7B68EE", 
              "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", 
              "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", 
              "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", 
              "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", 
              "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", 
              "#DEB887")



sinfo <- readr::read_tsv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/sampleinfo.txt")
batch <- sinfo[, c("sample_id", "dataset")]
batch <- batch[match(colnames(of_fil), batch$sample_id), ]
batch$dataset[is.na(batch$dataset)] <- "Unknown"
batch$sample_id[is.na(batch$sample_id)] <- setdiff(colnames(of_fil), batch$sample_id)
batch$dataset <- as.factor(batch$dataset)

sinfo2 <- left_join(batch, sinfo[, 1:4])
sinfo2[is.na(sinfo2)] <- "Unknown"


hax <- HeatmapAnnotation(disease = sinfo2$disease_type,
                         rna_type = sinfo2$fixed_rna_type,
                         dataset = sinfo2$dataset,
                         gender = sinfo2$fixed_gender,
                        col = list(disease = setNames(allcolour[1:length(unique(sinfo2$disease_type))], 
                                                      unique(sinfo2$disease_type)),
                                   rna_type = setNames(allcolour[1:length(unique(sinfo2$fixed_rna_type))], 
                                                       unique(sinfo2$fixed_rna_type)),
                                   dataset = setNames(allcolour[1:length(unique(sinfo2$dataset))], 
                                                      unique(sinfo2$dataset)), 
                                   gender = setNames(c("blue", "pink", "gray", "lightgray"), 
                                                      unique(sinfo2$fixed_gender))
                        ),
                        simple_anno_size = unit(0.2, "cm"),
                        annotation_name_gp = gpar(fontsize = 7))

hjh <- ComplexHeatmap::Heatmap(fpkm_mat[1:1000, ], 
                              name = " ",
                              show_row_names = FALSE, 
                              use_raster = T,
                              show_column_names = FALSE, 
                              clustering_method_columns = 'ward.D2',
                              clustering_method_rows = 'ward.D2',
                              show_row_dend = F,
                              show_column_dend = F,
                              column_title = "jh_heatmap",
                              heatmap_legend_param = list(
                                legend_direction = "horizontal", 
                                legend_width = unit(3, "cm")),
                              top_annotation = hax)


htjh <- draw(hjh, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh.png", width = 1800, height = 1500, units = "px")
htjh
dev.off()


# filter bad samples

column_order_h <- column_order(htjh)
odsamples <- colnames(fpkm_mat[, column_order_h])
testsamples <- odsamples[(length(column_order_h) - 150):length(column_order_h)]
badsamples <- testsamples[93:151]   ## observation

# 
# badsamples <- c(
#   "A7.x","A61", "A60.x", "A6",
#   "A59.x", "A58", "A57.x", "A56", "A55", "A54.x", "A53.x", "A52.x", "A51",
#   "A50.x", "A5.x", "A49", "A48.x", "A47.x", "A46.x", "A45.x", "A44.x", "A43",
#   "A42", "A41.x", "A40.x", "A4.x", "A39.x", "A38", "A37", "A36.x", "A35.x",
#   "A34.x", "A33.x", "A32.x", "A31", "A30.x", "A3.x", "A29.x", "A28.x", "A27.x",
#   "A26", "A25.x", "A24", "A23", "A22", "A21.x", "A20.x", "A2", "A19.x",
#   "A18.x", "A17.x", "A16.x", "A15", "A14.x", "A13.x", "A12.x", "A11.x", "A1.x", "A10"
# )




# 
# sinfo2f <- sinfo2[sinfo2$sample_id %notin% badsamples, ]
# 
# hax2 <- HeatmapAnnotation(disease = sinfo2f$disease_type,
#                           rna_type = sinfo2f$fixed_rna_type,
#                           dataset = sinfo2f$dataset,
#                           gender = sinfo2f$fixed_gender,
#                           col = list(disease = setNames(allcolour[1:length(unique(sinfo2f$disease_type))], 
#                                                         unique(sinfo2f$disease_type)),
#                                      rna_type = setNames(allcolour[1:length(unique(sinfo2f$fixed_rna_type))], 
#                                                          unique(sinfo2f$fixed_rna_type)),
#                                      dataset = setNames(allcolour[1:length(unique(sinfo2f$dataset))], 
#                                                         unique(sinfo2f$dataset)), 
#                                      gender = setNames(c("blue", "pink", "gray", "lightgray"), 
#                                                        unique(sinfo2f$fixed_gender))
#                           ),
#                           simple_anno_size = unit(0.2, "cm"),
#                           annotation_name_gp = gpar(fontsize = 7))
# 
# 
# hjh2 <- ComplexHeatmap::Heatmap(fpkm_mat[, colnames(fpkm_mat) %notin% badsamples][1:1000, ], 
#                                name = " ",
#                                show_row_names = FALSE, 
#                                use_raster = T,
#                                show_column_names = FALSE, 
#                                clustering_method_columns = 'ward.D2',
#                                clustering_method_rows = 'ward.D2',
#                                show_row_dend = F,
#                                show_column_dend = F,
#                                column_title = "jh_heatmap",
#                                heatmap_legend_param = list(
#                                  legend_direction = "horizontal", 
#                                  legend_width = unit(3, "cm")),
#                                top_annotation = hax2)
# 
# 
# htjh2 <- draw(hjh2, heatmap_legend_side="right")
# 
# png("/cluster/home/yjliu_jh/share/cluster_jh_filterbad.png", width = 1800, height = 1500, units = "px")
# htjh2
# dev.off()




# ======== get filtered and ordered matrix ========

hjh3 <- ComplexHeatmap::Heatmap(fpkm_mat[rownames(fpkm_mat) %notin% blackgenes,
                                         colnames(fpkm_mat) %notin% badsamples][1:1000, ], 
                                name = " ",
                                show_row_names = FALSE, 
                                use_raster = T,
                                show_column_names = FALSE, 
                                clustering_method_columns = 'ward.D2',
                                clustering_method_rows = 'ward.D2',
                                show_row_dend = F,
                                show_column_dend = F,
                                column_title = "jh_heatmap",
                                heatmap_legend_param = list(
                                  legend_direction = "horizontal", 
                                  legend_width = unit(3, "cm")),
                                top_annotation = hax2)


htjh3 <- draw(hjh3, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_filter_final.png", width = 1800, height = 1500, units = "px")
htjh3
dev.off()



ori_mat <- fpkm_mat[rownames(fpkm_mat) %notin% blackgenes, colnames(fpkm_mat) %notin% badsamples][1:1000, ]
column_order_h3 <- column_order(htjh3)
row_order_h3 <- row_order(htjh3)
ordered_mat <- ori_mat[row_order_h3, column_order_h3]
readr::write_rds(ordered_mat, "/cluster/home/yjliu_jh/share/ordered_leu_mat_1217.rds")

#odsamples <- colnames(fpkm_mat[, column_order_h])

ht_info <- data.frame(sample_id = colnames(ordered_mat))
ht_info <- left_join(ht_info, sinfo2)
which(ht_info$disease_type %in% "cll") # check cll clusters






# 
# 
# sub_mat <- ordered_mat[1:500, 7000:7500]
# 
# sub_sinfo <- sinfo2[sinfo2$sample_id %in% colnames(sub_mat), ]
# 
# hasub <- HeatmapAnnotation(disease = sub_sinfo$disease_type,
#                           rna_type = sub_sinfo$fixed_rna_type,
#                           dataset = sub_sinfo$dataset,
#                           gender = sub_sinfo$fixed_gender,
#                           col = list(disease = setNames(allcolour[1:length(unique(sub_sinfo$disease_type))], 
#                                                         unique(sub_sinfo$disease_type)),
#                                      rna_type = setNames(allcolour[1:length(unique(sub_sinfo$fixed_rna_type))], 
#                                                          unique(sub_sinfo$fixed_rna_type)),
#                                      dataset = setNames(allcolour[1:length(unique(sub_sinfo$dataset))], 
#                                                         unique(sub_sinfo$dataset)), 
#                                      gender = setNames(c("blue", "pink", "gray"), 
#                                                        unique(sub_sinfo$fixed_gender))
#                           ),
#                           simple_anno_size = unit(0.2, "cm"),
#                           annotation_name_gp = gpar(fontsize = 7))
# 
# 
# 
# hjhsub2 <- ComplexHeatmap::Heatmap(sub_mat, 
#                                 name = " ",
#                                 show_row_names = T, 
#                                 use_raster = T,
#                                 show_column_names = FALSE,
#                                 cluster_rows = F,
#                                 cluster_columns = F,
#                                 column_title = "jh_heatmap",
#                                 heatmap_legend_param = list(
#                                   legend_direction = "horizontal", 
#                                   legend_width = unit(3, "cm")),
#                                 top_annotation = hasub)
# 
# 
# htjhsub2 <- draw(hjhsub2, heatmap_legend_side="right")
# 
# png("/cluster/home/yjliu_jh/share/cluster_jh_sub4.png", width = 1800, height = 1500, units = "px")
# htjhsub2
# dev.off()
# 



# ==== filter cll and heatmap and check overlap of aml-like samples ===

mat_nocll <- ordered_mat[, setdiff(1:ncol(ordered_mat), c(7173:7448))]

aml_like_samples <- colnames(ordered_mat)[1:2755][ht_info[1:2755, ]$disease_type %notin% "aml"]


nocll_sinfo <- sinfo2[sinfo2$sample_id %in% colnames(mat_nocll), ]
nocll_sinfo <- nocll_sinfo[match(colnames(mat_nocll), nocll_sinfo$sample_id), ]


hanocll <- HeatmapAnnotation(disease = nocll_sinfo$disease_type,
                           rna_type = nocll_sinfo$fixed_rna_type,
                           dataset = nocll_sinfo$dataset,
                           gender = nocll_sinfo$fixed_gender,
                           col = list(disease = setNames(allcolour[1:length(unique(nocll_sinfo$disease_type))], 
                                                         unique(nocll_sinfo$disease_type)),
                                      rna_type = setNames(allcolour[1:length(unique(nocll_sinfo$fixed_rna_type))], 
                                                          unique(nocll_sinfo$fixed_rna_type)),
                                      dataset = setNames(allcolour[1:length(unique(nocll_sinfo$dataset))], 
                                                         unique(nocll_sinfo$dataset)), 
                                      gender = setNames(c("blue", "pink", "gray", "lightgray"), 
                                                        unique(nocll_sinfo$fixed_gender))
                           ),
                           simple_anno_size = unit(0.2, "cm"),
                           annotation_name_gp = gpar(fontsize = 7))


hjhsub4 <- ComplexHeatmap::Heatmap(mat_nocll[intersect(rownames(fpkm_mat)[1:800], rownames(mat_nocll)), ], 
                                   name = " ",
                                   show_row_names = F, 
                                   use_raster = T,
                                   show_column_names = FALSE,
                                   clustering_method_columns = 'ward.D2',
                                   clustering_method_rows = 'ward.D2',
                                   show_row_dend = F,
                                   show_column_dend = F,
                                   column_title = "jh_heatmap",
                                   heatmap_legend_param = list(
                                     legend_direction = "horizontal", 
                                     legend_width = unit(3, "cm")),
                                   top_annotation = hanocll)


htjhsub4 <- draw(hjhsub4, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_bigsub_800.png", width = 1800, height = 1500, units = "px")
htjhsub4
dev.off()


# ===== get another matrix for comparison ===

ori_mat2 <- mat_nocll[intersect(rownames(fpkm_mat)[1:800], rownames(mat_nocll)), ]
column_order_h4 <- column_order(htjhsub4)
row_order_h4 <- row_order(htjhsub4)
ordered_mat2 <- ori_mat2[row_order_h4, column_order_h4]


c_dend <- column_dend(htjh3)
col_dend <- as.hclust(c_dend)
ctree <- cutree(col_dend, k = 5)
ctree[names(ctree) %in% c("BA2381R", "BA3218R", "BA2341R", "BA2612R", "BA2973R")]  ## 2

clu_info <- data.frame(sample_id = names(ctree), cluster = ctree)
clu_info <- left_join(clu_info, sinfo2)
table(clu_info[clu_info$cluster == 2, ]$disease_type)
aml_like1 <- clu_info[clu_info$cluster == 2 & clu_info$disease_type != "aml", ]$sample_id

c_dend2 <- column_dend(htjhsub4)
col_dend2 <- as.hclust(c_dend2)
ctree2 <- cutree(col_dend2, k = 5)
ctree2[names(ctree2) %in% c("BA2381R", "BA3218R", "BA2341R", "BA2612R", "BA2973R")]  ## 4
#table(ctree2)

clu_info2 <- data.frame(sample_id = names(ctree2), cluster = ctree2)
clu_info2 <- left_join(clu_info2, sinfo2)
table(clu_info2[clu_info2$cluster == 4, ]$disease_type)
aml_like2 <- clu_info2[clu_info2$cluster == 4 & clu_info2$disease_type != "aml", ]$sample_id



robust_aml_like <- intersect(aml_like1, aml_like2)
readr::write_rds(robust_aml_like, "/cluster/home/yjliu_jh/share/robust_aml_like_samples_1217.rds")

































filter_black_genes <- function(dat, geneset = NULL, black_gene_list = T){
  stopifnot("please make sure gene_name is contained in the data" = "gene_name" %in% colnames(dat))
  filter_gene_list <- NULL
  if(black_gene_list){
    #read gene list
    black_path <- fs::path_package("jhdata", "extdata/genesets/heatmap_blackgeneset.xlsx")
    all_black_sheets <- readxl::excel_sheets(black_path)
    black_genes <- lapply(all_black_sheets, function(X) readxl::read_excel(black_path, sheet = X))
    blackgenes <- unique(bind_rows(black_genes)$gene_name)
    filter_gene_list <- c(blackgenes, geneset)
  }
  #filter filter_gene_list in dat
  dat_fil <- dat[dat$gene_name %notin% filter_gene_list, ]
  return(dat_fil)
}


























test_mat <- fpkm_mat[rownames(fpkm_mat) %notin% blackgenes,
                     colnames(fpkm_mat) %notin% badsamples][intersect(cl_genes, rownames(fpkm_mat))]






