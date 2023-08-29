# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# read data
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
meta_leu <- metadata(leu)
sampleinfo_r <- as.data.frame(colData(leu))


#assayNames(leu)

leu_tpm_raw <- assay(leu, "tpm_raw")
leu_rm_rna <- assay(leu, "tpm_remove_rnatype")
leu_rm_cohort <- assay(leu, "tpm_remove_cohort")
leu_tpm <- assay(leu, "tpm")


# plot heatmap and add disease type annotations

black_path <- fs::path_package("jhdata", "extdata/genesets/heatmap_blackgeneset.xlsx")
all_black_sheets <- readxl::excel_sheets(black_path)
black_genes <- lapply(all_black_sheets, function(X) readxl::read_excel(black_path, sheet = X))
blackgenes <- unique(bind_rows(black_genes)$gene_name)




# heatmap function 
quick_heat_n <- function(dat, anno){
  dat <- dat[rownames(dat) %notin% blackgenes, ]
  var <- matrixStats::rowVars(as.matrix(dat))
  dat_var <- cbind(dat, var)
  dat_ordered <- as.data.frame(dat_var) %>% dplyr::arrange(desc(var))
  mat <- dat_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
    as.matrix() 
  h <- ComplexHeatmap::Heatmap(mat[1:1000, ], 
                               name = " ",
                               show_row_names = F, 
                               use_raster = F,
                               show_column_names = FALSE, 
                               clustering_method_columns = 'ward.D2',
                               clustering_method_rows = 'ward.D2',
                               show_row_dend = F,
                               show_column_dend = T,
                               show_heatmap_legend = F,
                               column_dend_height = unit(6, "cm"),
                               column_title = "jh_heatmap",
                               #heatmap_legend_param = list(
                               #  legend_direction = "horizontal", 
                               #  legend_width = unit(3, "cm")),
                               top_annotation = anno)
  return(h)
}




# set jh color palette
allcolour = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", 
              "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", 
              "#FFFF00", "#808000", "#FF00FF", "#FA8072", "#7B68EE", 
              "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", 
              "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", 
              "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", 
              "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", 
              "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", 
              "#DEB887")

ha_fil <- HeatmapAnnotation(disease = sampleinfo_r$disease_type,
                            col = list(disease = setNames(allcolour[1:length(na.omit(unique(sampleinfo_r$disease_type)))], 
                                                          na.omit(unique(sampleinfo_r$disease_type)))),
                            simple_anno_size = unit(0.4, "cm"),
                            annotation_name_gp = gpar(fontsize = 9))
#identical(sampleinfo_r$sample_id, colnames(leu_tpm))   ## make sure the orders are right



h_fil <- quick_heat_n(leu_tpm_raw, ha_fil)
ht_fil <- draw(h_fil, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/dtype_tpm_raw.png", width = 2400, height = 1200, units = "px")
ht_fil
dev.off()


h_fil1 <- quick_heat_n(leu_rm_cohort, ha_fil)
ht_fil1 <- draw(h_fil1, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/dtype_rm_cohort.png", width = 2400, height = 1200, units = "px")
ht_fil1
dev.off()


h_fil2 <- quick_heat_n(leu_rm_rna, ha_fil)
ht_fil2 <- draw(h_fil2, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/dtype_rm_rna_type.png", width = 2400, height = 1200, units = "px")
ht_fil2
dev.off()


h_fil3 <- quick_heat_n(leu_tpm, ha_fil)
ht_fil3 <- draw(h_fil3, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/dtype_tpm.png", width = 2400, height = 1200, units = "px")
ht_fil3
dev.off()




# get that small cluster of samples for pca plot

quick_heat_s <- function(dat, anno){
  h <- ComplexHeatmap::Heatmap(dat, 
                               name = " ",
                               show_row_names = F, 
                               use_raster = F,
                               show_column_names = FALSE, 
                               cluster_rows = F,
                               cluster_columns = F,
                               column_title = "jh_heatmap",
                               heatmap_legend_param = list(
                                 legend_direction = "horizontal", 
                                 legend_width = unit(3, "cm")),
                               top_annotation = anno)
  return(h)
}


column_order_h3 <- column_order(ht_fil)
row_order_h3 <- row_order(ht_fil)
ordered_fil_mat <- tpm_raw_1000[row_order_h3, column_order_h3]
sub_mat <- ordered_fil_mat[, 1:200]

htinfo_sub3 <- sampleinfo_r[sampleinfo_r$sample_id %in% colnames(sub_mat), ]
htinfo_sub3 <- htinfo_sub3[match(colnames(sub_mat), htinfo_sub3$sample_id), ]

ha_fil_sub3 <- HeatmapAnnotation(disease = htinfo_sub3$disease_type,
                                 col = list(disease = setNames(allcolour[1:length(na.omit(unique(htinfo_sub3$disease_type)))], 
                                                               na.omit(unique(htinfo_sub3$disease_type)))),
                                 simple_anno_size = unit(0.4, "cm"),
                                 annotation_name_gp = gpar(fontsize = 9))

h_sub3 <- quick_heat_s(sub_mat, ha_fil_sub3)

ht_sub3 <- draw(h_sub3, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/smallcluster_test.png", width = 1200, height = 900, units = "px")
ht_sub3
dev.off()

xsamples <- htinfo_sub3$sample_id[1 : (which(htinfo_sub3$disease_type == "aml")[1] - 1)]
meta_leu <- metadata(leu)
pca_leu <- meta_leu[["pca"]]
sampleinfo_r2 <- sampleinfo_r[match(sampleinfo_r$sample_id, rownames(pca_leu)), ]
pca_leu$sample_id <- as.factor(sampleinfo_r2$sample_id)
pca_leu$colorgroup <- ifelse(pca_leu$sample_id %in% xsamples, "specical", "others")


temp_col <- c("red", "gray")
names(temp_col) <- c("specical", "others")

  pdf("/cluster/home/yjliu_jh/share/leu_clusterx.pdf", width = 7, height = 8)
  print(ggscatter(pca_leu, x = colnames(pca_leu)[1], y = colnames(pca_leu)[2], color = "colorgroup",
                  palette = temp_col, size = 0.2))
  dev.off()




  
  
  
# ===== re-draw the heatmaps =======
  
sampleinfo_rn <- sampleinfo_r[sampleinfo_r$remove %in% "no" & sampleinfo_r$disease_type %notin% "cll", ]

leu_tpm_raw <- leu_tpm_raw[, sampleinfo_rn$sample_id]

ha_fil2 <- HeatmapAnnotation(disease = sampleinfo_rn$disease_type,
                              col = list(disease = setNames(allcolour[1:length(na.omit(unique(sampleinfo_rn$disease_type)))], 
                                                            na.omit(unique(sampleinfo_rn$disease_type)))),
                              simple_anno_size = unit(0.4, "cm"),
                              annotation_name_gp = gpar(fontsize = 9))
  
h_filn <- quick_heat_n(leu_tpm_raw, ha_fil2)
ht_filn <- draw(h_filn, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/dtype_tpm_raw_removed.png", width = 2400, height = 1200, units = "px")
ht_filn
dev.off()


leu_rm_rna <- leu_rm_rna[, sampleinfo_rn$sample_id]
h_filn2 <- quick_heat_n(leu_rm_rna, ha_fil2)
ht_filn2 <- draw(h_filn2, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/dtype_rm_rna_removed.png", width = 2400, height = 1200, units = "px")
ht_filn2
dev.off()

readr::write_rds(ht_filn2, "/cluster/home/yjliu_jh/share/heatmap_rna_remove.rds")




# update the coldata information according to the above clustering
ht_rna <- readr::read_rds("/cluster/home/yjliu_jh/share/heatmap_rna_remove.rds")
column_order_ht <- column_order(ht_rna)
row_order_ht <- row_order(ht_rna)
data_mat <- ht_rna@ht_list$` `@matrix
data_ordered <- data_mat[row_order_ht, column_order_ht]


# read data 
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
col_leu <- as.data.frame(colData(leu))
colnames(colData(leu))[colnames(colData(leu)) %in% "disease_type"] <- "old_disease_type"
col_leu$disease_type <- colData(leu)$old_disease_type

c_dend <- column_dend(ht_rna)
col_dend <- as.hclust(c_dend)
ctree <- cutree(col_dend, k = 3)
col_leu$disease_type[col_leu$sample_id %in% names(ctree[ctree == 1])] <- "ball"
col_leu$disease_type[col_leu$sample_id %in% names(ctree[ctree == 2])] <- "aml"
col_leu$disease_type[col_leu$sample_id %in% names(ctree[ctree == 3])] <- "tall"
colData(leu)$disease_type <- col_leu$disease_type

readr::write_rds(leu, rds_fn)



# ======= seaparated PCA and heatmap plots

rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
meta_leu <- metadata(leu)
sampleinfo_r <- as.data.frame(colData(leu))
fusion_info <- meta_leu$fusion_info$fusion_groups
leu_rm_rna <- assay(leu, "tpm_remove_rnatype")

data_ordered_all <- leu_rm_rna[, colnames(data_ordered)]
dat_ball <- data_ordered_all[, 1:table(ctree)[1]]
dat_aml <- data_ordered_all[, (table(ctree)[1] + 1) : (table(ctree)[1] + table(ctree)[2])]
dat_tall <- data_ordered_all[, (ncol(data_ordered_all) - table(ctree)[3] + 1) : ncol(data_ordered_all)]

#sinfo_all <- left_join(sampleinfo_r, fusion_info)


fusion_info2 <- fusion_info
for(i in 2:ncol(fusion_info2)){
  fusion_info2[[colnames(fusion_info2)[i]]] <- as.character(fusion_info2[[colnames(fusion_info2)[i]]])
}


sinfo_all <- left_join(sampleinfo_r, fusion_info2)

sinfo_ball <- sinfo_all[sinfo_all$sample_id %in% colnames(dat_ball), ]
dat_ball <- dat_ball[, sinfo_ball$sample_id]
sinfo_aml <- sinfo_all[sinfo_all$sample_id %in% colnames(dat_aml), ]
dat_aml <- dat_aml[, sinfo_aml$sample_id]
sinfo_tall <- sinfo_all[sinfo_all$sample_id %in% colnames(dat_tall), ]
dat_tall <- dat_tall[, sinfo_tall$sample_id]

fusion_col <- circlize::colorRamp2(c(0, 1), c("white", "black"))
colxxx <- list()
for(i in 2:ncol(fusion_info)){
  colxxx[[colnames(fusion_info)[i]]] <- fusion_col
}

fusion_col2 <- c("white", "black")
names(fusion_col2) <- c("0", "1")
colxxx2 <- list()
for(i in 2:ncol(fusion_info)){
  colxxx2[[colnames(fusion_info)[i]]] <- fusion_col2
}








ha_fil_ball <- HeatmapAnnotation(df = sinfo_ball[, colnames(fusion_info)][, -1],
                            col = colxxx,
                            show_legend = F,
                            simple_anno_size = unit(0.4, "cm"),
                            annotation_name_gp = gpar(fontsize = 9))

ha_fil_tall <- HeatmapAnnotation(df = sinfo_tall[, colnames(fusion_info)][, -1],
                                 col = colxxx,
                                 show_legend = F,
                                 simple_anno_size = unit(0.4, "cm"),
                                 annotation_name_gp = gpar(fontsize = 9))

ha_fil_aml <- HeatmapAnnotation(df = sinfo_aml[, colnames(fusion_info)][, -1],
                                 col = colxxx,
                                show_legend = F,
                                 simple_anno_size = unit(0.4, "cm"),
                                 annotation_name_gp = gpar(fontsize = 9))




h_filna <- quick_heat_n(dat_ball, ha_fil_ball)
ht_filna <- draw(h_filna, show_heatmap_legend = FALSE)

png("/cluster/home/yjliu_jh/share/ball_only_heatmap.png", width = 2400, height = 1200, units = "px")
ht_filna
dev.off()

h_filna2 <- quick_heat_n(dat_tall, ha_fil_tall)
ht_filna2 <- draw(h_filna2, show_heatmap_legend = FALSE)

png("/cluster/home/yjliu_jh/share/tall_only_heatmap.png", width = 2400, height = 1200, units = "px")
ht_filna2
dev.off()

h_filna3 <- quick_heat_n(dat_aml, ha_fil_aml)
ht_filna3 <- draw(h_filna3, show_heatmap_legend = FALSE)

png("/cluster/home/yjliu_jh/share/aml_only_heatmap.png", width = 2400, height = 1200, units = "px")
ht_filna3
dev.off()


# ===== pca ========

mat_ball <- ht_filna@ht_list$` `@matrix
pca_ball <- prcomp(mat_ball, scale. = TRUE)
pca_ball1 <- pca_ball$rotation[, 1:3]
pca_ball1 <- as.data.frame(pca_ball1[sinfo_ball$sample_id, ])
pca_ball1$dataset <- sinfo_ball$dataset
# p_ball_temp <- autoplot(pca_ball, data = sinfo_ball, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/ball_cluster_pca_dataset_1000_new.pdf", width = 6, height = 6)
print(ggscatter(pca_ball1, x = colnames(pca_ball1)[1], y = colnames(pca_ball1)[2], color = "dataset", size = 0.2))
dev.off()

mat_tall <- ht_filna2@ht_list$` `@matrix
pca_tall <- prcomp(mat_tall, scale. = TRUE)
pca_tall1 <- pca_tall$rotation[, 1:3]
pca_tall1 <- as.data.frame(pca_tall1[sinfo_tall$sample_id, ])
pca_tall1$dataset <- sinfo_tall$dataset
# p_ball_temp <- autoplot(pca_ball, data = sinfo_ball, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/tall_cluster_pca_dataset_1000_new.pdf", width = 6, height = 6)
print(ggscatter(pca_tall1, x = colnames(pca_tall1)[1], y = colnames(pca_tall1)[2], color = "dataset", size = 0.2))
dev.off()

mat_aml <- ht_filna3@ht_list$` `@matrix
pca_aml <- prcomp(mat_aml, scale. = TRUE)
pca_aml1 <- pca_aml$rotation[, 1:3]
pca_aml1 <- as.data.frame(pca_aml1[sinfo_aml$sample_id, ])
pca_aml1$dataset <- sinfo_aml$dataset
# p_ball_temp <- autoplot(pca_ball, data = sinfo_ball, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/aml_cluster_pca_dataset_1000_new.pdf", width = 6, height = 6)
print(ggscatter(pca_aml1, x = colnames(pca_aml1)[1], y = colnames(pca_aml1)[2], color = "dataset", size = 0.2))
dev.off()










meta_leu <- metadata(leu)
rdim1 <- meta_leu[[1]]
sampleinfo_r2 <- sampleinfo_r[match(sampleinfo_r$sample_id, rownames(rdim1)), ]
dtype_col <- brewer.pal(n = 12, name = "Paired")
names(dtype_col) <- unique(sampleinfo_r2$disease_type)
#dtype_col <- dtype_col[!is.na(names(dtype_col))]
#dtype2 <- left_join(sampleinfo_r2, data.frame(disease_type = names(dtype_col), color = dtype_col))$color

# Note that print() is needed inside a loop!!
for (rdim in 1:length(meta_leu)){
  rdimx <- as.data.frame(meta_leu[[names(meta_leu)[rdim]]])
  rdimx$disease_type <- as.factor(sampleinfo_r2$disease_type)
  pdf(glue::glue("/cluster/home/yjliu_jh/share/{names(meta_leu)[rdim]}_leu_disease_typex.pdf"), width = 7, height = 8)
  print(ggscatter(rdimx, x = colnames(rdimx)[1], y = colnames(rdimx)[2], color = "disease_type",
                  palette = dtype_col, size = 0.2))
  dev.off()
}
















# ==== corr calculation ===

dat_manu <- function(dat){
  dat <- dat[rownames(dat) %notin% blackgenes, ]
  var <- matrixStats::rowVars(as.matrix(dat))
  dat_var <- cbind(dat, var)
  dat_ordered <- as.data.frame(dat_var) %>% dplyr::arrange(desc(var))
  mat <- dat_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
    as.matrix() 
  return(mat[1:1000, ])
}

tpm_raw_1000 <- dat_manu(leu_tpm_raw)


library(Hmisc)
correlations <- rcorr(tpm_raw_1000)
readr::write_rds(correlations, "/cluster/home/yjliu_jh/share/corr_tpm_raw_1000.rds")

#library(corrgram)
#corrgram(correlations)








# 
# all to aml: common       aml to all: uncommon!



