# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment", "ComplexHeatmap", "RColorBrewer")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# read data
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
meta_leu <- metadata(leu)
sampleinfo_r <- as.data.frame(colData(leu))




# construct a function to draw sub-heatmaps

sub_heat <- function(ht, ind, sinfo, cname = F, rname = F){
  
  column_order <- column_order(ht)
  row_order <- row_order(ht)
  data_mat <- ht@ht_list$` `@matrix
  ordered_fil_mat <- data_mat[row_order, column_order]
  sub_mat <- ordered_fil_mat[, ind]
  
  htinfo_sub <- sinfo[sinfo$sample_id %in% colnames(sub_mat), ]
  htinfo_sub <- htinfo_sub[match(colnames(sub_mat), htinfo_sub$sample_id), ]
  
  ha_fil_sub <- HeatmapAnnotation(dataset = htinfo_sub$dataset,
                                  col = list(dataset = setNames(allcolour[1:length(na.omit(unique(htinfo_sub$dataset)))], 
                                                                na.omit(unique(htinfo_sub$dataset)))),
                                  simple_anno_size = unit(0.4, "cm"),
                                  annotation_name_gp = gpar(fontsize = 9))
  
  h <- ComplexHeatmap::Heatmap(sub_mat, 
                               name = " ",
                               show_row_names = rname, 
                               use_raster = F,
                               show_column_names = cname, 
                               cluster_rows = F,
                               cluster_columns = F,
                               column_title = "jh_heatmap",
                               heatmap_legend_param = list(
                                 legend_direction = "horizontal", 
                                 legend_width = unit(3, "cm")),
                               top_annotation = ha_fil_sub)
  return(h)
}


get_info <- function(ht, ind, sinfo) {
  column_order <- column_order(ht)
  row_order <- row_order(ht)
  data_mat <- ht@ht_list$` `@matrix
  ordered_fil_mat <- data_mat[row_order, column_order]
  sub_mat <- ordered_fil_mat[, ind]
  
  htinfo_sub <- sinfo[sinfo$sample_id %in% colnames(sub_mat), ]
  htinfo_sub <- htinfo_sub[match(colnames(sub_mat), htinfo_sub$sample_id), ]
  return(htinfo_sub)
}



# ht_filna ball   ht_filna2 tall    ht_filna3 aml
# aml ball  cll tall
# 3369 4557  282  976



sp_aml1 <- sub_heat(ht_filna3, c(570:600, 3180:3368), sampleinfo_r)

htsp_aml1 <- draw(sp_aml1, show_heatmap_legend = FALSE)
png("/cluster/home/yjliu_jh/share/test_aml_sub1.png", width = 2400, height = 1200, units = "px")
htsp_aml1
dev.off()

info_aml1 <- get_info(ht_filna3, c(572:600, 3180:3368), sampleinfo_r)

sp_aml_samples <- info_aml1[info_aml1$dataset %in% c("beataml", "CRA003240", "EGAS00001003266"), ]$sample_id
sp_likely_samples <- info_aml1[info_aml1$dataset %in% c("GSE154790"), ]$sample_id


correlations$r[rownames(correlations$r) %in% sp_aml_samples, colnames(correlations$r) %in% sp_aml_samples]




sp_ball1 <- sub_heat(ht_filna, c(3900:4400), sampleinfo_r)

htsp_ball1 <- draw(sp_ball1, show_heatmap_legend = FALSE)
png("/cluster/home/yjliu_jh/share/test_ball_sub1.png", width = 2400, height = 1200, units = "px")
htsp_ball1
dev.off()

info_ball1 <- get_info(ht_filna, c(4000:4400), sampleinfo_r)
sp_ball_samples <- info_ball1[info_ball1$dataset %in% c("GSE148520"), ]$sample_id



sp_tall1 <- sub_heat(ht_filna2, c(150:250), sampleinfo_r)

htsp_tall1 <- draw(sp_tall1, show_heatmap_legend = FALSE)
png("/cluster/home/yjliu_jh/share/test_tall_sub1.png", width = 2400, height = 1200, units = "px")
htsp_tall1
dev.off()

info_tall1 <- get_info(ht_filna2, c(150:250), sampleinfo_r)
sp_tall_samples <- info_tall1[info_tall1$dataset %in% c("suosansan"), ]$sample_id



sampleinfo_x <- sampleinfo_r
sampleinfo_x$anno <- "others"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% sp_aml_samples] <- "aml_likely_dupe"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% sp_likely_samples] <- "aml_similar"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% sp_ball_samples] <- "ball_likely_dupe"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% sp_tall_samples] <- "tall_likely_dupe"


pca_leu <- meta_leu$reduced_dims$pca
sampleinfo_x <- sampleinfo_x[match(sampleinfo_x$sample_id, rownames(pca_leu)), ]
pca_leu$anno <- sampleinfo_x$anno
pca_leu$anno <- as.factor(pca_leu$anno)


pdf("/cluster/home/yjliu_jh/share/pca_leu_anno.pdf", width = 7, height = 8)
ggscatter(pca_leu, x = "PC1", y = "PC2", color = "anno", size = 0.4)
dev.off()

sanno <- sampleinfo_r[, c("sample_id", "dataset")][sampleinfo_r$sample_id %in% c(sp_aml_samples, sp_ball_samples, sp_tall_samples), ]


# ==== total heatmap from jh ===

# check distinct 3 samples left and that cluster right



sub_heat2 <- function(ht, ind, sinfo, cname = F, rname = F){
  
  column_order <- column_order(ht)
  row_order <- row_order(ht)
  data_mat <- ht@ht_list$` `@matrix
  ordered_fil_mat <- data_mat[row_order, column_order]
  sub_mat <- ordered_fil_mat[, ind]
  
  htinfo_sub <- sinfo[sinfo$sample_id %in% colnames(sub_mat), ]
  htinfo_sub <- htinfo_sub[match(colnames(sub_mat), htinfo_sub$sample_id), ]
  
  ha_fil_sub <- HeatmapAnnotation(dataset = htinfo_sub$dataset,
                                  disease_type = htinfo_sub$disease_type,
                                  disease_type_lyj = htinfo_sub$disease_type_lyj,
                                  col = list(dataset = setNames(allcolour[1:length(na.omit(unique(htinfo_sub$dataset)))], 
                                                                na.omit(unique(htinfo_sub$dataset))),
                                             disease_type = setNames(allcolour[1:length(na.omit(unique(htinfo_sub$disease_type)))], 
                                                                na.omit(unique(htinfo_sub$disease_type))),
                                             disease_type_lyj = setNames(allcolour[1:length(na.omit(unique(htinfo_sub$disease_type_lyj)))], 
                                                                na.omit(unique(htinfo_sub$disease_type_lyj)))),
                                  simple_anno_size = unit(0.4, "cm"),
                                  annotation_name_gp = gpar(fontsize = 9))
  
  h <- ComplexHeatmap::Heatmap(sub_mat, 
                               name = " ",
                               show_row_names = rname, 
                               use_raster = F,
                               show_column_names = cname, 
                               cluster_rows = F,
                               cluster_columns = F,
                               column_title = "jh_heatmap",
                               heatmap_legend_param = list(
                                 legend_direction = "horizontal", 
                                 legend_width = unit(3, "cm")),
                               top_annotation = ha_fil_sub)
  return(h)
}








total1 <- readr::read_rds("/cluster/home/flora_jh/projects/data_preprocess/1.fusiongenes/htmap/various_disease_types.rds")
png("/cluster/home/yjliu_jh/share/test_total1.png", width = 2400, height = 1200, units = "px")
total1
dev.off()

pdf("/cluster/home/yjliu_jh/share/test_total1.pdf", width = 24, height = 12)
total1
dev.off()


sp_tt1 <- sub_heat2(total1, c(20:120, 9000:9211), sampleinfo_r)

htsp_tt1 <- draw(sp_tt1, show_heatmap_legend = FALSE)
png("/cluster/home/yjliu_jh/share/test_total_sub1.png", width = 2400, height = 1200, units = "px")
htsp_tt1
dev.off()

info_tt1 <- get_info(total1, c(20:120, 9000:9211), sampleinfo_r)
three_samples <- info_tt1$sample_id[head(which(!is.na(info_tt1$disease_type)), 3)]
sp_cluster <- info_tt1$sample_id[245:313]


sampleinfo_x <- sampleinfo_r
sampleinfo_x$anno <- "others"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% three_samples] <- "left_three"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% clusterx] <- "right_cluster"

pca_leu <- meta_leu$reduced_dims$pca
sampleinfo_x <- sampleinfo_x[match(sampleinfo_x$sample_id, rownames(pca_leu)), ]
pca_leu$anno <- sampleinfo_x$anno
pca_leu$anno <- as.factor(pca_leu$anno)
palette_anno <- c("gray", "blue", "red")
names(palette_anno) <- c("others", "right_cluster", "left_three")

pdf("/cluster/home/yjliu_jh/share/pca_leu_annos2.pdf", width = 7, height = 8)
ggscatter(pca_leu, x = "PC1", y = "PC2", color = "anno", palette = palette_anno, size = 0.3)
dev.off()





# read data 
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
# samples likely from the same biospecimen
remove2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/remove_samples2.rds")
colData(leu)$remove[colData(leu)$sample_id %in% remove2] <- "yes"
readr::write_rds(leu, rds_fn)



clusterx <- c("mecom", "TRF067197a", "TRF058850", "SJMMBALL056740_D1", "SJMMBALL056742_D1", "SJMMBALL056741_D1",
  "SJMMBALL056743_D1", "SJMMALL057910_D1", "SJMMNORM056746_G1", "SJMMNORM056747_G1", "SJMMNORM056745_G1",
  "SRR14580584", "SRR14580598", "suzhou_aml_A17", "suzhou_aml_A21", "suzhou_aml_A29", "suzhou_aml_A16",
  "suzhou_aml_A12", "HRR224912", "suzhou_aml_A96", "SRR2239693", "suzhou_aml_A77", "suzhou_aml_A50", "HRR224911",
  "HRR178294", "HRR230574", "SRR8418476", "vhaa56_1", "SRR10734555", "SRR10734939", "vhaa3", "vhaa45", "vhaa12",
  "vhaa49", "HRR224854", "HRR224834", "HRR224669", "suzhou_aml_A11", "SRR14593315", "SRR12274321", "rarg6", "UPN34",
  "HRR050398", "HRR050397", "HRR050395", "HRR050396", "SRR12274339", "SRR12274337", "SRR12274338", "SRR12274333",
  "SRR12274335", "SRR12274334", "SRR12274336", "SRR12274340", "SRR12274341", "SRR12274342", "SRR12274331", "SRR12274324",
  "SRR12274322", "SRR12274319", "SRR12274323", "SRR12274320", "SRR12274332", "SRR12274326", "SRR12274330", "SRR12274329",
  "SRR12274325", "SRR12274328", "SRR12274327")





sampleinfo_x <- sampleinfo_r
sampleinfo_x$anno <- "others"
sampleinfo_x$anno[sampleinfo_x$sample_id %in% sp_sam2] <- "special_ukall"

pca_leu <- meta_leu$reduced_dims$pca
sampleinfo_x <- sampleinfo_x[match(sampleinfo_x$sample_id, rownames(pca_leu)), ]
pca_leu$anno <- sampleinfo_x$anno
pca_leu$anno <- as.factor(pca_leu$anno)
palette_anno <- c("gray", "red")
names(palette_anno) <- c("others", "special_ukall")

pdf("/cluster/home/yjliu_jh/share/pca_leu_annos3.pdf", width = 7, height = 8)
ggscatter(pca_leu, x = "PC1", y = "PC2", color = "anno", palette = palette_anno, size = 0.3)
dev.off()




total_mat <- total1@ht_list$` `@matrix



info_path <- "/cluster/home/yjliu_jh/projects/leu_j/data/leukemia_anno_new.xlsx"
collected_info_sheets <- readxl::excel_sheets(info_path)
collected_info <- lapply(collected_info_sheets, function(X) readxl::read_excel(info_path, sheet = X))
names(collected_info) <- collected_info_sheets


coln_info <- list()
dtst <- character()
for (i in 1:length(collected_info)) {
  coln_info[[i]] <- colnames(collected_info[[i]])
  dtst <- c(dtst, rep(collected_info_sheets[i], ncol(collected_info[[i]])))
}

col_match <- data.frame(dataset = dtst, col = unlist(coln_info))
colns <- unique(col_match$col)



# ============ age =================
# collect age columns and change days to years
tissue_related_colnames <- colns[grep("tissue", colns, ignore.case = T)]
# unify colnames to "tissue"
trc <- tissue_related_colnames[c(1, 2, 3, 4, 6, 8)]

for(i in 1:length(trc)){
  tissue_datasets <- col_match$dataset[col_match$col %in% trc[i]]
  for(j in 1:length(tissue_datasets)){
    colnames(collected_info[[tissue_datasets[j]]])[colnames(collected_info[[tissue_datasets[j]]]) %in% trc[i]] <- "tissue"
  }
}

cinfo_tissue <- list()
for (i in 1:length(collected_info)) {
  cinfo_tissue[[names(collected_info)[i]]] <- collected_info[[i]][, intersect(colnames(collected_info[[i]]), c("sample_id", "tissue"))]
}

test1 <- as.data.frame(bind_rows(cinfo_tissue))
write_rds(test1, "/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/rnaseq/samples/annotated_tissue.rds")



tissuex <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/rnaseq/samples/annotated_tissue.rds")
tissuex <- tissuex[!duplicated(tissuex$sample_id), ] 
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
col_leu <- as.data.frame(colData(leu))
col_leu <- left_join(col_leu, tissuex)
colData(leu)$tissue <- col_leu$tissue
readr::write_rds(leu, rds_fn)


