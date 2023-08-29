# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# read counts data for all datasets (cr: flora)
datasets <- dir_ls("/cluster/home/jhuang/projects/leukemia/analysis", type = "directory") |>
  path_file()
datasets <- datasets[!datasets %in% c("bak", "diagnosis")]
mat <- list()
merge_func <- function(x,y){left_join(x,y, by = c("gene_id", "gene_name"))}
for (dataset in datasets){
  csv <- glue("/cluster/home/jhuang/projects/leukemia/analysis/{dataset}/human/rnaseq/exp/tables/{dataset}_human.csv")
  if (file.exists(csv)){
    exp <- read_csv(csv)
    if (csv == "/cluster/home/jhuang/projects/leukemia/analysis/ebiomedicine/human/rnaseq/exp/tables/ebiomedicine_human.csv"){
      colnames(exp)[-c(1,2)] <- glue("ebio_{colnames(exp)[-c(1,2)]}") #change ebiomedicine sample names
    }
    mat[[dataset]]  <- exp
  }
}

data_n <- numeric()
for (i in 1:length(mat)){
  data_n[i] <- sum(ncol(mat[[i]]))
}


# output
overall_exp <- mat[which(data_n > 0)] %>% Reduce(f = merge_func) %>% dplyr::select(-starts_with("gene_name.")) %>%
  relocate(gene_id, gene_name, everything())
#sum(colnames(overall_exp) %in% c("gene_id", "gene_name"))  ## 2
overall_exp %>% write_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_exp_221217.csv")


# read sampleinfo
overall_exp <- readr::read_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_exp_221217.csv")
sinfo <- readr::read_tsv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/sampleinfo.txt")
batch <- sinfo[, c("sample_id", "dataset")]
batch <- batch[match(colnames(overall_exp)[-(1:2)], batch$sample_id), ]
batch$dataset[is.na(batch$dataset)] <- "Unknown"
batch$sample_id[is.na(batch$sample_id)] <- setdiff(colnames(overall_exp)[-(1:2)], batch$sample_id)
batch$dataset <- as.factor(batch$dataset)

# load hugo annotation of genes
hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/hgnc_complete_set_2022-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id")]
colnames(hugo_anno)[3] <- "gene_id"


# filter blacklist-genes and select for protein-coding genes
of_fil <- filter_black_heatmap_genes(overall_exp)
of_fil <- of_fil[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
of_fil <- of_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol") %>% as.matrix()
of_fil %>% readr::write_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_fil_221217.rds")

# order matrix according to variance (cr: flora)
var <- matrixStats::rowVars(of_fil)
fpkm_var <- cbind(of_fil, var)
fpkm_ordered <- as.data.frame(fpkm_var) %>% dplyr::arrange(desc(var))
fpkm_mat <- fpkm_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 


# get clustered cll samples from cl_jh.R
ordered_mat <- readr::read_rds("/cluster/home/yjliu_jh/share/ordered_leu_mat_1217.rds")
cllsamples <- colnames(ordered_mat)[7173:7448]
solidsamples <- setdiff(colnames(mat[["GSE142520"]]), c("gene_id", "gene_name"))
fil_mat <- fpkm_mat[, colnames(fpkm_mat) %notin% cllsamples]
fil_mat <- fil_mat[, colnames(fil_mat) %notin% solidsamples]




# ==== temp, before blacklist has been updated ===
nbgenes <- c("IGHV5−10−1", "IGHV1−69D", "IGHV1−2", "IGHV3−30", "IGKV3−11", "IGKV3−15",
             "IGKV1−5", "IGKV3−20", "IGLV2−23", "IGLV2−11", "IGLV2−14", "IGLV2−8", "IGKV1D−39",
             "IGKV1−33", "IGKV1−12", "IGLV1−44", "IGLV1−51", "IGLV1−40", "IGLV3−21", "IGKV4−1",
             "IGHV3−23", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLL5", "IGHG1", "IGHG3", "IGHV3−33",
             "IGHA2", "IGHA1", "IGHG4", "IGHG2", "IGKV2−28",
             "OPALIN",  "OR10Z1",  "CCDC200", "MASP2",   "SPDYE2",  "LHFPL5",  "MYH11",
             "MT-ND6",  "NRIP2",   "MTMR7",   "LHX4",    "LINGO3",  "PRR33"
)
fil_mat <- fil_mat[rownames(fil_mat) %notin% nbgenes, ]
# ==== end ===

# order annotation data to match heatmap data
sinfo2 <- left_join(batch, sinfo[, 1:4])
sinfo2[is.na(sinfo2)] <- "Unknown"
htinfo <- data.frame(sample_id = colnames(fil_mat))
htinfo <- left_join(htinfo, sinfo2)

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

# set annotation
ha_fil <- HeatmapAnnotation(disease = htinfo$disease_type,
                         rna_type = htinfo$fixed_rna_type,
                         dataset = htinfo$dataset,
                         gender = htinfo$fixed_gender,
                         col = list(disease = setNames(allcolour[1:length(unique(htinfo$disease_type))], 
                                                       unique(htinfo$disease_type)),
                                    rna_type = setNames(allcolour[1:length(unique(htinfo$fixed_rna_type))], 
                                                        unique(htinfo$fixed_rna_type)),
                                    dataset = setNames(allcolour[1:length(unique(htinfo$dataset))], 
                                                       unique(htinfo$dataset)), 
                                    gender = setNames(c("blue", "pink", "gray", "lightgray"), 
                                                      unique(htinfo$fixed_gender))
                         ),
                         simple_anno_size = unit(0.2, "cm"),
                         annotation_name_gp = gpar(fontsize = 7))


# heatmap
h_fil <- ComplexHeatmap::Heatmap(fil_mat[1:1000, ], 
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
                               top_annotation = ha_fil)


ht_fil <- draw(h_fil, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_fil_2.png", width = 1800, height = 1500, units = "px")
ht_fil
dev.off()

pdf("/cluster/home/yjliu_jh/share/cluster_jh_fil_2.pdf", width = 18, height = 15)
ht_fil
dev.off()





# ========== sth else =========================

# check clustered 

column_order_h3 <- column_order(ht_fil)
row_order_h3 <- row_order(ht_fil)
ordered_fil_mat <- fil_mat[1:1000, ][row_order_h3, column_order_h3]
info_reorder <- htinfo[column_order_h3, ]
tail(info_reorder)


# check which dataset that specific sample belongs to 

temp1 <- character()
for(i in 1:length(mat)){
  tempp <- rep(names(mat)[i], ncol(mat[[i]]) - 2)
  temp1 <- c(temp1, tempp)
}

temp1[which(htinfo$sample_id %in% "SRR10759155")]


# ===== functions for quick subsetting and plotting =======

new_fil <- function(htinfo){
  ha_fil <- HeatmapAnnotation(disease = htinfo$disease_type,
                              rna_type = htinfo$fixed_rna_type,
                              dataset = htinfo$dataset,
                              gender = htinfo$fixed_gender,
                              col = list(disease = setNames(allcolour[1:length(unique(htinfo$disease_type))], 
                                                            unique(htinfo$disease_type)),
                                         rna_type = setNames(allcolour[1:length(unique(htinfo$fixed_rna_type))], 
                                                             unique(htinfo$fixed_rna_type)),
                                         dataset = setNames(allcolour[1:length(unique(htinfo$dataset))], 
                                                            unique(htinfo$dataset)), 
                                         gender = setNames(c("blue", "pink", "gray", "lightgray")[1:length(unique(htinfo$fixed_gender))], 
                                                           unique(htinfo$fixed_gender))
                              ),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp = gpar(fontsize = 7))
  return(ha_fil)
}




quick_heat_h <- function(dat, anno){
  h <- ComplexHeatmap::Heatmap(dat, 
                               name = " ",
                               show_row_names = T, 
                               use_raster = F,
                               show_column_names = FALSE, 
                               clustering_method_columns = 'ward.D2',
                               clustering_method_rows = 'ward.D2',
                               show_row_dend = F,
                               show_column_dend = T,
                               column_dend_height = unit(8, "cm"),
                               column_title = "jh_heatmap",
                               heatmap_legend_param = list(
                                 legend_direction = "horizontal", 
                                 legend_width = unit(3, "cm")),
                               top_annotation = anno)
  return(h)
}



quick_heat_r <- function(dat, anno){
  h <- ComplexHeatmap::Heatmap(dat, 
                                   name = " ",
                                   show_row_names = T, 
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
                                   top_annotation = anno)
  return(h)
}


quick_heat_n <- function(dat, anno){
  h <- ComplexHeatmap::Heatmap(dat, 
                               name = " ",
                               show_row_names = T, 
                               use_raster = T,
                               cluster_rows = F,
                               cluster_columns = F,
                               show_column_names = FALSE, 
                               column_title = "jh_heatmap",
                               heatmap_legend_param = list(
                                 legend_direction = "horizontal", 
                                 legend_width = unit(3, "cm")),
                               top_annotation = anno)
  return(h)
}



# ==== temp 1 - try to locate clustered genes =====

sub_mat1 <- ordered_fil_mat[350:390, 1:555]
htinfo_sub1 <- htinfo[htinfo$sample_id %in% colnames(sub_mat1), ]
ha_fil_sub1 <- new_fil(htinfo_sub1)
h_sub1 <- quick_heat_n(sub_mat1, ha_fil_sub1)

ht_sub1 <- draw(h_sub1, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_fil_sub1.png", width = 1800, height = 1500, units = "px")
ht_sub1
dev.off()


# ==== temp2 - 13 genes =====

column_order_s1 <- column_order(ht_sub1)
row_order_s1 <- row_order(ht_sub1)
ordered_s1_mat <- sub_mat1[row_order_s1, column_order_s1]
which(rownames(ordered_s1_mat) %in% "OPALIN")
in_genes <- rownames(ordered_s1_mat)[22:34]


sub_mat2 <- ordered_fil_mat[in_genes, ]
htinfo_sub2 <- htinfo[htinfo$sample_id %in% colnames(sub_mat2), ]
htinfo_sub2 <- htinfo_sub2[match(colnames(ordered_fil_mat), htinfo_sub2$sample_id), ]
ha_fil_sub2 <- new_fil(htinfo_sub2)


h_sub2 <- quick_heat_n(sub_mat2, ha_fil_sub2)

h_sub2 <- quick_heat_n(sub_mat2, ha_fil)

ha_fil_test@anno_list$dataset <- ha_fil@anno_list$dataset
h_sub2 <- quick_heat_n(sub_mat2, ha_fil_test)

ht_sub2 <- draw(h_sub2, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_fil_sub2.png", width = 2200, height = 500, units = "px")
ht_sub2
dev.off()





# oe <- overall_exp[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
# oe <- oe %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
#   remove_rownames() %>% column_to_rownames(var = "symbol") %>% as.matrix()
# 
# oe1 <- as.matrix(overall_exp[, -(1:2)])
# var1 <- matrixStats::rowVars(oe1)
# fpkm_var1 <- cbind(oe1, var1)
# fpkm_ordered1 <- as.data.frame(fpkm_var1) %>% dplyr::arrange(desc(var1))
# fpkm_mat_all <- fpkm_ordered1 %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
#   as.matrix() 



# ===== histone 1 =====

his_genes1 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/target_pattern_extracted_genes.xlsx", sheet = 1, col_names = F)
his_genes1 <- his_genes1[[1]]

oe1 <- overall_exp[overall_exp$gene_name %in% his_genes1, -1] %>% as.data.frame() %>% 
  column_to_rownames(var = "gene_name") %>% as.matrix()
oe_mat1 <- oe1[, colnames(oe1) %in% colnames(fil_mat)]


htinfo_sub3 <- htinfo[htinfo$sample_id %in% colnames(oe_mat1), ]
htinfo_sub3 <- htinfo_sub3[match(colnames(oe_mat1), htinfo_sub3$sample_id), ]
ha_fil_sub3 <- new_fil(htinfo_sub3)
h_sub3 <- quick_heat_h(oe_mat1, ha_fil_sub3)

ht_sub3 <- draw(h_sub3, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_fil_sub3.png", width = 1800, height = 1500, units = "px")
ht_sub3
dev.off()



# ==== histone2 ======

his_genes2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/target_pattern_extracted_genes.xlsx", sheet = 2, col_names = F)
his_genes2 <- his_genes2[[1]]

oe2 <- overall_exp[overall_exp$gene_name %in% his_genes2, -1] %>% as.data.frame() %>% 
  column_to_rownames(var = "gene_name") %>% as.matrix()
sub_mat4 <- oe2[, colnames(oe2) %in% colnames(fil_mat)]


htinfo_sub4 <- htinfo[htinfo$sample_id %in% colnames(sub_mat4), ]
htinfo_sub4 <- htinfo_sub4[match(colnames(sub_mat4), htinfo_sub4$sample_id), ]
ha_fil_sub4 <- new_fil(htinfo_sub4)
h_sub4 <- quick_heat_h(sub_mat4, ha_fil_sub4)

ht_sub4 <- draw(h_sub4, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_fil_sub4.png", width = 1800, height = 1500, units = "px")
ht_sub4
dev.off()


# ==== sex genes (add xist?) ======


black_path <- fs::path_package("jhdata", "extdata/genesets/heatmap_blackgeneset.xlsx")
sex_genes <- readxl::read_excel(black_path, sheet = 3)  # , col_names = F
sex_genes <- c(sex_genes[[1]][c(1:4, 6:8)], "XIST")


overall_exp[ overall_exp$gene_id %in% "ENSG00000229807", 2] <- "XIST"
oe3 <- overall_exp[overall_exp$gene_name %in% sex_genes, -1] %>% as.data.frame() %>% 
  column_to_rownames(var = "gene_name") %>% as.matrix()
sub_mat5 <- oe3[, colnames(oe3) %in% colnames(fil_mat)]


htinfo_sub5 <- htinfo[htinfo$sample_id %in% colnames(sub_mat5), ]
htinfo_sub5 <- htinfo_sub5[match(colnames(sub_mat5), htinfo_sub5$sample_id), ]
ha_fil_sub5 <- new_fil(htinfo_sub5)
h_sub5 <- quick_heat_h(sub_mat5, ha_fil_sub5)

ht_sub5 <- draw(h_sub5, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_fil_sub5.png", width = 1800, height = 1500, units = "px")
ht_sub5
dev.off()



