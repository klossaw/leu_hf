# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# load data (generated from other scripts and noted below)
htinfo <- readr::read_rds("/cluster/home/yjliu_jh/share/sinfo_filmat_1219.rds")
fil_mat <- readr::read_rds("/cluster/home/yjliu_jh/share/filtered_leu_mat_1219.rds")
ht_fil <- readr::read_rds("/cluster/home/yjliu_jh/share/filtered_heatmap_object.rds")
fc <- readr::read_rds("/cluster/home/yjliu_jh/share/fusions_starfusion.rds")
ar <- readr::read_rds("/cluster/home/yjliu_jh/share/fusions_arriba.rds")
sf <- readr::read_rds("/cluster/home/yjliu_jh/share/fusions_fusioncatcher.rds")
fpkm_fil <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_fil_221217.rds")


# get t-all like samples and data, variables from clustering1217.R
column_order_htfil <- column_order(ht_fil)
sample_ordered <- colnames(fil_mat)[column_order_htfil]
ctree_ordered <- ctree[sample_ordered]
table(tail(ctree_ordered, 500))  ## confirm that cluster 4 represents t-all
tall_like_samples <- names(ctree[ctree == 4])
table(htinfo[htinfo$sample_id %in% tall_like_samples, ]$disease_type) ## double check the annotation
tall_mat <- fpkm_fil[, tall_like_samples]

# order matrix according to variance (cr: flora)
var <- matrixStats::rowVars(tall_mat)
fpkm_var <- cbind(tall_mat, var)
fpkm_ordered <- as.data.frame(fpkm_var) %>% dplyr::arrange(desc(var))
fpkm_mat_tall <- fpkm_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 


# collect fusion results (from fusion1219.R)
fc <- Filter(NROW, fc)
fc_all <- bind_rows(fc, .id = "sample_id") %>% as.data.frame()
fc_all$tool <- "fusioncatcher"
ar <- Filter(NROW, ar)
ar_all <- data.table::rbindlist(ar, idcol = "sample_id") %>% as.data.frame()
ar_all$tool <- "arriba"
sf <- Filter(NROW, sf)
sf_all <- data.table::rbindlist(sf, idcol = "sample_id") %>% as.data.frame()
sf_all$tool <- "starfusion"

# reshape the data and combine
fc_slim <- fc_all[, c("Gene_1_symbol(5end_fusion_partner)", "Gene_2_symbol(3end_fusion_partner)",
                      "Predicted_effect", "tool", "sample_id")]
colnames(fc_slim) <- c("gene1", "gene2", "anno", "tool", "sample_id")
ar_slim <- ar_all[, c("gene1", "gene2", "reading_frame", "tool", "sample_id")]
colnames(ar_slim)[3] <- "anno"
fusion_all <- rbind(fc_slim, ar_slim)

# star-fusion did not output ORF annotations
# need sth like --FusionInspector inspect --annotate --examine_coding_effect
sf_slim <- sf_all[, c("FusionName", "LargeAnchorSupport", "tool", "sample_id")] %>% 
  tidyr::separate(col = "FusionName", sep = "--", into = c("gene1", "gene2"))
colnames(sf_slim)[3] <- "anno"
fusion_all <- rbind(fusion_all, sf_slim)


# identified by at least 2 tools and are not annotated as out-of-frame, truncated or has no LDAS
fusion_fil1 <- fusion_all %>% dplyr::select(-anno) %>% unique() %>% 
  group_by(gene1, gene2, sample_id) %>% summarise(tool = toString(tool)) %>% 
  dplyr::filter(nchar(tool) > 13) %>% dplyr::select(-tool) %>% as.data.frame()
fusion_fil2 <- unique(fusion_all[fusion_all$anno %notin% c("out-of-frame", "NO_LDAS", "stop-codon") &
                                  !grepl("truncated", fusion_all$anno) ,][, c("gene1", "gene2", "sample_id")])
fusion_fil <- semi_join(fusion_fil1, fusion_fil2)
readr::write_rds(fusion_fil, "/cluster/home/yjliu_jh/share/leu_fusion_filtered.rds")


# manually select and annotate fusions (variable from clustering1217.R)
tall_info <- htinfo[htinfo$sample_id %in% colnames(tall_mat), ]
tall_info$dataset <- as.character(tall_info$dataset)

# construct a function to add fusion columns to data
process_fusions <- function(dat, fusion_data, fusion_list){
  for (i in 1:length(fusion_list)){
    fus <- strsplit(fusion_list[i], "_")[[1]] ## stringr::str_split_1(fusion_list[i], "_")
    fus <- rep(fus, 2)
    if (length(fus) > 2){
      test00 <- fusion_data[fusion_data$gene1 %in% fus[1] & fusion_data$gene2 %in% fus[2], ]$sample_id
    } else {
      test00 <- fusion_data[fusion_data$gene1 %in% fus[1] | fusion_data$gene2 %in% fus[2], ]$sample_id
    }
    dat[fusion_list[i]] <- ifelse(dat$sample_id %in% test00, "exist", "non-exist")
  }
  return(dat)
}


# generate info
fl <- c("SET_NUP214", "STIL_TAL1", "NUP98", "SPI1", "KMT2A", "MLLT10", "HOXA10",
        "TLX3", "TLX1", "NKX2-1", "TAL1", "LMO1", "LMO2", "ABL1", "LYL1", "ZFP36L2")

#tall_info_fus <- process_fusions(tall_info, fusion_fil, fl)
tall_info_fus <- process_fusions(tall_info, fusion_fil2, fl)


# generate annotation data
fus_col <- setNames(c("black", "white"), c("exist", "non-exist"))
test_fil <- HeatmapAnnotation(disease = tall_info_fus$disease_type,
                              rna_type = tall_info_fus$fixed_rna_type,
                              dataset = tall_info_fus$dataset,
                              gender = tall_info_fus$fixed_gender,
                              SET_NUP214 = tall_info_fus$SET_NUP214,
                              STIL_TAL1 = tall_info_fus$STIL_TAL1,
                              NUP98_fusions = tall_info_fus$NUP98,
                              SPI1_fusions = tall_info_fus$SPI1,
                              KMT2A_fusions = tall_info_fus$KMT2A,
                              MLLT10_fusions = tall_info_fus$MLLT10,
                              HOXA10_fusions = tall_info_fus$HOXA10,
                              ABL1_fusions = tall_info_fus$ABL1,
                              ZFP36L2_fusions = tall_info_fus$ZFP36L2,
                              TLX1_fusions = tall_info_fus$TLX1,
                              TLX2_fusions = tall_info_fus$TLX2,
                              TLX3_fusions = tall_info_fus$TLX3,
                              LMO1_fusions = tall_info_fus$LMO1,
                              LMO2_fusions = tall_info_fus$LMO2,
                              
                              col = list(disease = setNames(allcolour[1:length(unique(tall_info_fus$disease_type))], 
                                                            unique(tall_info_fus$disease_type)),
                                         rna_type = setNames(allcolour[1:length(unique(tall_info_fus$fixed_rna_type))], 
                                                             unique(tall_info_fus$fixed_rna_type)),
                                         dataset = setNames(allcolour[1:length(unique(tall_info_fus$dataset))], 
                                                            unique(tall_info_fus$dataset)), 
                                         gender = setNames(c("blue", "pink", "gray", "lightgray"), 
                                                           unique(tall_info_fus$fixed_gender)),
                                         SET_NUP214 = fus_col,
                                         STIL_TAL1 = fus_col,
                                         NUP98_fusions = fus_col,
                                         SPI1_fusions = fus_col,
                                         KMT2A_fusions = fus_col,
                                         MLLT10_fusions = fus_col,
                                         HOXA10_fusions = fus_col,
                                         ABL1_fusions = fus_col,
                                         ZFP36L2_fusions = fus_col,
                                         TLX1_fusions = fus_col,
                                         TLX2_fusions = fus_col,
                                         TLX3_fusions = fus_col,
                                         LMO1_fusions = fus_col,
                                         LMO2_fusions = fus_col
                                         
                              ),
                            simple_anno_size = unit(0.2, "cm"),
                            annotation_name_gp = gpar(fontsize = 7))


# heatmap
h_tall <- ComplexHeatmap::Heatmap(fpkm_mat_tall[1:1000, ], 
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
                                 top_annotation = test_fil)

# output heatmap
ht_tall <- draw(h_tall, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_tall_fus_nofilter.png", width = 1500, height = 2000, units = "px")
ht_tall
dev.off()

pdf("/cluster/home/yjliu_jh/share/cluster_jh_tall_fus_nofilter.pdf", width = 12, height = 18)
ht_tall
dev.off()



# extra

column_order_tall <- column_order(ht_tall)
row_order_tall <- row_order(ht_tall)
ordered_tall_mat <- fpkm_mat_tall[1:1000, ][row_order_tall, column_order_tall]
tall_info_reorder <- tall_info[column_order_tall, ]


tall_sub_mat1 <- ordered_tall_mat[909:1000, 1:555]
tall_info_sub1 <- tall_info_reorder[tall_info_reorder$sample_id %in% colnames(tall_sub_mat1), ]
ha_tall_sub1 <- new_fil(tall_info_sub1)
tall_sub1 <- quick_heat_n(tall_sub_mat1, ha_tall_sub1)

htall_sub1 <- draw(tall_sub1, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_tall_sub1.png", width = 1800, height = 1800, units = "px")
htall_sub1
dev.off()



sinfo3 <- sinfo2
sinfo3$dataset <- temp1  ## from clustering1217.R


samples_to_remove <- readr::read_csv("/cluster/home/ylxie_jh/projects/leukemia/analysis/jhuang/human/BALL/sams_remove_ylxie.csv")

