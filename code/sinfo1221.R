# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}




# keep newest update of sampleinfo
leu <- readr::read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")
sampleinfo_r <- colData(leu)
sampleinfo2 <- readxl::read_excel("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/sampleinfo2.xlsx")
sampleinfo3 <- sampleinfo2[match(sampleinfo_r$sample_id, sampleinfo2$sample_id), ]

new_cols_newname <- c("age", "old_gender", "old_rna_type", "subgroup", "sample_index",
                      "outcome", "fusion_genes", "mutations", "disease_type")
new_cols <- c("age", "gender", "rna_library", "subgroup", "sample_index",
              "outcome", "gene_fusion", "mutations", "disease_type")
sampleinfo_r[, new_cols_newname] <- sampleinfo3[, new_cols]



temp <- sampleinfo_r[sampleinfo_r$gender != sampleinfo_r$old_gender, ]
temp <- temp[temp$old_gender != "NA", ]

temp2 <- sampleinfo_r[sampleinfo_r$old_rna_type != sampleinfo_r$rna_type, ]
temp2 <- temp2[temp2$old_rna_type != "NA", ]


# heatmap using that particular genesets and two annotations

overall_exp <- as.matrix(assay(leu))

allcolour = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", 
              "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", 
              "#FFFF00", "#808000", "#FF00FF", "#FA8072", "#7B68EE", 
              "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", 
              "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", 
              "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", 
              "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", 
              "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", 
              "#DEB887")


# ===== histone 1 =====

his_genes1 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/target_pattern_extracted_genes.xlsx", sheet = 1, col_names = F)
his_genes1 <- his_genes1[[1]]

sinfo_r <- colData(leu)
oe_mat2 <- overall_exp[, colnames(overall_exp) %in% temp2$sample_id]
oe_mat2 <- oe_mat2[rownames(oe_mat2) %in% his_genes1, ]
sinfo_r2 <- sinfo_r[sinfo_r$sample_id %in% temp2$sample_id, ]
identical(sinfo_r2$sample_id, colnames(oe_mat2))


ha_fil2 <- HeatmapAnnotation(old_rna_type = sinfo_r2$old_rna_type,
                            rna_type = sinfo_r2$rna_type,
                            simple_anno_size = unit(0.2, "cm"),
                            annotation_name_gp = gpar(fontsize = 7))



h_sub3 <- quick_heat_h(oe_mat2, ha_fil2)

ht_sub3 <- draw(h_sub3, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_rna_type_diff.png", width = 1800, height = 1500, units = "px")
ht_sub3
dev.off()



# ==== sex genes (add xist?) ======


black_path <- fs::path_package("jhdata", "extdata/genesets/heatmap_blackgeneset.xlsx")
sex_genes <- readxl::read_excel(black_path, sheet = 4)  # , col_names = F
sex_genes <- c(sex_genes[[1]][c(1:4, 6:8)], "XIST")

oe_mat1 <- overall_exp[, colnames(overall_exp) %in% temp$sample_id]
oe_mat1 <- oe_mat1[rownames(oe_mat1) %in% sex_genes, ]
sinfo_r1 <- sinfo_r[sinfo_r$sample_id %in% temp$sample_id, ]
identical(sinfo_r1$sample_id, colnames(oe_mat1))

ha_fil1 <- HeatmapAnnotation(old_gender = sinfo_r1$old_gender,
                            gender = sinfo_r1$gender,
                            col = list(old_gender = setNames(c("pink", "blue", "gray", "lightgray")[1:length(unique(sinfo_r1$old_gender))], 
                                                         unique(sinfo_r1$old_gender)),
                                       gender = setNames(c("blue", "pink", "gray", "lightgray")[1:length(unique(sinfo_r1$gender))], 
                                                         unique(sinfo_r1$gender))
                            ),
                            simple_anno_size = unit(0.2, "cm"),
                            annotation_name_gp = gpar(fontsize = 7))

h_sub3 <- quick_heat_h(oe_mat1, ha_fil1)

ht_sub3 <- draw(h_sub3, heatmap_legend_side="right")

png("/cluster/home/yjliu_jh/share/cluster_jh_gender_diff.png", width = 500, height = 400, units = "px")
ht_sub3
dev.off()






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




















sum(sampleinfo_r$gender != sampleinfo_r$old_gender)

meta_leu <- metadata(leu)

sum(is.na(sampleinfo_r$disease_type))


# ============ reduceddim plots ==============
pca_leu <- meta_leu$pca
sampleinfo_r2 <- sampleinfo_r[match(sampleinfo_r$sample_id, rownames(pca_leu)), ]

pca_leu$disease_type <- sampleinfo_r2$disease_type 


tb_disease <- table(pca_leu$disease_type) 
r_diseases <- names(tb_disease)[tb_disease >20]
pca_leu$disease_type[pca_leu$disease_type %notin% r_diseases] <- "Others"


pca_leu$disease_type <- as.factor(pca_leu$disease_type)
pca_leu <- pca_leu[pca_leu$PC1 > -100, ]





pdf("/cluster/home/yjliu_jh/share/pca_leu_disease_type.pdf", width = 7, height = 8)
ggscatter(pca_leu, x = "PC1", y = "PC2", color = "disease_type", size = 0.4)
dev.off()


tsne_leu <- as.data.frame(meta_leu$tsne)
umap_leu <- as.data.frame(meta_leu$umap)
rownames(tsne_leu) <- rownames(umap_leu)
#identical(rownames(tsne_leu), sampleinfo_r2$sample_id)

tsne_leu$disease_type <- sampleinfo_r2$disease_type 
umap_leu$disease_type <- sampleinfo_r2$disease_type 
tsne_leu$disease_type[tsne_leu$disease_type %notin% r_diseases] <- "Others"
umap_leu$disease_type[umap_leu$disease_type %notin% r_diseases] <- "Others"


pdf("/cluster/home/yjliu_jh/share/tsne_leu_disease_type.pdf", width = 7, height = 8)
ggscatter(tsne_leu, x = "tSNE_1", y = "tSNE_2", color = "disease_type", size = 0.4)
dev.off()


pdf("/cluster/home/yjliu_jh/share/umap_leu_disease_type.pdf", width = 7, height = 8)
ggscatter(umap_leu, x = "UMAP_1", y = "UMAP_2", color = "disease_type", size = 0.4)
dev.off()





# =========== end ===========================


# tsne x < 30 and x > 12  y < -30    ->> cll
tsne_leu[tsne_leu$tSNE1 < 30 & tsne_leu$tSNE1 > 12 & tsne_leu$tSNE2 < -30, ]$disease_type <- "cll"

tsne_leu[tsne_leu$tSNE1 < -5 & tsne_leu$tSNE1 > -20 & tsne_leu$tSNE2 > 50, ]

all_diseases <- unique(sampleinfo_r2$disease_type)

all_diseases <- c(r_diseases, setdiff(all_diseases, r_diseases))


c28 <- c(
  "gray", "gray", "gray", "gray", "gray", "gray", "gray", 
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1"
)
#pie(rep(1, 25), col = c25)
names(c28) <- all_diseases


pdf("/cluster/home/yjliu_jh/share/tsne_leu_disease_type_new.pdf", width = 7, height = 8)
ggscatter(tsne_leu, x = "tSNE1", y = "tSNE2", color = "disease_type", palette = c28, size = 0.4)
dev.off()
















sinfo3 <- readr::read_rds("/cluster/home/yjliu_jh/share/sinfo_8486.rds")
info_path <- "/cluster/home/yjliu_jh/projects/leu_j/data/leukemia_pure_rnaseq.xlsx"
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

age_related_colnames <- colns[grep("age", colns)]
age_day <- col_match$dataset[col_match$col %in% "age(days)"]
collected_info <- 
  
  # needed column: age, subgroup, index? fusion, mutation, outcome, disease_type, gender, rna_library
  
  
  mut_related_colnames <- colns[grep("mutat", colns, ignore.case=TRUE), ]
fus_related_colnames <- colns[grep("fusion", colns)]

arc <- age_related_colnames[c(1, 4, 6, 10, 11)]

head(as.data.frame(collected_info[[col_match$dataset[col_match$col %in% "RT-PCR for \r\nSPI1 fusion"]]]))


tempdt <- setdiff(unique(sinfo3$dataset), collected_info_sheets)


sinfo3[sinfo3$dataset %in% "EGAS00001002217", ]$sample_id     ## mixed diseases, such as LGG samples

sinfo3[sinfo3$dataset %in% "HRA000122", ]$sample_id

sinfo3[sinfo3$dataset %in% tempdt[7], ]$sample_id



#HRA000122       HRA000489       HRA000789       phs000218       phs000873       pnas_tall 
#124              11             292             375               2             130 





# 1 make sure dataset and sample_id exists in sinfo

# HRA000122 match sample      
# HRA000789 already have anno, ignore
# HRA000489 match anno from S5
# phs000218 output from target_aml and 218 extra and add to excel
# phs000873 
# pnas_tall 



# ====== deleted =====
# SRR1909157 should be deleted, as it's bone marrow




# needed column: age, subgroup, index? fusion, mutation, outcome, age_group, karyotype
# already had:  disease_type, gender, rna_library







newly_added_datasets <- c("sj_tall", "yj_tall", "scmcaml", "sih_m2b", "scmc_m2b", "ng_2016", "m2b_pnas")





# phs000218
targetaml_all1 <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/TARGET_AML_mRNA-seq_20170609.sdrf.txt")
phs218extra_r2s <- readr::read_csv("/cluster/home/yjliu_jh/projects/leu_j/data/218extra.csv")

targetaml_r2s <- targetaml_all1[, c("Source Name", "Comment[SRA_RUN]...46")]
phs218extra_r2s <- phs218extra_r2s[, c("submitted_subject_id", "Run")]

colnames(targetaml_r2s) <- c("case_submitter_id", "sample_id")
colnames(phs218extra_r2s) <- c("case_submitter_id", "sample_id")
tog218_r2s <- rbind(targetaml_r2s, phs218extra_r2s)

# match run to sample_id to annotation

tar218_2 <- readr::read_tsv("/cluster/home/yjliu_jh/projects/leu_j/data/218_target_all_P2.tsv")
tar218_1 <- readr::read_tsv("/cluster/home/yjliu_jh/projects/leu_j/data/218_target_aml.tsv")
tar218 <- rbind(tar218_1, tar218_2)
tar218x <- left_join(tar218, tog218_r2s)

tar218_sub <- tar218x[, c("sample_id", "ethnicity", "gender", "race", "vital_status", "age_at_diagnosis", 
                          "days_to_last_follow_up", "primary_diagnosis", "site_of_resection_or_biopsy",
                          "year_of_diagnosis")]
tar218_sub <- tar218x[tar218x$sample_id %in% sinfo3[sinfo3$dataset %in% "phs000218", ]$sample_id, ]

readr::write_csv(tar218_sub, "/cluster/home/yjliu_jh/projects/leu_j/data/phs000218_sinfo.csv")



# HRA000122



























# ======== pca plot ============


# data from other scripts
unscaled_tall_mat <- t(of_fil[, colnames(of_fil) %in% colnames(tall_mat)])
pca_tall <- prcomp(unscaled_tall_mat, scale. = TRUE)
p_tall_temp <- autoplot(pca_tall, data = tall_info2, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/tall_cluster_pca.pdf", width = 6, height = 6)
p_tall_temp
dev.off()

p_tall_temp2 <- autoplot(pca_tall, data = tall_info2, colour = "fixed_rna_type")

pdf("/cluster/home/yjliu_jh/share/tall_cluster_pca2.pdf", width = 6, height = 6)
p_tall_temp2
dev.off()


unscaled_tall_mat2 <- t(of_fil[, colnames(of_fil) %in% colnames(tall_mat)][rownames(fpkm_mat_tall)[1:1000], ])

pca_tall <- prcomp(unscaled_tall_mat2, scale. = TRUE)

p_tall_temp <- autoplot(pca_tall, data = tall_info2, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/tall_cluster_pca_dataset_1000.pdf", width = 6, height = 6)
p_tall_temp
dev.off()

p_tall_temp2 <- autoplot(pca_tall, data = tall_info2, colour = "fixed_rna_type")

pdf("/cluster/home/yjliu_jh/share/tall_cluster_pca_rnatype_1000.pdf", width = 6, height = 6)
p_tall_temp2
dev.off()









rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
fusion_tab <- readxl::read_excel("/cluster/home/flora_jh/projects/data_preprocess/1.fusiongenes/docs/fusion.xlsx")
rds <- readRDS(rds_fn)
fusion_tab <- fusion_tab %>% dplyr::select(sample_id, fusion, fusion_type)
dat <- colData(rds) %>% as.data.frame() %>% left_join(fusion_tab)
colData(rds)$fusion <- dat$fusion
colData(rds)$fusion_type <- dat$fusion_type

readr::write_rds(rds, rds_fn)

