# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# read data generated from other scripts
htinfo <- readr::read_rds("/cluster/home/yjliu_jh/share/sinfo_filmat_1219.rds")
fil_mat <- readr::read_rds("/cluster/home/yjliu_jh/share/filtered_leu_mat_1219.rds")
ht_fil <- readr::read_rds("/cluster/home/yjliu_jh/share/filtered_heatmap_object.rds")
sinfo3 <- readr::read_rds("/cluster/home/yjliu_jh/share/sinfo_8486.rds")
of_fil <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_fil_221217.rds")

# get t-all like samples and data, variables from clustering1217.R
column_order_htfil <- column_order(ht_fil)
sample_ordered <- colnames(fil_mat)[column_order_htfil]
ctree_ordered <- ctree[sample_ordered]
table(tail(ctree_ordered, 500))  ## confirm that cluster 4 represents t-all
tall_like_samples <- names(ctree[ctree == 4])
table(htinfo[htinfo$sample_id %in% tall_like_samples, ]$disease_type) ## double check the annotation
tall_mat <- fpkm_fil[, tall_like_samples]
tall_info2 <- sinfo3[sinfo3$sample_id %in% colnames(tall_mat), ]


# order matrix according to variance (cr: flora)
var <- matrixStats::rowVars(tall_mat)
fpkm_var <- cbind(tall_mat, var)
fpkm_ordered <- as.data.frame(fpkm_var) %>% dplyr::arrange(desc(var))
fpkm_mat_tall <- fpkm_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 

# plot PCA
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




# ==== ly PCA ===

# read data
aml_mat <- readr::read_rds("/cluster/home/yliang_jh/projects/mRNA/leukemia/mat_expr_top_filter2.rds")
aml_info <- readr::read_rds("/cluster/home/yliang_jh/projects/mRNA/leukemia/df_anno.rds")
identical(colnames(aml_mat), rownames(aml_info))

# read overall expression
overall_fpkm <- readr::read_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_short_221216.csv")
of_all <- overall_fpkm %>% dplyr::select(-c(gene_id)) %>% as.data.frame() %>% na.omit() %>%
  remove_rownames() %>% column_to_rownames(var = "gene_name") %>% as.matrix()

rm(list = "overall_fpkm")
gc()


# select datasets to remove re-naming of identical sample ids among datasets
sinfo4 <- sinfo3[sinfo3$dataset %in% aml_info$dataset, ]
of_aml_fil <- of_all[, sinfo4$sample_id]
colnames(of_aml_fil) <- sub("\\.y", "", colnames(of_aml_fil))
aml_mat_unscaled <- t(of_aml_fil[rownames(aml_mat), colnames(aml_mat)])


# plot PCA
pca_aml <- prcomp(aml_mat_unscaled, scale. = TRUE)
p_aml_temp <- autoplot(pca_aml, data = aml_info, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/aml_cluster_pca_dataset_1000.pdf", width = 6, height = 6)
p_aml_temp
dev.off()

p_aml_temp2 <- autoplot(pca_aml, data = aml_info, colour = "fixed_rna_type")

pdf("/cluster/home/yjliu_jh/share/aml_cluster_pca_rnatype_1000.pdf", width = 6, height = 6)
p_aml_temp2
dev.off()



# ==== ball PCA ===

ball_mat <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/jhuang/human/ball/left_mat_var_fil.rds")
ball_mat <- ball_mat[, -1]  ## remove the column for ordering
ball_info <- sinfo3[sinfo3$sample_id %in% colnames(ball_mat), ]
ball_info <- ball_info[match(colnames(ball_mat), ball_info$sample_id), ]

# plot PCA
pca_ball <- prcomp(t(ball_mat[1:1000, ]), scale. = TRUE)
p_ball_temp <- autoplot(pca_ball, data = ball_info, colour = "dataset")

pdf("/cluster/home/yjliu_jh/share/ball_cluster_pca_dataset_1000.pdf", width = 8, height = 6)
p_ball_temp
dev.off()

p_ball_temp2 <- autoplot(pca_ball, data = ball_info, colour = "fixed_rna_type")

pdf("/cluster/home/yjliu_jh/share/ball_cluster_pca_rnatype_1000.pdf", width = 6, height = 6)
p_ball_temp2
dev.off()




