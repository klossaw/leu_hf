pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



# kinase-activating alterations summarized in a previous study
phlike_alter <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_nejm1_25207766.xlsx", 3)



leu_si_path <- "/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/leukemia_pure_rnaseq.xlsx"
tab_names <- readxl::excel_sheets(path = leu_si_path)
leu_si_all <- lapply(tab_names, function(x) readxl::read_excel(path = leu_si_path, sheet = x))

mcols <- c("sample_id", "disease_type", "dataset")

for (i in 1:length(leu_si_all)){
  if ("dataset" %notin% colnames(leu_si_all[[i]])) {
    leu_si_all[[i]]$dataset <- tab_names[i]
  }
}

leu_s_all <- list()
for (i in 1:length(leu_si_all)){
  leu_s_all[[i]] <- leu_si_all[[i]][, mcols]
}

leu_s <- data.table::rbindlist(leu_s_all)


temp8 <- get_classes(rlb["ATC:skmeans"], k = 8)
temp8$sample_id <- rownames(temp8)
temp8 <- left_join(temp8, leu_s)
table(temp8[, c("class", "disease_type")])








overall_exp <- readr::read_csv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/overall_exp_counts.csv")

# annotate ensembl gene id to symbols of coding genes
hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/hgnc_complete_set_2022-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id", "entrez_id")]
colnames(hugo_anno)[3] <- "gene_id"
oe_fil <- overall_exp[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
oe_fil <- oe_fil %>% dplyr::select(-c(gene_id, locus_group, entrez_id)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
oe_fil <- as.matrix(round(oe_fil))


# annotation and exp data didn't match, get intersect
leu_c <- leu_s[leu_s$sample_id %in% intersect(leu_s$sample_id, colnames(oe_fil)), ]
oe_c <- oe_fil[, intersect(leu_s$sample_id, colnames(oe_fil))]

aml_samples <- leu_c$sample_id[leu_c$disease_type %in% "aml"]
ball_samples <- leu_c$sample_id[leu_c$disease_type %in% "ball"]
mpal_samples <- leu_c$sample_id[leu_c$disease_type %in% "mpal"]
tall_samples <- leu_c$sample_id[leu_c$disease_type %in% "tall"]

aml_exp <- oe_c[, aml_samples]
ball_exp <- oe_c[, ball_samples]
mpal_exp <- oe_c[, mpal_samples]
tall_exp <- oe_c[, tall_samples]



# ======= clustering using overall data ======= 



adjusted_c <- ComBat_seq(counts = oe_c, batch = leu_c$dataset)


# 1 test the correlation   2 get statistics like mean and sd   3 substract and check summary for values
#> cor(adjusted_c[1, ], oe_c[1, ])
#[1] 0.7408264
#> cor(adjusted_c[, 1], oe_c[, 1])
#[1] 0.9289143
#> cor(adjusted_c[, 4], oe_c[, 4])
#[1] 0.4969646
#> cor(adjusted_c[4, ], oe_c[4, ])
#[1] 0.903339

# summary also big differences, but median always zero
# next try if normalized expression changes much





# if the difference is low, may build centroids for subsequent new data? without 














# construct dds object and get DESeq-normalized matrix
ddsc <- DESeqDataSetFromMatrix(countData = oe_c,
                               colData = leu_c,
                               design = ~ dataset) 

ddsa <- DESeqDataSetFromMatrix(countData = adjusted_c,
                               colData = leu_c,
                               design = ~ dataset) 

vsd1 <- vst(ddsc, blind = FALSE)
vsd2 <- vst(ddsa, blind = FALSE)



norm_exp_all <- assay(vsd1)
norm_exp_adj <- assay(vsd2)



oe_mat <- adjust_matrix(norm_exp_all)
osd <- transform(oe_mat, SD=apply(oe_mat, 1, sd, na.rm = TRUE))
oe_mat2 <- oe_mat[order(osd$SD, decreasing = T), ][1:3000, ]

oe_mat3 <- adjust_matrix(norm_exp_adj)
osd <- transform(oe_mat3, SD=apply(oe_mat3, 1, sd, na.rm = TRUE))
oe_mat4 <- oe_mat3[order(osd$SD, decreasing = T), ][1:3000, ]




rlb1_aml = run_all_consensus_partition_methods(oe_mat2[, aml_samples], max_k = 12,
                                          top_value_method = c("SD", "ATC"),
                                          partition_method = c("kmeans", "skmeans", "hclust"),
                                          cores = 40)
readr::write_rds(rlb1_aml, "/cluster/home/yjliu_jh/projects/leu_j/data/aml_exp_clu_3000.rds")

rlb1_ball = run_all_consensus_partition_methods(oe_mat2[, ball_samples], max_k = 12,
                                               top_value_method = c("SD", "ATC"),
                                               partition_method = c("kmeans", "skmeans", "hclust"),
                                               cores = 40)
readr::write_rds(rlb1_ball, "/cluster/home/yjliu_jh/projects/leu_j/data/ball_exp_clu_3000.rds")

rlb1_mpal = run_all_consensus_partition_methods(oe_mat2[, mpal_samples], max_k = 12,
                                               top_value_method = c("SD", "ATC"),
                                               partition_method = c("kmeans", "skmeans", "hclust"),
                                               cores = 40)
readr::write_rds(rlb1_mpal, "/cluster/home/yjliu_jh/projects/leu_j/data/mpal_exp_clu_3000.rds")

rlb1_tall = run_all_consensus_partition_methods(oe_mat2[, tall_samples], max_k = 12,
                                               top_value_method = c("SD", "ATC"),
                                               partition_method = c("kmeans", "skmeans", "hclust"),
                                               cores = 40)
readr::write_rds(rlb1_tall, "/cluster/home/yjliu_jh/projects/leu_j/data/tall_exp_clu_3000.rds")


rlb2_aml = run_all_consensus_partition_methods(oe_mat4[, aml_samples], max_k = 12,
                                               top_value_method = c("SD", "ATC"),
                                               partition_method = c("kmeans", "skmeans", "hclust"),
                                               cores = 40)
readr::write_rds(rlb2_aml, "/cluster/home/yjliu_jh/projects/leu_j/data/aml_adj_clu_3000.rds")

rlb2_ball = run_all_consensus_partition_methods(oe_mat4[, ball_samples], max_k = 12,
                                                top_value_method = c("SD", "ATC"),
                                                partition_method = c("kmeans", "skmeans", "hclust"),
                                                cores = 40)
readr::write_rds(rlb2_ball, "/cluster/home/yjliu_jh/projects/leu_j/data/ball_adj_clu_3000.rds")

rlb2_mpal = run_all_consensus_partition_methods(oe_mat4[, mpal_samples], max_k = 12,
                                                top_value_method = c("SD", "ATC"),
                                                partition_method = c("kmeans", "skmeans", "hclust"),
                                                cores = 40)
readr::write_rds(rlb2_mpal, "/cluster/home/yjliu_jh/projects/leu_j/data/mpal_adj_clu_3000.rds")

rlb2_tall = run_all_consensus_partition_methods(oe_mat4[, tall_samples], max_k = 12,
                                                top_value_method = c("SD", "ATC"),
                                                partition_method = c("kmeans", "skmeans", "hclust"),
                                                cores = 40)
readr::write_rds(rlb2_tall, "/cluster/home/yjliu_jh/projects/leu_j/data/tall_adj_clu_3000.rds")




# only 1 day needed for the run











#ddso <- estimateSizeFactors(ddso)
#normalized_counts1 <- counts(ddso, normalized=TRUE)
#readr::write_rds(ddso, "/cluster/home/yjliu_jh/projects/leu_j/data/dds_overall.rds")


# use cola on DESeq-altered exp matrix
overall_exp <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/exp_norm_all.rds")
oe_mat <- adjust_matrix(overall_exp)
osd <- transform(oe_mat, SD=apply(oe_mat, 1, sd, na.rm = TRUE))
oe_mat2 <- oe_mat[order(osd$SD, decreasing = T), ][1:3000, ]

rlb = run_all_consensus_partition_methods(oe_mat2, max_k = 12,
                                          top_value_method = c("SD", "ATC"),
                                          partition_method = c("kmeans", "skmeans", "hclust"),
                                          cores = 40)
readr::write_rds(rlb, "/cluster/home/yjliu_jh/projects/leu_j/data/new_exp_clu_3000.rds")


# use cola on batch-corrected and DESeq-altered exp matrix
overall_exp <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/exp_corr_norm_all.rds")
oe_mat <- adjust_matrix(overall_exp)
osd <- transform(oe_mat, SD=apply(oe_mat, 1, sd, na.rm = TRUE))
oe_mat2 <- oe_mat[order(osd$SD, decreasing = T), ][1:3000, ]

rlb2 = run_all_consensus_partition_methods(oe_mat2, max_k = 12,
                                           top_value_method = c("SD", "ATC"),
                                           partition_method = c("kmeans", "skmeans", "hclust"),
                                           cores = 40)
readr::write_rds(rlb2, "/cluster/home/yjliu_jh/projects/leu_j/data/new_exp_corr_clu_3000.rds")







# test if count-based clustering is comparable with exp-based
#rls1 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/rls.rds")
rls2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/rls_exp.rds")
cla1 <- get_classes(rls2, k = 10)
cla2 <- get_classes(rls2, k = 9)
cla3 <- get_classes(rls2, k = 8) 
cla1$sample_id <- rownames(cla1)
cla2$sample_id <- rownames(cla2)
cla3$sample_id <- rownames(cla3)
# clusters using exp obviously better than count in concordance, according to tests (not shown here)
# note that it's single source data, batch correction is not needed


# now test if exp from salmon and deseq-version exp can have affects on clustering













cbg <- left_join(groupings, cla1)
igraph::compare(cbg$class, as.numeric(as.factor(cbg$subgroups)), method = "nmi")




cbg$class <- paste0("X", cbg$class)
cbg$comb <- paste0(cbg$subgroups, "-", cbg$class)
cbg %>% group_by(subgroups) %>% mutate(ac = n()) %>% group_by(subgroups, class) %>% 
  summarise(percentage = n() / ac) %>% unique() %>% as.data.frame()
# clusters G1 G2 G3 G4 G8 had high concordance, and these clusters are EXACTLY 
# the clusters who have distinct genomic patterns!
cbg %>% group_by(class) %>% mutate(ac = n()) %>% group_by(class, subgroups) %>% 
  summarise(percentage = n() / ac) %>% unique() %>% as.data.frame()






## use heatmaps and a sankey plot to draw the figure 




















# 
tog_g5 <- cbg[cbg$subgroups %in% "G5", ]
table(tog_g5[tog_g5$characteristic %in% "KMT2A-r", "class"])
table(tog_g5[tog_g5$characteristic %in% "NPM1", "class"])

# npm1 X2    kmt2a-r X6  (G4 NPM1 - all X3)    G3-X5 identical, PML-RARA specific




table(tog[tog$characteristic %in% "NPM1", "class"])



# batch-corrected 

adjusted_mat <- ComBat_seq(counts = oe_fil, batch = batch$dataset)
readr::write_rds(adjusted_mat, "/cluster/home/yjliu_jh/projects/leu_j/data/count_corrected_all.rds")
ddso2 <- DESeqDataSetFromMatrix(countData = adjusted_mat,
                                colData = batch,
                                design = ~ class) 
vsd2 <- vst(ddso2, blind = FALSE)
normalizedExp2 <- assay(vsd2)
readr::write_rds(normalizedExp2, "/cluster/home/yjliu_jh/projects/leu_j/data/exp_corr_norm_all.rds")



# compare batch-corrected count clustering results   count clustering results  






# compare normalized exp with TPM or DESeq2 with exp file 


# compare batch-corrected and normalized exp with exp file


# use a small set (maybe leu) to check differences


#




# compare clustering results 
# check concordance of leukemia set? fusion as gold standard?
















