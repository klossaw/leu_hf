pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# read integrated exp data and dataset information
overall_exp <- read_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_exp_count.csv")
batch <- read_tsv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/sampleinfo.txt")
batch <- batch[, c("sample_id", "dataset")]
batch <- batch[match(colnames(overall_exp)[-(1:2)], batch$sample_id), ]
batch$dataset <- as.factor(batch$dataset)

# annotate ensembl gene id to symbols of coding genes
hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/hgnc_complete_set_2022-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id")]
colnames(hugo_anno)[3] <- "gene_id"
oe_fil <- overall_exp[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
oe_fil <- oe_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
oe_fil <- as.matrix(round(oe_fil))
# readr::write_rds(oe_fil, "/cluster/home/yjliu_jh/projects/leu_j/data/readcounts.rds")

# construct dds data and vst() to draw PCA plot
# first use smaller set to see the patterns


batch_s <- batch[batch$dataset %in% names(table(batch$dataset)[table(batch$dataset) > 100]), ]
oe_fil_s <- oe_fil[, colnames(oe_fil) %in% batch_s$sample_id]
ddss <- DESeqDataSetFromMatrix(countData = oe_fil_s,
                               colData = batch_s,
                               design = ~ dataset)
#dds <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/readcounts_dds.rds")
vsd <- varianceStabilizingTransformation(ddss)

# check PCA before batch correction
pcadata <- DESeq2::plotPCA(vsd, intgroup = c("dataset"), returnData = TRUE)
pca_before <- ggplot(pcadata, aes(x = PC1, y = PC2, color = group)) + geom_point(alpha = 1/3, size = 1.2)


# use sva::ComBat_seq to correct batch effects  ## some said it's better than limma?
batch_s$dataset <- as.character(batch_s$dataset)
adjusted_s <- ComBat_seq(counts = oe_fil_s, batch = batch_s$dataset)

ddss_adj <- DESeqDataSetFromMatrix(countData = adjusted_s,
                                   colData = batch_s,
                                   design = ~ dataset)
vsd_adj <- varianceStabilizingTransformation(ddss_adj)
pcadata2 <- DESeq2::plotPCA(vsd_adj, intgroup = c("dataset"), returnData = TRUE)
pca_after <- ggplot(pcadata2, aes(x = PC1, y = PC2, color = group)) + geom_point(alpha = 1/3, size = 1.2)

readr::write_rds(pcadata, "/cluster/home/yjliu_jh/projects/leu_j/data/temp_pca1.rds")
readr::write_rds(pcadata2, "/cluster/home/yjliu_jh/projects/leu_j/data/temp_pca2.rds")

pcadata <- read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/temp_pca1.rds")
pcadata2 <- read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/temp_pca2.rds")





png("pcabefore.png", width = 720, height = 720)
pca_before
dev.off()

png("pcaafter.png", width = 720, height = 720)
pca_after
dev.off()

# we confirmed that the batch correction works well

# before clustering results are out, 
# check differential gene expression between groups / inter-group comparison with similar characteristics? 
# and also before/after correction?


ddss <- DESeqDataSetFromMatrix(countData = oe_fil,
                               colData = batch,
                               design = ~ dataset)







# use cola
overall_exp <- read_csv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/overall_exp.csv")
oe_fil <- overall_exp[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
oe_fil <- oe_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
oe_mat <- adjust_matrix(as.matrix(oe_fil))
# 
rl = run_all_consensus_partition_methods(oe_mat, max_k = 12,
                                         top_value_method = c("SD", "MAD", "CV", "ATC"),
                                         partition_method = c("hclust", "kmeans", "skmeans", "mclust", "pam"),
                                         cores = 20)
cola_report(rl, output_dir = "/cluster/home/yjliu_jh/projects/leu_j/output/")


# (use cola on the small subset of leukemia)
count_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human_counts.csv"))
exp_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human.csv"))
s_fil <- exp_s[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
s_fil <- s_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
s_fil <- adjust_matrix(as.matrix(s_fil))

ssd <- transform(s_fil, SD=apply(s_fil, 1, sd, na.rm = TRUE))
s_fil2 <- s_fil[order(ssd$SD, decreasing = T), ][1:3000, ]
# stuck at CV-skmeans  500 rows are still too hard for the algorithm. 250 top rows (from 3000 rows) are at regular speed


# try again:
rls = run_all_consensus_partition_methods(s_fil2, max_k = 12,
                                         top_value_method = c("SD", "MAD", "CV", "ATC"),
                                         partition_method = c("hclust", "kmeans", "skmeans", "mclust", "pam"),
                                         cores = 20)
readr::write_rds(rls, "/cluster/home/yjliu_jh/projects/leu_j/data/rls_exp.rds")
# cola_report(rls, output_dir = "/cluster/home/yjliu_jh/projects/leu_j/output/")


suggest_best_k(rls)


png("what.png", width = 1280, height = 1280)
collect_plots(rls, k = 10)
dev.off()

png("what12.png", width = 1280, height = 1280)
collect_plots(rls, k = 12)
dev.off()



png("whatisthis.png", width = 1280, height = 1280)
collect_stats(rls, k = 10, layout_nrow = 3)
dev.off()


ress <- rls["SD", "skmeans"]
png("whatheatmap.png", width = 1280, height = 1280)
top_rows_heatmap(ress, k = 9)
dev.off()

rls_sd <- rls["SD", ]

png("classes8_sd.png", width = 1280, height = 720)
collect_classes(rls_sd, k = 8)
dev.off()




png("whatheatmapall.png", width = 1280, height = 1280)
top_rows_heatmap(rls, k = 9)
dev.off()

png("classes.png", width = 1280, height = 1280)
collect_classes(rls, k = 9)
dev.off()

png("classes10.png", width = 1280, height = 1280)
collect_classes(rls, k = 10)
dev.off()

png("classes3.png", width = 1280, height = 1280)
collect_classes(rls, k = 3)
dev.off()

png("classes4.png", width = 1280, height = 1280)
collect_classes(rls, k = 4)
dev.off()

png("classes8.png", width = 1280, height = 1280)
collect_classes(rls, k = 8)
dev.off()

png("classes6.png", width = 1280, height = 1280)
collect_classes(rls, k = 6)
dev.off()


png("whatoverlap.png", width = 1280, height = 1280)
top_rows_overlap(rls, k = 10)
dev.off()


cola_report(rls, output_dir = "/cluster/home/yjliu_jh/projects/leu_j/output/rls")

res = consensus_partition(oe_mat, top_value_method = "ATC", partition_method = "kmeans", ...)


rls1 = run_all_consensus_partition_methods(mat1, max_k = 12,
                                           top_value_method = c("SD", "MAD", "ATC"),
                                           partition_method = c("hclust", "kmeans", "skmeans", "mclust", "pam"),
                                           cores = 20)
# select methods for large dataset
osd <- transform(oe_mat, SD=apply(oe_mat, 1, sd, na.rm = TRUE))
oe_mat2 <- oe_mat[order(osd$SD, decreasing = T), ][1:3000, ]
rlb = run_all_consensus_partition_methods(oe_mat2, max_k = 12,
                                          top_value_method = c("SD", "ATC"),
                                          partition_method = c("kmeans", "skmeans"),
                                          cores = 40)
readr::write_rds(rlb, "/cluster/home/yjliu_jh/projects/leu_j/data/rlb_exp_3000.rds")

oe_mat2 <- oe_mat[order(osd$SD, decreasing = T), ][1:5000, ]
rlb = run_all_consensus_partition_methods(oe_mat2, max_k = 12,
                                          top_value_method = c("SD", "ATC"),
                                          partition_method = c("kmeans", "skmeans"),
                                          cores = 40)
readr::write_rds(rlb, "/cluster/home/yjliu_jh/projects/leu_j/data/rlb_exp_5000.rds")




readr::write_rds(oe_mat, "/cluster/home/yjliu_jh/projects/leu_j/data/oe_mat.rds")



count1 <- read_csv("/cluster/home/jhuang/projects/leukemia/analysis/linxiangjie/human/rnaseq/exp/tables/linxiangjie_human_counts.csv")
counts_fil <- count1[count1$gene_name %in% heatgenes$gene_name, ]
counts_fil <- counts_fil[, -1] %>% as.data.frame() %>% unique() %>% remove_rownames() %>% column_to_rownames(var = "gene_name") 
counts_fil <- round(counts_fil[, -1])
counts_fil <- as.matrix(counts_fil)

grouping1 <- ct2[, 1:2]
# detailed groupings should be curated after
grouping1 <- grouping1[match(colnames(counts_fil), grouping1$sample_id), ]
grouping1$subgroups <- as.factor(grouping1$subgroups)

dds1 <- DESeqDataSetFromMatrix(countData = counts_fil,
                               colData = grouping1,
                               design = ~ subgroups)
keep <- rowSums(counts(dds1) >= 2) > ceiling(0.2 * ncol(counts_fil))
dds1 <- dds1[keep,]
keep2 <- rowSums(counts(dds1) == 0) < ceiling(0.5 * ncol(counts_fil))
dds1 <- dds1[keep2,]
suppressMessages(dds <- DESeq(dds1)) 





