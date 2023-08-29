pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq",
          "DESeq2", "umap", "Rtsne", "FactoMineR", "factoextra", "matrixStats")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
} 

# --- after apll_analysis_02 ---
# - dat: AML exp matrix  counts_dat: AML counts matrix  meta_plot5: annotation data 
# - anno_col: colors for subtype annotations

# load counts and exp matrices
# dat <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/aml_dat.rds")
# counts_dat <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/aml_dat_counts.rds")

# subset dat
small_dat <- dat[, rownames(meta_plot5)[meta_plot5$datasets %in% c("suzhou", "beataml", "aplsz") | meta_plot5$pheno != "AML"]]
small_counts_dat <- counts_dat[, colnames(small_dat)]
meta_small <- meta_plot5[colnames(small_dat), 1:8]
meta_small$sample_id <- rownames(meta_small)

smallest_dat <- dat[, rownames(meta_plot5)[meta_plot5$datasets %in% c("suzhou", "aplsz") | meta_plot5$pheno %notin% c("AML", "APL")]]
meta_smallest <- meta_small[colnames(smallest_dat), ]
counts_smallest <- counts_dat[, colnames(smallest_dat)]
counts_smallest <- counts_smallest[-which(apply(counts_smallest, 1, var) == 0), ]

more_candids <- c(rownames(meta_plot2[meta_plot2$pheno %in% "APL" & meta_plot2$dataset %in% "HRA002693", ]))
counts_4g <- counts_dat[, c(colnames(smallest_dat), more_candids)]
counts_4g <- counts_4g[-which(apply(counts_4g, 1, var) == 0), ]
exp_4g_dat <- dat[, colnames(counts_4g)]
meta_4g <- meta_small[colnames(counts_4g), ]
readr::write_rds(meta_4g, "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/meta_small_used.rds")

# output dat and annotation for Dr. Wen
jhtools::list2excel(list(meta = meta_4g), "~/projects/leu_j/output/leu/aml/metadata_apll_4groups.xlsx")
readr::write_tsv(as.data.frame(counts_4g), "~/projects/leu_j/output/leu/aml/counts_apll_4groups.tsv")

# UMAP
exp_fil_small <- small_dat %>% dplyr::arrange(desc(rowVars(as.matrix(small_dat))))
exp_fil_small <- exp_fil_small[1:1500, ]
umap_results <- umap::umap(t(exp_fil_small))
umap_plot_df <- data.frame(umap_results$layout) %>% rownames_to_column("sample_id") %>% left_join(meta_small)
p_umap <- ggplot(umap_plot_df, aes(x = X1, y = X2, color = pheno)) + geom_point() + theme_classic() +
  xlab("UMAP_1") + ylab("UMAP_2") + scale_color_manual(values = anno_col, name = "groups") + 
  theme(legend.position = c(0.12, 0.85))
ggsave("~/projects/leu_j/output/leu/aml/umap_small_01.pdf", p_umap, width = 6, height = 6)

exp_fil_smallest <- smallest_dat %>% dplyr::arrange(desc(rowVars(as.matrix(smallest_dat))))
exp_fil_smallest <- exp_fil_smallest[1:2500, ]
umap_results_2 <- umap::umap(t(exp_fil_smallest))
umap_plot_df2 <- data.frame(umap_results_2$layout) %>% rownames_to_column("sample_id") %>% left_join(meta_smallest)
p_umap <- ggplot(umap_plot_df2, aes(x = X1, y = X2, color = anno)) + geom_point() + theme_classic() +
  xlab("UMAP_1") + ylab("UMAP_2") + scale_color_manual(values = anno_col, name = "groups") + 
  theme(legend.position = c(0.83, 0.83))
ggsave("~/projects/leu_j/output/leu/aml/umap_small_03_detailed2500.pdf", p_umap, width = 5, height = 5)


exp_fil_4g <- exp_4g_dat %>% dplyr::arrange(desc(rowVars(as.matrix(exp_4g_dat))))
exp_fil_4g <- exp_fil_4g[1:2500, ]
umap_results_15 <- umap::umap(t(exp_fil_4g)) 
umap_plot_df15 <- data.frame(umap_results_15$layout) %>% rownames_to_column("sample_id") %>% left_join(meta_4g)
p_umap <- ggplot(umap_plot_df15, aes(x = X1, y = X2, color = anno)) + geom_point() + theme_classic() +
  xlab("UMAP_1") + ylab("UMAP_2") + scale_color_manual(values = anno_col, name = "groups") + 
  theme(legend.position = c(0.16, 0.83))
ggsave("~/projects/leu_j/output/leu/aml/umap_small_04_detailed2500.pdf", p_umap, width = 5, height = 5)


aml_anno <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx")
aml_anno <- left_join(aml_anno, as.data.frame(metadata[, c("sample_id", "subgroups_230712")]))
# temp_id <- umap_plot_df15$sample_id[umap_plot_df15$X2 > 6]


# defb1 krt7 krt19 anxa4 




# PCA
pca_results <- prcomp(t(smallest_dat), scale = TRUE)
p_pca <- fviz_pca_ind(pca_results, geom = "point", col.ind = meta_smallest$pheno, addEllipses = TRUE, ellipse.level = 0.95) +
  labs(title = "PCA", x = "PC1", y = "PC2") + 
  theme_classic() + scale_shape_manual(values = c(16, 16, 16, 16)) +
  scale_color_manual(values = anno_4)
ggsave("~/projects/leu_j/output/leu/aml/pca_small_03.pdf", p_pca, width = 5, height = 4)


pca_results <- prcomp(t(exp_4g_dat), scale = TRUE)
p_pca <- fviz_pca_ind(pca_results, geom = "point", col.ind = meta_4g$pheno, addEllipses = TRUE, ellipse.level = 0.95) +
  labs(title = "PCA", x = "PC1", y = "PC2") + 
  theme_classic() + scale_shape_manual(values = c(16, 16, 16, 16)) +
  scale_color_manual(values = anno_4)
ggsave("~/projects/leu_j/output/leu/aml/pca_small_04.pdf", p_pca, width = 5, height = 4)



# DESeq2 PCA
dds <- DESeqDataSetFromMatrix(countData = counts_smallest, 
                              colData = meta_smallest, 
                              design = ~ pheno)
keep <- rowSums(counts(dds) >= 2) > ceiling(0.2 * ncol(counts_smallest))
dds <- dds[keep,]
keep2 <- rowSums(counts(dds) == 0) < ceiling(0.5 * ncol(counts_smallest))
dds <- dds[keep2,]
suppressMessages(dds <- DESeq(dds)) 
vst <- vst(dds)
pca_deseq2 <- plotPCA(vst, intgroup = "pheno")
ggsave("~/projects/leu_j/output/leu/aml/pca_small_03_deseq2_vst.pdf", pca_deseq2, width = 6, height = 6)




dds <- DESeqDataSetFromMatrix(countData = counts_4g, 
                              colData = meta_4g, 
                              design = ~ pheno)
keep <- rowSums(counts(dds) >= 2) > ceiling(0.2 * ncol(counts_4g))
dds <- dds[keep,]
keep2 <- rowSums(counts(dds) == 0) < ceiling(0.5 * ncol(counts_4g))
dds <- dds[keep2,]
suppressMessages(dds <- DESeq(dds)) 
vst <- vst(dds)
pca_deseq2 <- plotPCA(vst, intgroup = "pheno")
ggsave("~/projects/leu_j/output/leu/aml/pca_small_04_deseq2_vst.pdf", pca_deseq2, width = 6, height = 6)









# Tsne
tSNE_fit <- Rtsne(t(smallest_dat), perplexity = floor((ncol(smallest_dat) - 1) / 3), dims = 2)
tSNE_df <- data.frame(tSNE_fit$Y)
colnames(tSNE_df) <- c("tSNE1", "tSNE2")
tSNE_df$sample_id <- colnames(smallest_dat)
tSNE_df <- left_join(tSNE_df, meta_small)

p_tsne <- tSNE_df %>% ggplot(aes(x = tSNE1, y = tSNE2, color = pheno)) +
  geom_point() + theme_classic() + scale_color_manual(values = anno_4, name = "") + 
  theme(legend.position = c(0.3, 0.8))
ggsave("~/projects/leu_j/output/leu/aml/tsne_small_03_deseq2_vst.pdf", p_tsne, width = 5, height = 5)




tSNE_fit <- Rtsne(t(exp_4g_dat), perplexity = floor((ncol(exp_4g_dat) - 1) / 3), dims = 2)
tSNE_df <- data.frame(tSNE_fit$Y)
colnames(tSNE_df) <- c("tSNE1", "tSNE2")
tSNE_df$sample_id <- colnames(exp_4g_dat)
tSNE_df <- left_join(tSNE_df, meta_4g)

p_tsne <- tSNE_df %>% ggplot(aes(x = tSNE1, y = tSNE2, color = anno)) +
  geom_point() + theme_classic() + scale_color_manual(values = anno_col, name = "") + 
  theme(legend.position = c(0.85, 0.84))
ggsave("~/projects/leu_j/output/leu/aml/tsne_small_04_detailed.pdf", p_tsne, width = 5, height = 5)









# go enrichment and gsea

normalizedExp1 <- assay(vst)
# corResult <- cor(normalizedExp1)
# pdf("~/projects/leu_j/output/leu/aml/tempcor001.pdf", width = 33, height = 33)
# corrplot(corResult, is.corr = F, tl.col = "black",
#          method = "color", col = colorRampPalette(c("blue","yellow","red"))(200))
# dev.off()


#=======================================================================
# Data preparation


meta_aapl <- meta_4g
meta_aapl$group <- ifelse(meta_aapl$pheno %in% "Atypical_APL", "AAPL", "not-AAPL")
meta_apll <- meta_4g
meta_apll$group <- ifelse(meta_apll$pheno %in% "APL-like", "APLL", "not-APLL")




dds <- DESeqDataSetFromMatrix(countData = counts_4g, 
                              colData = meta_aapl, 
                              design = ~ group)
keep <- rowSums(counts(dds) >= 2) > ceiling(0.2 * ncol(counts_4g))
dds <- dds[keep,]
keep2 <- rowSums(counts(dds) == 0) < ceiling(0.5 * ncol(counts_4g))
dds <- dds[keep2,]
suppressMessages(dds_aapl <- DESeq(dds)) 




dds <- DESeqDataSetFromMatrix(countData = counts_4g, 
                              colData = meta_apll, 
                              design = ~ group)
keep <- rowSums(counts(dds) >= 2) > ceiling(0.2 * ncol(counts_4g))
dds <- dds[keep,]
keep2 <- rowSums(counts(dds) == 0) < ceiling(0.5 * ncol(counts_4g))
dds <- dds[keep2,]
suppressMessages(dds_apll <- DESeq(dds)) 


dds_aapl$condition <- colData(dds_aapl)$group
dds_aapl$condition <- relevel(dds_aapl$condition, ref = "not-AAPL")

dds_apll$condition <- colData(dds_apll)$group
dds_apll$condition <- relevel(dds_apll$condition, ref = "not-APLL")




GSEAInput_a <- data.frame(rownames(results(dds_aapl)), results(dds_aapl)$log2FoldChange)
colnames(GSEAInput_a) <- c("Symbol", "logFC")
GSEA2 <- GSEAInput_a$logFC
names(GSEA2) <- as.character(GSEAInput_a$Symbol)
GSEA2 <- sort(GSEA2, decreasing = TRUE)
gseGO_res_a <- gseGO(GSEA2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                     minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.2)


ego <- pairwise_termsim(gseGO_res_a, method="Wang", semData = godata('org.Hs.eg.db', ont="BP"))


GSEAInput_b <- data.frame(rownames(results(dds_apll)), results(dds_apll)$log2FoldChange)
colnames(GSEAInput_b) <- c("Symbol", "logFC")
GSEA3 <- GSEAInput_b$logFC
names(GSEA3) <- as.character(GSEAInput_b$Symbol)
GSEA3 <- sort(GSEA3, decreasing = TRUE)
gseGO_res_b <- gseGO(GSEA3, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                     minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.2)





pdf("~/projects/leu_j/output/leu/aml/gsea_dot_aapl_01.pdf", width = 8, height = 9)
dotplot(gseGO_res_a, showCategory = 10)
dev.off()



pdf("~/projects/leu_j/output/leu/aml/gsea_emap_aapl_01.pdf", width = 8, height = 9)
emapplot(ego, showCategory = 10)
dev.off()


pdf("~/projects/leu_j/output/leu/aml/gsea_ridge_aapl_01.pdf", width = 8, height = 12)
ridgeplot(gseGO_res_a, 14) 
dev.off()


pdf("~/projects/leu_j/output/leu/aml/gsea_line_aapl_01.pdf", width = 8, height = 9)
gseaplot2(gseGO_res_a, 1:4) 
dev.off()




ego2 <- pairwise_termsim(gseGO_res_b, method="Wang", semData = godata('org.Hs.eg.db', ont="BP"))


pdf("~/projects/leu_j/output/leu/aml/gsea_dot_apll_01.pdf", width = 8, height = 9)
dotplot(gseGO_res_b, showCategory = 10)
dev.off()



pdf("~/projects/leu_j/output/leu/aml/gsea_emap_apll_01.pdf", width = 8, height = 9)
emapplot(ego2, showCategory = 10)
dev.off()


pdf("~/projects/leu_j/output/leu/aml/gsea_ridge_apll_01.pdf", width = 8, height = 12)
ridgeplot(gseGO_res_b, 14) 
dev.off()


pdf("~/projects/leu_j/output/leu/aml/gsea_line_apll_01.pdf", width = 8, height = 9)
gseaplot2(gseGO_res_b, 1:4) 
dev.off()









# function to order and add columns to the results


hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leukemia/data/hgnc_complete_set_2023-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id", "entrez_id")]
colnames(hugo_anno)[3] <- "gene_id"


add_col <- function(data){
  data$gene_name <- rownames(data)
  data$absFC <- abs(data$log2FoldChange)
  data <- data[order(data$pvalue), ]
  data2 <- left_join(as.data.frame(data), hugo_anno[, c("symbol", "entrez_id")], by = c("gene_name" = "symbol"))
  data$entrez_id <- data2$entrez_id
  data
}

# check functions
check_ego_gene <- function(data){
  ego <- enrichGO(gene = data, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego <- ego[order(ego$p.adjust), ]
  ego
}

check_ego_gene2 <- function(data){
  ego <- enrichGO(gene = data, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego
}




# auto-de functions (no need for QC in current dataset)
auto_de <- function(dds, compare){
  as.data.frame(add_col(results(dds, c("group", compare))))
}

test <- auto_de(dds_aapl, c("AAPL", "not-AAPL"))
test2 <- auto_de(dds_apll, c("APLL", "not-APLL"))

aapl_up <- rownames(test[test$log2FoldChange > 1, ])[1:200]
aapl_down <- rownames(test[test$log2FoldChange < -1, ])[1:200]
apll_up <- rownames(test2[test2$log2FoldChange > 1, ])[1:500]
apll_down <- rownames(test2[test2$log2FoldChange < -1, ])[1:500]


temp <- dotplot(enrichGO(gene = apll_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                         ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
pdf("~/projects/leu_j/output/leu/aml/go_dotplot_apll_up_02.pdf", height = 8, width = 6)
print(temp)
dev.off()














hsa_sets <- msigdbr(species = "Homo sapiens")
hsa_sets <- hsa_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
gseMsigDB_res1 <- GSEA(gene = GSEA1, TERM2GENE = hsa_sets)
gseaplot2(gseMsigDB_res1, 1:4) 

gseKEGG_res1 <- gseKEGG(geneList = GSEA1, organism = 'hsa', minGSSize = 120,
                        pvalueCutoff = 0.05, verbose = FALSE)
# no internet connection
















#=====================================================
gtest1 <- detailedgroups[, c("sample_id", "fusions")]
levels(gtest1$fusions) <- c("non-RARG","non-RARG","CPSF6-RARG","other-RARG","other-RARG","other-RARG")
gtest2 <- detailedgroups[, c("sample_id", "fusions")]
levels(gtest2$fusions) <- c("non-RARG","non-RARG","other-RARG","HNRNP-RARG","HNRNP-RARG","other-RARG")
gtest3 <- detailedgroups[, c("sample_id", "fusions")]
levels(gtest3$fusions) <- c("non-RARG","non-RARG","other-RARG","other-RARG","other-RARG","NUP98-RARG")

gtest1 <- unique(gtest1)
gtest2 <- unique(gtest2)
gtest3 <- unique(gtest3)


dds_CPSF6 <- DESeqDataSetFromMatrix(countData = all_exp_fil,
                                    colData = gtest1,
                                    design = ~ fusions)
dds_HNRNP <- DESeqDataSetFromMatrix(countData = all_exp_fil,
                                    colData = gtest2,
                                    design = ~ fusions)
dds_NUP98 <- DESeqDataSetFromMatrix(countData = all_exp_fil,
                                    colData = gtest3,
                                    design = ~ fusions)


keep <- rowSums(counts(dds_CPSF6) >= 2) > 9
dds_CPSF6 <- dds_CPSF6[keep,]
keep <- rowSums(counts(dds_HNRNP) >= 2) > 9
dds_HNRNP <- dds_HNRNP[keep,]
keep <- rowSums(counts(dds_NUP98) >= 2) > 9
dds_NUP98 <- dds_NUP98[keep,]



suppressMessages(dds_CPSF6 <- DESeq(dds_CPSF6)) 
suppressMessages(dds_HNRNP <- DESeq(dds_HNRNP)) 
suppressMessages(dds_NUP98 <- DESeq(dds_NUP98)) 


res_CPSF6 <- results(dds_CPSF6, c("fusions", "CPSF6-RARG", "non-RARG"))
GSEAInput_CPSF6 <- data.frame(rownames(res_CPSF6), res_CPSF6$log2FoldChange)
colnames(GSEAInput_CPSF6) <- c("Symbol", "logFC")
GSEA_CPSF6 <- GSEAInput_CPSF6$logFC
names(GSEA_CPSF6) = as.character(GSEAInput_CPSF6$Symbol)
GSEA_CPSF6 = sort(GSEA_CPSF6, decreasing = TRUE)

gseGO_res_CPSF6 <- gseGO(GSEA_CPSF6, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",
                         minGSSize = 5, maxGSSize = 1000, pvalueCutoff=0.2)


ego_CPSF6 <- pairwise_termsim(gseGO_res_CPSF6)
treeplot(ego_CPSF6) + ggtitle("CPSF6 vs non-RARG")


dotplot(gseGO_res_CPSF6, showCategory=10)
emapplot(ego_CPSF6, showCategory = 10)

ridgeplot(gseGO_res_CPSF6, 18) 
gseaplot2(gseGO_res_CPSF6, 1:4) 



