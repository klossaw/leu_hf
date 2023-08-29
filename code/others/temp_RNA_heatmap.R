metadata <- read_tsv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/sampleinfo.txt")
s_totalRNA <- metadata$sample_id[metadata$rna_type == "total_RNA"]
s_mRNA <- metadata$sample_id[metadata$rna_type == "mRNA"]

overall_exp <- read_csv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/overall_exp.csv")

zt_genes <- read_csv("/cluster/home/ztao_jh/jianhuiRNA.csv")
zt_genes <- zt_genes$...1[1:50]

zt_predict <- readxl::read_excel("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/zt_predict.xlsx")
zt_predict <- zt_predict[, c("sample_id", "total_RNA_percent")]
zt_predict$predict <- ifelse(zt_predict$total_RNA_percent > 0.6, "total",
                             ifelse(zt_predict$total_RNA_percent < 0.4, "mRNA", "amb"))
zt_predict_sure <- zt_predict[zt_predict$predict != "amb", ]
sure_samples <- zt_predict_sure$sample_id
exp_mat <- overall_exp[overall_exp$gene_name %in% zt_genes, ]
exp_mat <- exp_mat[, sure_samples]
Var_gene_mat <- exp_mat %>% as.matrix() 
metadata <- metadata %>% left_join(zt_predict_sure) %>% na.omit()



data <- overall_exp %>% dplyr::filter(gene_name %in% zt_genes) 
data <- data[,c("gene_id", "gene_name", sure_samples)]
fpkm_mat <- data %>% as.data.frame() %>% dplyr::select(-c("gene_name", "gene_id"))
rownames(fpkm_mat) <- data$gene_name
Variations <- matrixStats::rowVars(as.matrix(fpkm_mat))
fpkm_with_var <- cbind(fpkm_mat, Variations)
Var_genes <- fpkm_with_var %>% dplyr::arrange(desc(Variations))
Var_gene_mat <- Var_genes %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 




col_dend = dendsort(hclust(dist(t(Var_gene_mat))))

ha <-  HeatmapAnnotation(dataset = metadata$dataset,
                         result = metadata$predict
)
h <- ComplexHeatmap::Heatmap(as.matrix(Var_gene_mat), name = "FPKM",
                             show_row_names = FALSE, use_raster = F, cluster_columns = F,
                             show_column_names = FALSE, top_annotation = ha)
od=column_order(h)
sample_order = data.frame(sample_id = colnames(Var_gene_mat)[od])
metadata <- left_join(sample_order, metadata, by = "sample_id")
meta1 <- metadata[metadata$predict %in% "mRNA", ]
meta2 <- metadata[metadata$predict %notin% "mRNA", ]
metadata <- bind_rows(meta1, meta2)
new_order <- metadata$sample_id
Var_gene_mat <- Var_gene_mat[, new_order]


metadata %>% write_tsv("htmap/RNA_type/overall_rna_type4.txt")
pdf("htmap/RNA_type/overall_rna_type4.pdf", width = 10, height = 6)
print(h)
dev.off()


