# cr: flora
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# read counts data for all datasets
datasets <- dir_ls("/cluster/home/jhuang/projects/leukemia/analysis", type = "directory") |>
  path_file()
datasets <- datasets[!datasets %in% c("bak", "diagnosis")]
i=1
mat <- vector("list", length(datasets))
merge_func <- function(x,y){left_join(x,y, by = c("gene_id", "gene_name"))}
for (dataset in datasets){
  csv <- glue("/cluster/home/jhuang/projects/leukemia/analysis/{dataset}/human/rnaseq/exp/tables/{dataset}_human_counts.csv")
  if (file.exists(csv)){
    exp <- read_csv(csv)
    if (csv == "/cluster/home/jhuang/projects/leukemia/analysis/ebiomedicine/human/rnaseq/exp/tables/ebiomedicine_human.csv"){
      colnames(exp)[-c(1,2)] <- glue("ebio_{colnames(exp)[-c(1,2)]}") #change ebiomedicine sample names
    }
    mat[[i]]  <- exp
  }
  i=i+1
}

# test if empty
data_n <- numeric()
for (i in 1:length(mat)){
  data_n[i] <- sum(ncol(mat[[i]]))
}

# output
overall_exp <- mat[-which(data_n == 0)] %>% Reduce(f = merge_func) %>% dplyr::select(-starts_with("gene_name.")) %>%
  relocate(gene_id, gene_name, everything())
#sum(colnames(overall_exp) %in% c("gene_id", "gene_name"))  ## 2
overall_exp %>% write_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_exp_count_short_221216.csv")


# read sampleinfo
overall_exp <- readr::read_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_exp_count_short_221216.csv")
sinfo <- readr::read_tsv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/sampleinfo.txt")
batch <- sinfo[, c("sample_id", "dataset")]
batch <- batch[match(colnames(overall_exp)[-(1:2)], batch$sample_id), ]
batch$dataset[is.na(batch$dataset)] <- "Unknown"
batch$sample_id[is.na(batch$sample_id)] <- setdiff(colnames(overall_exp)[-(1:2)], batch$sample_id)
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
oe_fil %>% readr::write_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_count_221216.rds")

# select 3000 most variable genes to plotPCA and ensure the significance of correction
oe_mat <- cola::adjust_matrix(oe_fil)
osd <- transform(oe_mat, SD=apply(oe_mat, 1, sd, na.rm = TRUE))
oe_fil_s <- round(oe_mat[order(osd$SD, decreasing = T), ][1:3000, ])

# original
ddss <- DESeqDataSetFromMatrix(countData = oe_fil_s,
                               colData = batch,
                               design= ~ dataset)
vsd <- varianceStabilizingTransformation(ddss)
pcadata <- DESeq2::plotPCA(vsd, intgroup = c("dataset"), returnData = TRUE)

# adjusted
adjusted_s <- sva::ComBat_seq(counts = oe_fil_s, batch = batch$dataset)
ddss_adj <- DESeqDataSetFromMatrix(countData = adjusted_s,
                                   colData = batch,
                                   design= ~ dataset)
vsd_adj <- varianceStabilizingTransformation(ddss_adj)
pcadata2 <- DESeq2::plotPCA(vsd_adj, intgroup = c("dataset"), returnData = TRUE)

readr::write_rds(vsd_adj, "/cluster/home/yjliu_jh/projects/leu_j/data/vsd_adj_221216.rds")


# use the adjusted exp data for clustering

DESeq2::plotPCA(vsd[, colnames(vsd) %in% bigsamples], intgroup = c("dataset"))
DESeq2::plotPCA(vsd_adj[, colnames(vsd_adj) %in% bigsamples], intgroup = c("dataset"))





vsd_adj <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/vsd_adj_221216.rds")

bigsets <- names(sort(table(batch$dataset), decreasing = T)[sort(table(batch$dataset), decreasing = T) > 200])
bigsamples <- batch$sample_id[batch$dataset %in% bigsets]

exp_vst_adj <- assay(vsd_adj)



sinfo2 <- left_join(batch, sinfo[, 1:4])
sinfo2[is.na(sinfo2)] <- "Unknown"

ha <- HeatmapAnnotation(
  df = as.data.frame(sinfo2),
  show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE)
  #gp = gpar(border = "gray", col = "white", lwd = 0.2),
  #height = unit(0.4, "cm"), simple_anno_size_adjust = TRUE,

)

col = list(
  fusions = fusions_col
)



cl_genes <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/chen_aml_clustering_genes.xlsx", skip = 1)
cl_genes <- cl_genes[[1]]


rm_genes <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/target_pattern_extracted_genes.xlsx", sheet = 2, col_names = F)
rm_genes <- rm_genes[[1]]



length(intersect(cl_genes, rownames(exp_vst_adj)[row_order_h][1:1000]))


vadj <- matrixStats::rowVars(exp_vst_adj)
exp_vst_adj_var <- cbind(exp_vst_adj, vadj)
var_adj <- as.data.frame(exp_vst_adj_var) %>% dplyr::arrange(desc(vadj))
var_adj <- var_adj[rownames(var_adj) %notin% rm_genes, ]
vst_mat <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 

length(intersect(cl_genes, rownames(var_adj)[1:1000]))





h <- ComplexHeatmap::Heatmap(vst_mat[1:700, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = FALSE, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = F,
                             column_title = "test_heatmap",
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha)


ht <- draw(h)
png("/cluster/home/yjliu_jh/share/cluster_sub1.png", width = 1800, height = 1500, units = "px")
ht
dev.off()





h2 <- ComplexHeatmap::Heatmap(vst_mat[intersect(cl_genes, rownames(vst_mat)), ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = FALSE, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = F,
                             column_title = "test_heatmap",
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha)


ht2 <- draw(h2)
png("/cluster/home/yjliu_jh/share/cluster_sub2.png", width = 1800, height = 1500, units = "px")
ht2
dev.off()












row_order_h <- row_order(ht)
column_order_h <- column_order(ht)



temp_mat <- fpkm_mat[1:1000, ][row_order_h, ]

odgenes <- rownames(fpkm_mat[1:1000, ][row_order_h, ])
odsamples <- colnames(fpkm_mat[, column_order_h])






mat_row_od <- exp_vst_adj[1:1000, ][row_order_h, ]
mat_sub1 <- mat_row_od[1:150, ]

h_sub1 <- ComplexHeatmap::Heatmap(mat_sub1, 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = FALSE, 
                             cluster_rows = F,
                             clustering_method_columns = 'ward.D2',
                             show_column_dend = F,
                             column_title = "what",
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha)

png("/cluster/home/yjliu_jh/share/cluster_sub1.png", width = 1800, height = 1500, units = "px")
draw(h_sub1)
dev.off()






ordered_mat <- exp_vst_adj[1:1000, ][row_order_h, column_order_h]

h@column_order


pdf("/cluster/home/yjliu_jh/share/clusterpdf.pdf", width = 15, height = 12)
draw(h)
dev.off()



png("/cluster/home/yjliu_jh/share/clusterng.png", width = 1800, height = 1500, units = "px")
draw(h)
dev.off()





# ====== fpkm part ===========

i=1
mat <- vector("list", length(datasets))
merge_func <- function(x,y){left_join(x,y, by = c("gene_id", "gene_name"))}
for (dataset in datasets){
  csv <- glue("/cluster/home/jhuang/projects/leukemia/analysis/{dataset}/human/rnaseq/exp/tables/{dataset}_human.csv")
  if (file.exists(csv)){
    exp <- read_csv(csv)
    if (csv == "/cluster/home/jhuang/projects/leukemia/analysis/ebiomedicine/human/rnaseq/exp/tables/ebiomedicine_human.csv"){
      colnames(exp)[-c(1,2)] <- glue("ebio_{colnames(exp)[-c(1,2)]}") #change ebiomedicine sample names
    }
    mat[[i]]  <- exp
  }
  i=i+1
}

# output
overall_fpkm <- mat[-which(data_n == 0)] %>% Reduce(f = merge_func) %>% dplyr::select(-starts_with("gene_name.")) %>%
  relocate(gene_id, gene_name, everything())
#sum(colnames(overall_exp) %in% c("gene_id", "gene_name"))  ## 2
overall_fpkm %>% write_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_short_221216.csv")


of_fil <- overall_fpkm[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
of_fil <- of_fil %>% dplyr::select(-c(gene_id, locus_group, entrez_id)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
of_fil <- as.matrix(of_fil)
of_fil %>% readr::write_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_221216.rds")


of_fil <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/overall_fpkm_221216.rds")


var <- matrixStats::rowVars(of_fil)
fpkm_var <- cbind(of_fil, var)
fpkm_ordered <- as.data.frame(fpkm_var) %>% dplyr::arrange(desc(var))
fpkm_ordered <- fpkm_ordered[rownames(fpkm_ordered) %notin% rm_genes, ]
fpkm_mat <- fpkm_ordered %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 

length(intersect(cl_genes, rownames(fpkm_ordered)[1:1000]))


ha <- HeatmapAnnotation(
  df = as.data.frame(sinfo2),
  show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE)
)

h <- ComplexHeatmap::Heatmap(fpkm_mat[1:1000, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = FALSE, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = F,
                             column_title = "refiltered_heatmap",
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha)

png("/cluster/home/yjliu_jh/share/cluster_newfilter.png", width = 1800, height = 1500, units = "px")
draw(h)
dev.off()



hx <- ComplexHeatmap::Heatmap(fpkm_mat[intersect(rownames(fpkm_mat), cl_genes), ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = FALSE, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = F,
                             column_title = "refiltered_heatmap",
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha)

htx <- draw(hx)
png("/cluster/home/yjliu_jh/share/cluster_newfilter2.png", width = 1800, height = 1500, units = "px")
htx
dev.off()




h02 <- ComplexHeatmap::Heatmap(fpkm_mat[1:1000, ][row_order_h, ][500:1000, ], 
                              name = " ",
                              show_row_names = FALSE, 
                              use_raster = T,
                              show_column_names = FALSE, 
                              clustering_method_columns = 'ward.D2',
                              clustering_method_rows = 'ward.D2',
                              show_row_dend = F,
                              show_column_dend = F,
                              column_title = "refiltered_heatmap",
                              heatmap_legend_param = list(
                                legend_direction = "horizontal", 
                                legend_width = unit(3, "cm")),
                              top_annotation = ha)

ht02 <- draw(h02)
png("/cluster/home/yjliu_jh/share/cluster_newfilter02.png", width = 1800, height = 1500, units = "px")
ht02
dev.off()











sinfo <- readr::read_tsv("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/sampleinfo.txt")
batch <- sinfo[, c("sample_id", "dataset")]
batch <- batch[match(colnames(of_fil), batch$sample_id), ]
batch$dataset[is.na(batch$dataset)] <- "Unknown"
batch$sample_id[is.na(batch$sample_id)] <- setdiff(colnames(of_fil), batch$sample_id)
batch$dataset <- as.factor(batch$dataset)

sinfo2 <- left_join(batch, sinfo[, 1:4])
sinfo2[is.na(sinfo2)] <- "Unknown"
















# ======= clustering using overall data ======= 

# construct dds object and get DESeq-normalized matrix
ddso <- DESeqDataSetFromMatrix(countData = oe_fil,
                               colData = batch,
                               design = ~ class) 

#ddso <- estimateSizeFactors(ddso)
#normalized_counts1 <- counts(ddso, normalized=TRUE)
readr::write_rds(ddso, "/cluster/home/yjliu_jh/projects/leu_j/data/dds_overall.rds")
readr::write_rds(normalizedExp1, "/cluster/home/yjliu_jh/projects/leu_j/data/exp_norm_all.rds")
vsd1 <- vst(ddso, blind = FALSE)
normalizedExp1 <- assay(vsd1)


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








overall_exp <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/exp_norm_all.rds")
oe_mat <- adjust_matrix(overall_exp)
osd <- transform(oe_mat, SD=apply(oe_mat, 1, sd, na.rm = TRUE))
oe_mat2 <- oe_mat[order(osd$SD, decreasing = T), ][1:3000, ]



names(mat) <- datasets
bad_ind <- (1:length(colnames(of_fil)))[colnames(of_fil) %in% badsamples]
bad_data <- data.frame(sample_id = colnames(of_fil)[bad_ind])
data_info <- character()
data_info2 <- character()
for(i in 1:length(datasets)){
  #tempcolni <- colnames(mat[[i]])[-(1:2)]
  tempdataset <- rep(datasets[i], length(tempcolni))
  #data_info <- c(data_info, tempcolni)
  data_info2 <- c(data_info2, tempdataset)
}
data_infoall <- data.frame(sample_id = colnames(of_fil), dataset = data_info2)



tempin <- left_join(bad_data, data_infoall)

bad_mat <- of_fil[, badsamples]



base_path <- "/cluster/home/jhuang/projects/leukemia/analysis/pnas_tall/human/rnaseq/salmon"
samples_to_rm <- setdiff(paste0("A", 1:61), c("A8", "A9"))
rm_paths <- paste0(base_path, "/", samples_to_rm, "/")
unlink(rm_paths, recursive = T)

