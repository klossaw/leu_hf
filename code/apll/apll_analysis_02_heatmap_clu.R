pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
} 

# --- after apll_analysis_02 ---
# - dat: AML exp matrix  dat_mat: scaled AML exp matrix  meta_plot: annotation data 
# - mutations_new: mutation annotation    fusions_new: fusion annotation

# set colors
# set colors
all_colors <- c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", 
                "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00",
                "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3",
                "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB",
                "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", "#008B8B",
                "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4",
                "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347",
                "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#DEB887")
col_fun <- circlize::colorRamp2(c(0, 10, 100), c("#FFEEEE", "#FFBBBB", "#FF0000"))
anno_col <- c("RARG" = "#FFE119", "RARB" = "#FF8534", "TTMV-RARA" = "#E25385",
              "pml-others" = "#000075", "RARA-others" = "#AAFFC3",
              "APL" = "#912EB4", "APL-like" = "#3CB44B", "AML" = "#96B2F3")
anno_4 <- c("APL" = "#912EB4", "AML" = "#96B2F3", "Atypical_APL" = "#F58231", "APL-like" = "#3CB44B")
anno_def <- ggsci::pal_nejm("default")(8)

# load itd data
itd_data <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/itd_data_hamlet.rds")
itd_data2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/itd_data_filt3r.rds")



# ========= functions ============

# add new columns to annotation data.frame (exist or not)
enrich_anno <- function(anno, data, extra = NULL, sample_col = "sample_id") {
  if (is.null(extra)) { extra <- names(data) }
  for (col in extra){
    anno[[col]] <- ifelse(anno[[sample_col]] %in% data[[col]], "exist", "non-exist")
  }
  anno
}

# set colors for columns
auto_anno <- function(data, cols, colors) {
  anno_list <- list()
  for (col in cols) {
    elements <- sort(unique(data[[col]]), decreasing = T)
    anno_list[[col]] <- setNames(colors[1:length(elements)], elements)
  }
  anno_list
}

# add new columns to annotation data.frame (based on annotation col)
enrich_anno2 <- function(anno, data, extra, sample_col = "sample_id", anno_col = "anno") {
  for (col in extra){
    new_col <- left_join(anno, data[[col]][, c(sample_col, anno_col)])
    anno[[col]] <- new_col[[anno_col]]
  }
  anno
}

# select hotspots for mutations
mut_hotspot <- function (gene, mut_data, col = "HGVSp_Short", hotspot = NULL, k = 10, sum = T, agg = F) {
  mut_gene <- mut_data[mut_data$Hugo_Symbol %in% gene, ] %>% as.data.frame()
  if (is.null(hotspot)) {
    hotspot <- names(sort(table(mut_gene[[col]]), decreasing = T))[1:k]
  } 
  mut_gene_fil <- mut_gene[mut_gene[[col]] %in% hotspot, ] 
  if (sum == T) {
    mut_gene_fil <- mut_gene_fil[, c("sample_id", col)] %>% unique() %>% group_by(sample_id) %>% 
      summarise(mutation = paste(.data[[col]], collapse = ", "))
    if (agg == T) { mut_gene_fil$mutation[grepl(",", mut_gene_fil$mutation)] <- "mixed" }
  } else { mut_gene_fil$mutation <- paste0(mut_gene_fil$Hugo_Symbol, "_", mut_gene_fil[[col]]) }
  mut_gene_fil
}



# select interested gene patterns
mut_genes <- c("KRAS", "WT1", "CEBPA", "DNMT3A", "CYB5R2", "ADAMTS5", "CEBPD", "COL6A5", "COBL", "SLC22A1")
fus_genes <- c("PML-RARA", "KMT2A-r", "RUNX1-RUNX1T1", "CBFB-MYH11", "NUP98-r")
mut_genes_sep <- c("NPM1", "IGH1", "IGH2", "TET2")
itd_genes <- c("FLT3-ITD", "KMT2A-PTD")
itd2_genes <- c("filt3r_FLT3-ITD")

# list construction
mutall_list <- list()
for (i in 1:length(mut_genes)) {
  mutall_list[[mut_genes[i]]] <- unique(mutations_new$sample_id[mutations_new$Hugo_Symbol %in% mut_genes[i]])
}

fusion_list <- list()
for (i in 1:length(fus_genes)) {
  fusion_list[[fus_genes[i]]] <- unique(fusions_new$sample_id[fusions_new$fusion_anno %in% fus_genes[i]])
}

itd_list <- list()
for (i in 1:length(itd_genes)) {
  itd_list[[itd_genes[i]]] <- unique(itd_data$sample_id[itd_data$itd %in% itd_genes[i]])
}

itd2_list <- list()
for (i in 1:length(itd2_genes)) {
  itd2_list[[itd2_genes[i]]] <- unique(itd_data2$sample_id)
}

# hotspot list construction
idh1_hotspots <- c("p.R132H", "p.R132C", "p.R132G", "p.R132S", "p.R132L", "p.R132V",
                   "p.R132I", "p.P33S", "p.R100Q")
idh2_hotspots <- c("p.R140Q", "p.R172K", "p.R172S", "p.R172M", "p.R172G", "p.R172W",
                   "p.R140W", "p.R140L", "p.R140G")
idh1_mut <- mut_hotspot("IDH1", mutations_new, hotspot = idh1_hotspots)
idh2_mut <- mut_hotspot("IDH2", mutations_new, hotspot = idh2_hotspots)
npm1_mut <- mut_hotspot("NPM1", mutations_new)
tet2_mut <- mut_hotspot("TET2", mutations_new, k = 7, agg = T)

mut_sep_list <- list(npm1_mut, idh1_mut, idh2_mut, tet2_mut)
names(mut_sep_list) <- mut_genes_sep


# annotation data for plot
meta_plot2 <- as.data.frame(meta_plot)[, c("gender", "age", "rna_type", "disease_type", "datasets", "anno", "call", "pheno")]
meta_plot2$sample_id <- rownames(meta_plot2)
meta_plot2$age <- as.numeric(meta_plot2$age)
meta_plot2$pheno[meta_plot2$pheno %in% "RARA-others"] <- "Atypical_APL"
meta_plot3x <- enrich_anno(meta_plot2, mutall_list, mut_genes)
meta_plot4x <- enrich_anno(meta_plot3x, itd_list, itd_genes)
meta_plot4x <- enrich_anno(meta_plot4x, itd2_list, itd2_genes)
meta_plot4x <- enrich_anno(meta_plot4x, fusion_list, fus_genes)
meta_plot5 <- enrich_anno2(meta_plot4x, mut_sep_list, mut_genes_sep, anno_col = "mutation")


construct_anno <- function(anno_plot) { ## contains global variables
  HeatmapAnnotation(
    df = as.data.frame(anno_plot[, -c(3, 4, 9)]),
    show_legend = c(TRUE, TRUE, FALSE, T, T, T, 
                    rep(F, length(mut_genes) + length(itd_genes) + 1 + length(fus_genes)),
                    rep(T, length(mut_genes_sep))),
    col = c(list(age = col_fun, 
                 anno = anno_col,
                 call = anno_col,
                 pheno = anno_4,
                 gender = setNames(anno_def[1:2], c("female", "male"))),
            auto_anno(anno_plot, c(mut_genes, itd_genes, itd2_genes, fus_genes), c("ivory", "magenta")),
            auto_anno(anno_plot, c(mut_genes_sep, "sub_groups", "subgroups_230712"), all_colors))
  )
}


# subset 1: include as much apl/apl-like samples and datasets from suzhou hospitals and the beatAML project
small_dat <- dat[, rownames(meta_plot5)[meta_plot5$datasets %in% c("suzhou", "beataml", "aplsz") | meta_plot5$pheno != "AML"]]
meta_small <- meta_plot5[colnames(small_dat), ]
ha_small <- construct_anno(meta_small)

# may also start here and load the functions
mutations_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/mutations_new.rds")
fusions_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/fusions_new.rds")
meta_plot5 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/meta_plot_aml.rds")
ha_small <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/ha_small.rds")
dat <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/aml_dat.rds")
counts_dat <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/rds/aml_dat_counts.rds")


setwd("/cluster/home/yjliu_jh/projects/leu_j/output")
quick_heatmap(small_dat, ha = ha_small, width = 28, height = 16, outdir = "leu/aml/apll", column_split = 26, top_var_percent = 0.95)

# mark outlier apls and apl-likes, and narrow the scope
cinfo_1 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/output/leu/aml/apll/sampleinfo_0.95.xlsx")
c1anno <- left_join(cinfo_1, data.frame(sample_id = rownames(meta_plot5), anno = meta_plot5$pheno))
c1anno$ct <- ifelse(c1anno$anno %notin% "AML", 1, 0)
perc_c1 <- c1anno %>% group_by(sub_groups) %>% summarise(perc = sum(ct) / n())
amb_aml_samples <- c1anno$sample_id[c1anno$anno %in% "AML" & c1anno$sub_groups %in% perc_c1$sub_groups[perc_c1$perc > 0.8]]
amb_apl_samples <- c1anno$sample_id[c1anno$anno %notin% "AML" & c1anno$sub_groups %in% perc_c1$sub_groups[perc_c1$perc < 0.2]]
mix_samples <- c1anno$sample_id[c1anno$sub_groups %in% perc_c1$sub_groups[perc_c1$perc > 0.3 & perc_c1$perc < 0.7]]

# subset 2: only include apl samples and samples from suzhou
smaller_dat <- dat[, rownames(meta_plot5)[meta_plot5$datasets %in% c("suzhou", "aplsz") | meta_plot5$pheno != "AML"]]
meta_smaller <- meta_plot5[colnames(smaller_dat), ]
ha_smaller <- construct_anno(meta_smaller)

quick_heatmap(smaller_dat, ha = ha_smaller, width = 28, height = 16, outdir = "leu/aml/apll_small", column_split = 13, top_var_percent = 0.95)


# subset 3: only include rare apl-like samples and samples mentioned from Dr. Wen in suzhou
smallest_dat <- dat[, rownames(meta_plot5)[meta_plot5$datasets %in% c("suzhou", "aplsz") | meta_plot5$pheno %notin% c("AML", "APL")]]
meta_smallest <- meta_plot5[colnames(smallest_dat), ]
ha_smallest <- construct_anno(meta_smallest)

quick_heatmap(smallest_dat, ha = ha_smallest, width = 28, height = 16, outdir = "leu/aml/apll_smallest", column_split = 8, top_var_percent = 0.95)


# subset 4: no aml 
min_dat <- smallest_dat[, meta_smallest$pheno %notin% "AML"]
meta_min <- meta_plot5[colnames(min_dat), ]
ha_min <- construct_anno(meta_min)

# modify to add colnames
quick_hm2 <- function(dat, top_var_percent = 0.9, outdir = "tall", ha = NULL,  
              heatmap_title = "leukemia", use_raster = F, column_split = 10, 
              height = 6, width = 16, cluster_columns = TRUE, show_col = FALSE) {
  dat_var <- apply(dat, 1, var)
  exp_plot <- data.matrix(dat[dat_var >= quantile(dat_var, probs = top_var_percent), ])
  hm <- jhHeatmap(as.matrix(exp_plot), glue("{outdir}/jhHeatmap_{top_var_percent}.pdf"), 
                  width = width, height = height, Colv = cluster_columns)
  pdf(glue("{outdir}/heatmap_{top_var_percent}.pdf"), height = height, 
      width = width)
  mat_scaled <- t(hm$carpet)[row.names(exp_plot), colnames(exp_plot)]
  ht <- ComplexHeatmap::Heatmap(mat_scaled, name = heatmap_title, 
                                use_raster = use_raster, col = circlize::colorRamp2(c(-1.8, 0, 1.8), c("blue", "white", "red")),
                                row_order = rev(hm$rowInd), column_order = hm$colInd, 
                                cluster_rows = rev(hm$rowDendrogram), cluster_columns = hm$colDendrogram,
                                show_row_names = FALSE, show_column_names = show_col,
                                column_names_side = "top", top_annotation = ha, column_dend_height = unit(1.5, "cm"),
                                column_split = column_split, column_gap = unit(1, "mm"), column_title = NULL
  )
  ht <- ComplexHeatmap::draw(ht, annotation_legend_side = "left", 
                             heatmap_legend_side = "right")
  dev.off()
}

quick_hm2(min_dat, ha = ha_min, width = 20, height = 18, outdir = "leu/aml/apll_min2", column_split = 6, top_var_percent = 0.95, show_col = T)


