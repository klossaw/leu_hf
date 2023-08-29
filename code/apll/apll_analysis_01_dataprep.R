pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
} 

# main data
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
cld <- colData(rds) |> as.data.frame() %>% replace(., is.na(.), "")
fil <- cld$disease_type  %in% c("AML")
se <- rds[,fil]
metadata <- colData(se)
fpkm <- assay(se, "tpm_remove_rnatype") %>% as.data.frame() %>% dplyr::select(metadata$sample_id) # exp matrix
dat <- fpkm |> as.data.frame() |> tibble::rownames_to_column("gene_name") |>
  filter_black_heatmap_genes() |> tibble::remove_rownames() |> tibble::column_to_rownames("gene_name")

# add new data
setwd("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/exp")
new_exp <- read_csv("datExpBeforeSVA.csv")  ## by jhuanglabRNAseq::getexpSalmon
new_exp <- new_exp[, -1] %>% as.data.frame() %>% na.omit() %>%
  remove_rownames() %>% column_to_rownames(var = "gene_name")
dat <- cbind(dat, new_exp[rownames(dat), ])

# load mutation data
mutations_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/temp/all_mutations_0716.rds")
mut_int <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/apll_new_mutations.rds")
mut_apll_new <- mut_int[, colnames(mut_int) %in% c(colnames(mutations_new), "tool"), with = F]
temp <- mut_apll_new[, toString(tool), by = setdiff(names(mut_apll_new), "tool")]
colnames(temp)[5] <- "tools"
mutations_new <- rbind(mutations_new, as.data.frame(temp)[, colnames(mutations_new)])

# load fusion data 
fusions_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/fusions_new.rds")
fu_all <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/apll_new_fusions.rds")
fu_all$fusion_anno <- fu_all$fusion
fu_all$fusion_anno[grepl("KMT2A", fu_all$fusion_anno)] <- "KMT2A-r"
fu_all$fusion_anno[grepl("NUP98", fu_all$fusion_anno)] <- "NUP98-r"
fu_all$fusion_anno[grepl("MECOM", fu_all$fusion_anno)] <- "MECOM-r"
fu_all$fusion_anno[grepl("DUX4", fu_all$fusion_anno)] <- "DUX4-r"
fu_all$fusion_anno[grepl("PAX5", fu_all$fusion_anno)] <- "PAX5-r"
fu_all$fusion_anno[grepl("MLLT10", fu_all$fusion_anno)] <- "MLLT10-r"
fu_all$fusion_anno[grepl("CRLF2", fu_all$fusion_anno)] <- "CRLF2-r"
fu_all$fusion_anno[grepl("LMO1", fu_all$fusion_anno)] <- "LMO1-r"
fu_all$fusion_anno[grepl("LMO2", fu_all$fusion_anno)] <- "LMO2-r"
fusions_new <- rbind(fusions_new, fu_all)

# load other data
ttmv_data <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/ttmv_data.rds")
ttmv1 <- unique(ttmv_data$sample_id[grepl("RAR", ttmv_data$gene)])
ttmv2 <- unique(ttmv_data$sample_id[grepl("chr17:40", ttmv_data$Interval)])
ttmv_samples <- unique(c(ttmv1, ttmv2, "aplsz@s202014392_344698"))

# load counts data and map current gene symbol names to new names
all_leu_counts <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/counts8608.rds")
hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leukemia/data/hgnc_complete_set_2023-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id")]
colnames(hugo_anno)[3] <- "gene_id"
counts_old_anno <- read_csv("datExpBeforeSVA_counts.csv")
new_counts <- counts_old_anno[, -2] %>% left_join(hugo_anno) %>% na.omit() %>% 
  dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>% remove_rownames() %>% 
  column_to_rownames(var = "symbol") %>% as.matrix() %>% round()
counts_dat <- cbind(all_leu_counts, new_counts[rownames(all_leu_counts), ])[, colnames(dat)]

# change exp names to newest annotation
old_new <- left_join(counts_old_anno[, 1:2], hugo_anno)
old_new <- old_new %>% mutate(final_symbol = coalesce(symbol, gene_name))
gene_index <- match(rownames(dat), old_new$gene_name)
rownames(dat) <- old_new$final_symbol[gene_index]


# modify metadata for plot
meta_plot <- as.data.frame(metadata)[, c("gender", "age", "rna_type", "disease_type", "datasets")]
new_samples <- c(paste0("B00", 1:7), "B009", "B010")
meta_add <- data.frame(gender = NA_character_, age = NA_integer_, rna_type = "unknown",
                       disease_type = "AML", datasets = "aplsz_new", sample_id = new_samples)
meta_add <- meta_add %>% remove_rownames() %>% column_to_rownames(var = "sample_id")
meta_plot <- rbind(meta_plot, meta_add)

# 1st round annotation based on clusterings and datasets
# load clustering results
clusterings <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx")
clusterings <- as.data.frame(clusterings)
colnames(clusterings)[2] <- "rna_group"
g1g2samples <- clusterings$sample_id[clusterings$rna_group %in% c("G1", "G2")]
pmlsetsamples <- cld$sample_id[cld$datasets %in% c("aplsz", "GSE172057", "liuhongxing", "rarg")]
all_apl_samples <- unique(c(g1g2samples, pmlsetsamples, new_samples))
all_apl_samples <- all_apl_samples[!grepl("mecom", all_apl_samples)]

# add 1 layer annotation to plot - by primary clustering and dataset annotations
fusions_rar <- fusions_new[grepl("RAR", fusions_new$fusion), ]
pml_rara1 <- unique(fusions_new$sample_id[fusions_new$fusion %in% c("PML-RARA", "RARA-PML")])
pml_rara2 <- setdiff(unique(fusions_new$sample_id[fusions_new$gene3 %in% c("RARA", "RARB", "RARG")]), pml_rara1)
apl_like_samples <- rownames(meta_plot)[meta_plot$ident %in% "APL-like"]
fusions1 <- fusions_new[fusions_new$sample_id %in% apl_like_samples, ]
fusions2 <- fusions1[grepl("RAR", fusions1$fusion), ]
fusions1 <- fusions1[grepl("PML", fusions1$fusion), ]
rarg_samples <- unique(fusions2[grepl("RARG", fusions2$fusion), ]$sample_id)
rarb_samples <- unique(fusions2[grepl("RARB", fusions2$fusion), ]$sample_id)
rara_other1 <- unique(fusions2[grepl("RARA", fusions2$fusion), ]$sample_id)  ## RARA-AS1 also count here
rara_other_samples <- setdiff(rara_other1, c(rarg_samples, rarb_samples))
pml_other1 <- unique(fusions1[grepl("PML", fusions1$fusion), ]$sample_id)
pml_other_samples <- setdiff(pml_other1, c(rarg_samples, rarb_samples, rara_other_samples))
# set annotation 
meta_plot$anno <- meta_plot$ident
meta_plot$anno[rownames(meta_plot) %in% rarg_samples] <- "RARG"
meta_plot$anno[rownames(meta_plot) %in% rarb_samples] <- "RARB"
meta_plot$anno[rownames(meta_plot) %in% rara_other_samples] <- "RARA-others"
meta_plot$anno[rownames(meta_plot) %in% pml_other_samples] <- "pml-others"
meta_plot$anno[rownames(meta_plot) %in% ttmv_samples] <- "TTMV-RARA"


# add 2 layer annotation to plot - by clustering and callers
vadj <- matrixStats::rowVars(as.matrix(dat))
dat_mat <- dat %>% dplyr::arrange(desc(vadj))  %>% t() %>% scale() %>% t()
ht <- draw(ComplexHeatmap::Heatmap(dat_mat[1:1000, ], column_km = 26, use_raster = F,
                                   clustering_method_columns = 'ward.D2'))
cl <- map_int(column_order(ht), length)
ci <- NULL
for (i in 1:length(cl)) { ci <- c(ci, paste0(rep("G", cl[i]), i)) }
ordered_samples <- data.frame(sample_id = colnames(dat_mat)[unlist(column_order(ht))], 
                              subgroup = ci)
ods <- left_join(ordered_samples, data.frame(sample_id = rownames(meta_plot), anno = meta_plot$anno))
ods$ct <- ifelse(ods$anno %notin% "AML", 1, 0)
perc_ods <- ods %>% group_by(subgroup) %>% summarise(perc = sum(ct) / n())
ambiguous_apl_samples <- ods$sample_id[ods$anno %in% "AML" & ods$subgroup %in% perc_ods$subgroup[perc_ods$perc > 0.8]]
rara_other2 <- unique(fusions2[fusions2$gene5 %in% "RARA" | fusions2$gene3 %in% "RARA", ]$sample_id) 
rara_other_samples2 <- setdiff(rara_other2, c(rarg_samples, rarb_samples))
# three samples: "suzhou@suzhou_aml_A267" "HRA002693@HRR719193"    "scmc_m2b@KC"
# set annotation
meta_plot$age <- as.numeric(meta_plot$age)
meta_plot$call <- "AML"
meta_plot$call[rownames(meta_plot) %in% pml_rara1] <- "APL"
meta_plot$call[rownames(meta_plot) %in% ambiguous_apl_samples] <- "APL-like"
meta_plot$call[rownames(meta_plot) %in% rarg_samples] <- "RARG"
meta_plot$call[rownames(meta_plot) %in% ttmv_samples] <- "TTMV-RARA"
meta_plot$call[rownames(meta_plot) %in% rarb_samples] <- "RARB"
meta_plot$call[rownames(meta_plot) %in% rara_other_samples2] <- "RARA-others"
meta_plot$call[rownames(meta_plot) %in% pml_other_samples] <- "pml-others"


# add 3 layer annotation to plot - based on clinical annotation and callers
apll_samples <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/anno/apll_samples_clinical.rds")
apll_rar_samples <- unique(c(rarg_samples, rarb_samples, rara_other_samples2, ttmv_samples))
meta_plot$pheno <- "AML"
meta_plot$pheno[rownames(meta_plot) %in% pml_rara1] <- "APL"
meta_plot$pheno[rownames(meta_plot) %in% apll_samples] <- "APL-like"
meta_plot$pheno[rownames(meta_plot) %in% apll_rar_samples] <- "RARA-others"


# --- now we have plot matrix and annotations for 02 part




