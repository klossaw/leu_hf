
sample_RNA_type <- function(dat){
  stopifnot("please make sure gene_name is contained in the data" = "gene_name" %in% colnames(dat))
  dat_mat <- dat %>% dplyr::select(-gene_id) %>% as.data.frame() %>% na.omit() %>%
    remove_rownames() %>% column_to_rownames(var = "gene_name")
  his_genes1 <- c("H4C11",  "H2BC3",  "H4C4",   "H3C12",  "H2BC14", "H4C6",   "H2AC16", "H2AC11",
                  "H2AC12", "H1−5",   "H3C2",   "H3C3",   "H3C11",  "H2AC17", "H2AC13", "H2BC13",
                  "H2AC4",  "H2AC14", "H4C9",   "H4C1",   "H4C2",   "H4C13",  "H3C7",   "H2BC9",
                  "H2BC10", "H3C8",   "H3C13",  "H2BC7",  "H2BC18", "H2BC15", "H4C14",  "H4C8",
                  "H2BC17", "H3C10",  "H2BC11", "H4C12",  "H3C15",  "H2BC8",  "H2AC8",  "H2BC6",
                  "H4C5",   "H2AC15", "H1−3",   "H3C4",   "H4−16",  "H2AC20", "H2AC21", "H3C1",
                  "H2BC5",  "H1−4",   "H2AC7",  "H2BC4",  "H1−2",   "H2AC6")
  inter_genes <- intersect(his_genes1, rownames(dat_mat))
  stopifnot("please make sure histone genes is contained in the data" = length(inter_genes) > 0)
  dat_mat <- dat_mat[rownames(dat_mat) %in% inter_genes, ]
  clust <- data.frame(cluster = cutree(hclust(dist(t(dat_mat))), k = 2))
  diff <- mean(unlist(dat_mat[, clust$cluster %in% 1])) - mean(unlist(dat_mat[, clust$cluster %in% 2]))
  clust$rna_type <- "total_RNA"
  clust$rna_type[clust$cluster == (2 - sum(diff > 0))] <- "mRNA"
  new_dat <- data.frame(sample_id = rownames(clust), rna_type = clust$rna_type)
  return(new_dat)
}





sample_gender <- function(dat){
  stopifnot("please make sure gene_name is contained in the data" = "gene_name" %in% colnames(dat))
  dat_mat <- dat %>% dplyr::select(-gene_id) %>% as.data.frame() %>% na.omit() %>%
    remove_rownames() %>% column_to_rownames(var = "gene_name")
  gender_genes <- c("KDM5D", "UTY", "USP9Y", "DDX3Y", "RPS4Y1", "ZFY", "EIF1AY", "ENSG00000229807", "XIST")
  inter_genes <- intersect(gender_genes, rownames(dat_mat))
  stopifnot("please make sure gender-linked genes is contained in the data" = length(inter_genes) > 0)
  dat_mat <- dat_mat[rownames(dat_mat) %in% inter_genes, ]
  clust <- data.frame(cluster = cutree(hclust(dist(t(dat_mat))), k = 2))
  if (inter_genes[1] %notin% c("ENSG00000229807", "XIST")){
    diff <- mean(unlist(dat_mat[inter_genes[1], clust$cluster %in% 1])) - mean(unlist(dat_mat[inter_genes[1], clust$cluster %in% 2]))
  } else {
    diff <- mean(unlist(dat_mat[inter_genes[1], clust$cluster %in% 2])) - mean(unlist(dat_mat[inter_genes[1], clust$cluster %in% 1]))
  }
  clust$gender <- "female"
  clust$gender[clust$cluster == (2 - sum(diff > 0))] <- "male"
  new_dat <- data.frame(sample_id = rownames(clust), gender = clust$gender)
  return(new_dat)
}



















