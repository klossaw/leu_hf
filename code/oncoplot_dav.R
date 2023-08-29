pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


hamlet_dav <- read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/hematology_dav/snv-indels_variants_all.csv")
hamlet_dav_fil <- hamlet_dav[is.na(hamlet_dav$filters), ]
#head(sort(table(unique(hamlet_dav_fil[, 1:2])$gene_symbol), decreasing = T), 20)

rna_edit <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/TABLE1_hg38.txt")
editdata <- rna_edit[, c(1:7)]
colnames(editdata)[1:3] <- c("CHROM", "POS", "REF")
temph <- inner_join(hamlet_dav_fil[, 1:10], editdata)
temph$flag <- mapply(grepl, temph$Ed, temph$alleles)
temph <- temph[temph$flag == TRUE, 1:10]
tempa <- anti_join(hamlet_dav_fil, temph)
tempa1 <- tempa[tempa$IMPACT %notin% c("MODIFIER"), ]
tempa2 <- tempa1[, c("sample_id", "gene_symbol", "VARIANT_CLASS")]
tempa2$VARIANT_CLASS <- ifelse(tempa2$VARIANT_CLASS == "SNV", "SNV", "INDEL")

indel_part <- tempa2[tempa2$VARIANT_CLASS == "INDEL", ]
indel_part$VARIANT_CLASS <- 1
indel_part <- unique(indel_part)
indel_m <- indel_part %>% pivot_wider(names_from = sample_id, values_from = VARIANT_CLASS, values_fill = 0) %>% 
  as.data.frame() %>% remove_rownames() %>% column_to_rownames(var = "gene_symbol") %>% as.matrix()

snv_part <- tempa2[tempa2$VARIANT_CLASS == "SNV", ]
snv_part$VARIANT_CLASS <- 1
snv_part <- unique(snv_part)
snv_m <- snv_part %>% pivot_wider(names_from = sample_id, values_from = VARIANT_CLASS, values_fill = 0) %>% 
  as.data.frame() %>% remove_rownames() %>% column_to_rownames(var = "gene_symbol") %>% as.matrix()

mat_list = list(snv = snv_m,
                indel = indel_m)
mat_list = unify_mat_list(mat_list)
oncoPrint(mat_list)





hamlet_leu <- read_csv("/cluster/home/yliang_jh/projects/leukemia_linxiangjie/hamlet/snv-indels_variants_all.csv", 
                       col_types = cols(CDS_position = col_character(), Protein_position = col_character()))



temph <- inner_join(hamlet_leu[, 1:10], editdata)
temph$flag <- mapply(grepl, temph$Ed, temph$alleles)
temph <- temph[temph$flag == TRUE, 1:10]
tempa <- anti_join(hamlet_leu, temph)
tempa1 <- tempa[tempa$IMPACT %notin% c("MODIFIER"), ]
tempa2 <- tempa1[, c("sample_id", "gene_symbol", "VARIANT_CLASS")]















