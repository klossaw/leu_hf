# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# check if sample_ids match
leu <- readr::read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")
sampleinfo_r <- colData(leu)
sampleinfo2 <- readr::read_csv("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo2.csv")

all(sampleinfo_r$sample_id %in% sampleinfo2$sample_id)  ## true

# reorder the rows and add new cols
sampleinfo3 <- sampleinfo2[match(sampleinfo_r$sample_id, sampleinfo2$sample_id), ]
# manage different annotations:
#sampleinfo3$rna_type_r <- sampleinfo_r$rna_type
#sampleinfo3$gender_r <- sampleinfo_r$gender
#sampleinfo3 <- sampleinfo3 %>% mutate(rna_type = coalesce(rna_type_r, rna_library),
#                                      gender = coalesce(gender_r, gender))


new_cols_newname <- c("age", "old_gender", "old_rna_type", "subgroup", "sample_index",
                      "outcome", "fusion_genes", "mutations", "disease_type")
new_cols <- c("age", "gender", "rna_library", "subgroup", "sample_index",
              "outcome", "gene_fusion", "mutations", "disease_type")

sampleinfo_r[, new_cols_newname] <- sampleinfo3[, new_cols]

# save data
colData(leu) <- sampleinfo_r
readr::write_rds(leu, "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")

