# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

toadd <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_add.rds")

rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
sampleinfo_r <- as.data.frame(colData(leu))
# should not be like this if properly ran quick_revise_sampleinfo.R ...
sampleinfo_r[, c("karyotype", "age_group")] <- NA
sampleinfo_r$age <- as.numeric(sampleinfo_r$age)
sampleinfo_r$karyotype <- as.character(sampleinfo_r$karyotype)
sampleinfo_r$age_group <- as.character(sampleinfo_r$age_group)
sampleinfo_r$subgroup <- as.character(sampleinfo_r$subgroup)
sampleinfo_r$mutations <- as.character(sampleinfo_r$mutations)
sampleinfo_r$outcome <- as.character(sampleinfo_r$outcome)

test2 <- dplyr::rows_patch(as.data.frame(sampleinfo_r), toadd, by = "sample_id", unmatched = "ignore")

for (i in 1:ncol(test2)) {
  colData(leu)[[colnames(test2)[i]]] <- test2[[colnames(test2)[i]]]
}


# check before overwrite!
readr::write_rds(leu, rds_fn)


rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
sampleinfo_r <- as.data.frame(colData(leu))
temp_sinfo <- sampleinfo_r
temp_sinfo$new_sample_id <- temp_sinfo$sample_id
ind <- temp_sinfo$dataset %in% "EGAS00001001952"
temp_sinfo$new_sample_id[ind] <- substr(temp_sinfo$new_sample_id[ind], 1, 12)

info_path <- "/cluster/home/yjliu_jh/projects/leu_j/data/leukemia_anno_new.xlsx"
info952 <- readxl::read_excel(info_path, sheet = "EGAS00001001952")
info952 <- info952[, c("sample_id", "Age_Dx_yrs")]
temp_join <- data.frame(sample_id = temp_sinfo$new_sample_id[ind])
temp_join <- left_join(temp_join, info952)
sampleinfo_r$age[ind] <- as.numeric(temp_join$Age_Dx_yrs)
colData(leu)$age <- sampleinfo_r$age
readr::write_rds(leu, rds_fn)


sinfo_beataml <- readxl::read_excel(info_path, sheet = "beataml")
ba_bld <- sinfo_beataml$sample_id[sinfo_beataml$sample_type %in% "Blood Derived Cancer - Peripheral Blood"]
ba_norm <- sinfo_beataml$sample_id[sinfo_beataml$sample_type %in% "Blood Derived Normal"]
readr::write_rds(ba_bld, "/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/rnaseq/samples/remove_samples_0104_ba1.rds")
readr::write_rds(ba_norm, "/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/rnaseq/samples/remove_samples_0104_ba2.rds")
readr::write_rds(c(ba_norm, ba_bld), "/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/rnaseq/samples/remove_samples_0104_ba_all.rds")




rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
remove_samples <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/rnaseq/samples/remove_samples_0104_ba_all.rds")
col_leu <- as.data.frame(colData(leu))
col_leu$remove[col_leu$sample_id %in% remove_samples] <- "yes"
colData(leu)$remove <- col_leu$remove
readr::write_rds(leu, rds_fn)






