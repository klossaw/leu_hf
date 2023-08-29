# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# keep <subgroups> column newest
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
new_labels <- read_csv("/cluster/home/ylxie_jh/projects/leukemia/analysis/meta/jhuang/human/figures/heatmap/need_to_add_info_230823.csv")
colData(rds)$subgroups_230712 <- new_labels$subgroups_230712  ## orders are identical
colData(rds)$subgroups <- colData(rds)$subgroups_230712
readr::write_rds(rds, rds_fn)

