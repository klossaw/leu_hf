# cr: flora
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


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
overall_exp <- mat[-which(data_n == 0)] %>% Reduce(f = merge_func) %>% dplyr::select(-starts_with("gene_name.")) %>%
  relocate(gene_id, gene_name, everything())
overall_exp %>% write_csv("/cluster/home/yjliu_jh/projects/leu_j/data/overall_exp_count_short.csv")



data_n <- numeric()
for (i in 1:length(mat)){
  data_n[i] <- sum(ncol(mat[[i]]))
}





