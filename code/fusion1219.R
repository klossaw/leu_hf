# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


base_path <- "/cluster/home/jhuang/projects/leukemia/analysis"



# fusion data for all datasets
datasets <- dir_ls("/cluster/home/jhuang/projects/leukemia/analysis", type = "directory") |>
  path_file()


# use variables from clustering1217.R
data_sets <- data.frame(sample_id = colnames(of_fil), dataset = temp1)




fc <- list()
for (i in 1:nrow(data_sets)){
  tryCatch({
    temp <- readr::read_delim(glue("{base_path}/{data_sets[i, 2]}/human/rnaseq/fusion/fusioncatcher/{data_sets[i, 1]}/final-list_candidate-fusion-genes.txt"))
    fc[[data_sets[i, 1]]] <- temp
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

readr::write_rds(fc, "/cluster/home/yjliu_jh/share/fusions_fusioncatcher.rds")

ar <- list()
for (i in 1:nrow(data_sets)){
  tryCatch({
    temp <- readr::read_delim(glue("{base_path}/{data_sets[i, 2]}/human/rnaseq/fusion/arriba/{data_sets[i, 1]}/fusions.tsv"))
    colnames(temp)[1] <- "gene1"
    ar[[data_sets[i, 1]]] <- temp
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

readr::write_rds(ar, "/cluster/home/yjliu_jh/share/fusions_arriba.rds")


sf <- list()
for (i in 1:nrow(data_sets)){
  tryCatch({
    temp <- readr::read_delim(glue("{base_path}/{data_sets[i, 2]}/human/rnaseq/fusion/starfusion/{data_sets[i, 1]}/star-fusion.fusion_predictions.tsv"))
    colnames(temp)[1] <- "FusionName"
    sf[[data_sets[i, 1]]] <- temp
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

readr::write_rds(sf, "/cluster/home/yjliu_jh/share/fusions_starfusion.rds")




