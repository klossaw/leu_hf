pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


datasets <- dir_ls("/cluster/home/yliang_jh/projects/leukemia_linxiangjie/hamlet", type = "directory") |>
  path_file()
datasets <- datasets[!datasets %in% c("logs", "rnaseq", "tmp")]
itdlist <- list()
fusionlist <- list()
for (i in 1:length(datasets)){
  tryCatch({
  fus <- read_delim(glue::glue("/cluster/home/yliang_jh/projects/leukemia_linxiangjie/hamlet/{datasets[i]}/fusion/{datasets[i]}.fuma"))
  fusionlist[[datasets[i]]]  <- fus
  # current itd output of hamlet is not right, so we just distinguish EXIST or not
  itd1 <- read_csv(glue::glue("/cluster/home/yliang_jh/projects/leukemia_linxiangjie/hamlet/{datasets[i]}/itd/{datasets[i]}.flt3.csv"))
  itd2 <- read_csv(glue::glue("/cluster/home/yliang_jh/projects/leukemia_linxiangjie/hamlet/{datasets[i]}/itd/{datasets[i]}.kmt2a.csv"))
  itd <- c(nrow(itd1), nrow(itd2))
  itdlist[[datasets[i]]]  <- itd
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# setdiff(datasets, names(itdlist))    ## may need re-run?
# "amlzy137" "amlzy139" "amlzy202" "amlzy249" "amlzy315"


fusion_all <- rbindlist(fusionlist, use.names = F,idcol = "sample_id")
colnames(fusion_all) <- c("sample_id", "left", "right", "large200k", "fc", "sf")
fusion_all <- separate_rows(fusion_all, left, sep = ":")
fusion_all <- separate_rows(fusion_all, right, sep = ":")
fusion_all$fusion <- paste0(fusion_all$left, "-", fusion_all$right)
fusion_all <- fusion_all[, c("sample_id", "fusion", "large200k", "fc", "sf")]
readr::write_rds(fusion_all, "/cluster/home/yjliu_jh/projects/leu_j/data/fusion_all_hamlet.rds")

itd_all <- data.frame(matrix(ncol = 2))
for(i in 1:length(itdlist)){
  itd_all[i, ] <- t(matrix(itdlist[[i]]))
}
colnames(itd_all) <- c("FLT3", "KMT2A")
rownames(itd_all) <- names(itdlist)
readr::write_rds(itd_all, "/cluster/home/yjliu_jh/projects/leu_j/data/itd_all_hamlet.rds")


# subset fusions for leukemia drivers
fusion_all <- tidyr::separate(fusion_all, col = "fusion", sep = "-", into = c("left", "right"))
cols_6 <- as.character(readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 
                                          6, skip = 3, n_max = 1, col_names = FALSE))
leu_drivers <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 6, 
                                  skip = 5, col_names = cols_6)
drivers <- c(leu_drivers$Gene, "PML", "NPM1", "CBFB", "BCR")
fusion_fil <- fusion_all[fusion_all$left %in% drivers | fusion_all$right %in% drivers, ]
fusion_fil$fusion <- paste0(fusion_fil$left, "-", fusion_fil$right)
head(sort(table(fusion_fil$fusion), decreasing = T), 29)

# check if characteristics of subgroups are matched
sinfo2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_20220809.xlsx")
ffil <- left_join(fusion_fil, sinfo2[, 1:3])

# G1 characteristic CEBPA    G2 AML-ETO   G3 PML-RARA   G4 NPM1  G5 NPM1 KMT2A with FLT3-ITD
# G6 not obvious (NO)   G7 NPM1 KMT2A  similar to G5?   G8 CBFB-MYH11   G9 NO    G10 NO maybe BCR-ABL1

# head(sort(table(ffil[ffil$subgroups %in% "G2", ]$fusion), decreasing = T), 29)
# after checking, the fusion characteristics in sampleinfo can be well matched to the results called by hamlet
# mutation needs further filtering






