pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", 
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


meta <- read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds") %>%
  colData() %>% as.data.frame()
samples <- dir_ls(path = "/cluster/home/jhuang/projects/leukemia/data/", glob = "*/human/rnaseq/*_R1.fq.gz",
                  type = "file", recurse = T) %>% .[!str_detect(.,"bak|dup")] %>% path_file() %>% str_remove("_R1.fq.gz")
datasets <- dir_ls(path = "/cluster/home/jhuang/projects/leukemia/data/", glob = "*/human/rnaseq/*_R1.fq.gz",
                   type = "file", recurse = T) %>% .[!str_detect(.,"bak|dup|upn")] %>% strsplit(split = "/") %>% 
  lapply(FUN = function(x){x[8]}) %>% unlist() %>% unname()
tab <- data.frame(sample_id = samples, datasets = datasets) %>% 
  dplyr::filter(!datasets == "bak")
metadata <- tab %>% left_join(meta, by = "sample_id") %>% dplyr::select("sample_id", "datasets", "age", "gender",
                                                                        starts_with("subgroups_"), starts_with("subtype_"), starts_with("PMID"), "fusion_genes", "mutations", "disease_type", 
                                                                        "rna_type") %>% mutate(sample_index = 1:nrow(tab)) ## %>% rename("rna_library" = "rna_type")

metadata %>% df2excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_leu_jh2.xlsx")

# already removed these data of other cancers
# temp <- metadata[metadata$dataset %in% "EGAS00001002217", ]
# other_samples <- temp[is.na(temp$gender), c("sample_id", "datasets")]
# other_samples %>% df2excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_other_diseases.xlsx")





datasets <- dir_ls(path = "/cluster/home/jhuang/projects/leukemia/data/", glob = "*/human/rnaseq/*_R1.fq.gz",
                   type = "file", recurse = T) %>% .[!str_detect(.,"bak|dup|PRJNA852777")] %>% strsplit(split = "/") %>% 
  lapply(FUN = function(x){x[8]}) %>% unlist() %>% unname() %>% unique()
i=1
mat <- vector("list", length(datasets))
merge_func <- function(x,y){left_join(x,y, by = c("gene_id", "gene_name"))}
for (dataset in datasets){
  csv <- glue("/cluster/home/jhuang/projects/leukemia/analysis/{dataset}/human/rnaseq/exp/tables/{dataset}_human.csv")
  exp <- read_csv(csv)
  mat[[i]]  <- exp
  i=i+1
}
names(mat) <- datasets
overall_exp <- mat %>% Reduce(f = merge_func) %>% dplyr::select(-starts_with("gene_name.")) %>%
  relocate(gene_id, gene_name, everything())
overall_exp %>% write_csv("/cluster/home/yjliu_jh/projects/temp/overall_exp515.csv")

# count matrix
i=1
mat <- vector("list", length(datasets))
merge_func <- function(x,y){left_join(x,y, by = c("gene_id", "gene_name"))}
for (dataset in datasets){
  csv <- glue("/cluster/home/jhuang/projects/leukemia/analysis/{dataset}/human/rnaseq/exp/tables/{dataset}_human_counts.csv")
  exp <- read_csv(csv)
  mat[[i]]  <- exp
  i=i+1
}
names(mat) <- datasets
overall_exp <- mat %>% Reduce(f = merge_func) %>% dplyr::select(-starts_with("gene_name.")) %>%
  relocate(gene_id, gene_name, everything()) 
overall_exp %>% write_csv("/cluster/home/yjliu_jh/projects/temp/overall_exp_counts515.csv")




sinfo <- readxl::read_excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_leu_jh2.xlsx")

# get these collected annotations
#
anno1 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/GSE172057_anno.xlsx", sheet = 22, skip = 1)
colnames(anno1)[1] <- "patient_id"
anno1 <- anno1[-((nrow(anno1) - 2) : nrow(anno1)), ]
anno1$patient_id <- gsub("APL(\\d{1}AC)", "APL0\\1", anno1$patient_id)  ## change sample names
anno1$patient_id <- gsub("APL(\\d{2}AC)", "APL0\\1", anno1$patient_id) 
anno1$patient_id <- sub("AC", "", anno1$patient_id)
anno1p <- anno1[, c(1, 8:15)] %>%  pivot_longer(-patient_id) %>%
  group_by(patient_id) %>% summarise(mutation = paste(name[value > 0L], collapse = '; '))

ids <- intersect(anno1$patient_id, sinfo$sample_id)
ind1 <- which(sinfo$sample_id %in% ids)
ind1x <- match(ids, anno1$patient_id)

sinfo[ind1, "mutations"] <- anno1p$mutation[ind1x]
sinfo[ind1, "gender"] <- tolower(anno1$Gender[ind1x])
sinfo[ind1, "age"] <- anno1$Age[ind1x]



# HRA00113
anno2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA001135_anno.xlsx")


# 
anno3 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA002693_anno.xlsx")




sinfo %>% jhtools::df2excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_leu_jh4.xlsx")


