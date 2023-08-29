pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2", "ggthemes", "EnsDb.Hsapiens.v86",
          "jhuanglabRNAseq", "clusterProfiler", "DOSE", "org.Hs.eg.db", "rstatix", "maftools")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

new_samples <- c(paste0("B00", 1:7), "B009", "B010")

#  ======= mutations =======
# read in new mutation data
sel_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
              "Variant_Classification", "Variant_Type", "Consequence", "PUBMED", "FILTER",
              "HGVSc", "HGVSp_Short", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
              "SIFT", "PolyPhen", "IMPACT", "CLIN_SIG", "Existing_variation", "Protein_position")
maf_loc <- "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/mutations/apll/maf"

# read
speedseq <- list()
for(i in 1:length(new_samples)){
  tryCatch({
    temp <- read.maf(glue("{maf_loc}/{new_samples[i]}/{new_samples[i]}.speedseq.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    speedseq[[new_samples[i]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

hc <- list()
for(i in 1:length(new_samples)){
  tryCatch({
    temp <- read.maf(glue("{maf_loc}/{new_samples[i]}/{new_samples[i]}.HaplotypeCaller.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    hc[[new_samples[i]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ri <- list()
for(i in 1:length(new_samples)){
  tryCatch({
    temp <- read.maf(glue("{maf_loc}/{new_samples[i]}/{new_samples[i]}.rnaindel.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    ri[[new_samples[i]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

vs2s <- list()
for(i in 1:length(new_samples)){
  tryCatch({
    temp <- read.maf(glue("{maf_loc}/{new_samples[i]}/{new_samples[i]}.snp.varscan2.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    vs2s[[new_samples[i]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

vs2i <- list()
for(i in 1:length(new_samples)){
  tryCatch({
    temp <- read.maf(glue("{maf_loc}/{new_samples[i]}/{new_samples[i]}.indel.varscan2.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    vs2i[[new_samples[i]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

mut_hc <- dplyr::bind_rows(hc, .id = "sample_id")
mut_ss <- dplyr::bind_rows(speedseq, .id = "sample_id")
mut_ri <- dplyr::bind_rows(ri, .id = "sample_id")
mut_vs <- dplyr::bind_rows(vs2s, .id = "sample_id")
indel_vs <- dplyr::bind_rows(vs2i, .id = "sample_id")
mut_hc$tool <- "HC"
mut_ss$tool <- "SS"
mut_vs$tool <- "VS"
mut_ri$tool <- "RI"
indel_vs$tool <- "VSI"
mut_int <- bind_rows(mut_hc, mut_ss, mut_ri, mut_vs, indel_vs)
readr::write_rds(mut_int, "/cluster/home/yjliu_jh/projects/leu_j/data/apll_new_mutations.rds")


#  ======= fusions =======
get_comb <- function(path, name, samples, tool) {
  temp_list <- list()
  for (i in 1:length(samples)) {
    temp_list[[samples[i]]] <- read_tsv(glue("{path}/{samples[i]}/{name}"))
  }
  temp <- bind_rows(temp_list, .id = "sample_id")
  temp$tool <- tool
  temp
}

fu_ar <- get_comb("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/fusion/arriba", "fusions.tsv", new_samples, "AR")
fu_fc <- get_comb("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/fusion/fusioncatcher", "final-list_candidate-fusion-genes.txt", new_samples, "FC")
fu_sf <- get_comb("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/fusion/starfusion", "star-fusion.fusion_predictions.tsv", new_samples, "SF")

fu_ar_fil <- unique(fu_ar[, c(1:3, ncol(fu_ar))])
colnames(fu_ar_fil) <- c("sample_id", "gene5", "gene3", "tool")
fu_ar_fil$fusion <- paste0(fu_ar_fil$gene5, "-", fu_ar_fil$gene3)
fu_fc_fil <- unique(fu_fc[, c(1:3, ncol(fu_fc))])
colnames(fu_fc_fil) <- c("sample_id", "gene5", "gene3", "tool")
fu_fc_fil$fusion <- paste0(fu_fc_fil$gene5, "-", fu_fc_fil$gene3)
fu_sf_fil <- unique(fu_sf[, c(1:2, ncol(fu_sf))])
fu_sf_fil <- fu_sf_fil %>% separate(`#FusionName`, into = c("gene5", "gene3"), sep = "--")
fu_sf_fil$fusion <- paste0(fu_sf_fil$gene5, "-", fu_sf_fil$gene3)
fu_all <- bind_rows(fu_ar_fil, fu_fc_fil, fu_sf_fil)
fu_all <- cbind(data.frame(dataset = "apll_new"), as.data.frame(fu_all))

readr::write_rds(fu_all, "/cluster/home/yjliu_jh/projects/leu_j/data/apll_new_fusions.rds")






