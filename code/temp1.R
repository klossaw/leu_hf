pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



leuloc <- "/cluster/home/jhuang/projects/leukemia/analysis/linxiangjie/human/rnaseq"
fusions <- readxl::read_excel(glue::glue("{leuloc}/fusion/tables/linxiangjie_fusions.xlsx"))
fusions2 <- readxl::read_excel(glue::glue("{leuloc}/fusion/tables/linxiangjie_fusions.xlsx"), 2)
count_mat <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human.csv"))
count1 <- read_csv("/cluster/home/jhuang/projects/leukemia/analysis/linxiangjie/human/rnaseq/exp/tables/linxiangjie_human_counts.csv")
count1 <- count1 %>% as.data.frame() %>% column_to_rownames(var = "gene_id") %>% as.matrix()
sinfo1 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/SampleInfo_0.9.xlsx")
sinfo2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_20220809.xlsx")
heatgenes <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/SampleInfo_0.9.xlsx", 2)
heatgenes <- heatgenes$gene_name


ct2 <- sinfo2[, c("sample_id", "subgroups", "characteristic", "age", "gender", "mutation", "fusion", "fusion_type")]
count_mat <- count_mat[, -1] %>% na.omit() %>% as.data.frame() %>% column_to_rownames(var = "gene_name")
count_filtered <- count_mat[rownames(count_mat) %in% heatgenes, ]
cf <- count_filtered
ct1 <- cbind(rownames(count_filtered), count_filtered)
colnames(ct1)[1] <- "gene_name"

count_extended <- left_join(ct1, ct2)




cols_6 <- as.character(readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 
                                           6, skip = 3, n_max = 1, col_names = FALSE))
leu_drivers <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 6, 
                                  skip = 5, col_names = cols_6)

gene_path <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 14, skip = 3)







ct2e <- separate_rows(ct2, mutation, sep = ",")
#table(ct2e[ct2e$subgroups %in% "G2", "mutation"])
samples <- colnames(cf)
clusters <- unique(ct2e$subgroups)
ct2e <- mutate_all(ct2e, .funs = toupper)
gperc <- list()
for (i in 1:length(clusters)){
  gperc[[clusters[i]]] <- data.frame(table(ct2e[ct2e$subgroups %in% clusters[i], "characteristic"]))
  colnames(gperc[[clusters[i]]])[1] <- "Characteristic"
  gperc[[clusters[i]]]$Percentage <- gperc[[i]]$Freq / sum(gperc[[i]]$Freq)
}

perc <- dplyr::bind_rows(gperc, .id = 'Group')
perc[perc$Percentage > 0.1, ]

ct2e$characteristic[is.na(ct2e$characteristic)] <- "unknown"
curr_drivers <- unique(c(ct2e$characteristic, ct2e$mutation, ct2e$fusion))

setdiff(curr_drivers, leu_drivers$Gene)



mafloc <- "/cluster/home/jhuang/projects/leukemia/analysis/linxiangjie/human/rnaseq/mutations/linxiangjie/maf"
temp <- read.maf(paste0(mafloc, "/amlzy151/amlzy151.speedseq.maf"))
slotNames(temp)

temp2 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                       "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                       "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
colnames(temp2)[14] <- "VEP_Impact"
#as.numeric(object.size(temp2)) / as.numeric(object.size(temp))

mut_ss$SIFT <- gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", mut_ss$SIFT, perl = T)
mut_ss$PolyPhen <- gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", mut_ss$PolyPhen, perl = T)
mut_ss$SIFT <- as.numeric(mut_ss$SIFT)
mut_ss$PolyPhen <- as.numeric(mut_ss$PolyPhen)

mut_ss


group_by(across(c(-tool)))


# vcf_qual: Phred Quality Score  bad if too low   threshold 20: 1% no variant; 50: 1 in 1e5 chance
# gnomAD_AF: allele freq (population) that may not be accurate in some specific occasions?
# Variant allele freq (reads) = t_alt_count / (t_ref_count + t_alt_count)?




temp <- read.maf(paste0(mafloc, "/amlzy151/amlzy151.HaplotypeCaller.maf"))
temp <- read.maf(paste0(mafloc, "/amlzy151/amlzy151.DeepVariant.maf"))
dim(temp@data)
temp3 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                       "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                       "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
colnames(temp3)[14] <- "VEP_Impact"

temp <- read.maf(paste0(mafloc, "/amlzy151/amlzy151.DeepVariant.maf"))
dim(temp@data)
temp3 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                       "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                       "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
colnames(temp3)[14] <- "VEP_Impact"


temp1 <- read.maf(paste0(mafloc, "/amlzy165/amlzy165.speedseq.maf"))




speedseq <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/speedseq.rds")
mut_ss <- dplyr::bind_rows(speedseq, .id = 'Samples')
mut_ss <- mut_ss %>% mutate(vaf_sample = t_alt_count / (t_ref_count + t_alt_count))

# some filtering process...

ss1 <- mut_ss %>% group_by(Samples, Hugo_Symbol) %>% summarise(count = n())
ss1$ct <- 1
ss2 <- ss1 %>% group_by(Hugo_Symbol) %>% summarise(times = sum(ct))
ss2 <- ss2[order(ss2$times, decreasing = T), ]
ss2 <- ss2[ss2$Hugo_Symbol %notin% unique(flags$FLAGS), ]

tgenes <- ss2$Hugo_Symbol[1:20]

#[1] "ACIN1"      "ADAR"       "AIP"        "ASH1L"      "BAG1"       "BIVM-ERCC5" "CAP1"       "CCDC66"     "CD44"      
#[10] "CEBPZ"      "CGAS"       "COG4"       "CRELD1"     "CTSC"       "CYBA"       "CYC1"       "DDX60L"     "DENND3"    
#[19] "DIDO1"      "DNHD1"      "EHBP1L1"    "EIF2AK3"    "EIF4G1"     "FAAP100"    "FAM193B"    "FAM91A1"    "FANCA"     
#[28] "FMNL1"      "FXYD5"      "GART"    



mut_ss %>% 



head(match(leu_drivers$Gene, ss2$Hugo_Symbol), 180)


ss3 <- mut_ss %>% group_by(Samples, Hugo_Symbol, HGVSc) %>% summarise(count = n())
ss3$count <- 1  ## multiple transcripts of same gene? 
ss3 <- ss3 %>% group_by(Hugo_Symbol, HGVSc) %>% summarise(tcount = sum(count)) %>% 
  group_by(Hugo_Symbol) %>% mutate(allcount = sum(tcount), freq = tcount/allcount) 

ss4 <- ss3[ss3$freq < 0.8, ]
ss4 <- ss4[ss4$tcount < 100, ]
ss4 <- ss4 %>% group_by(Hugo_Symbol) %>% summarise(newcount = sum(tcount))
ss4 <- ss4[ss4$Hugo_Symbol %notin% unique(flags$FLAGS), ]
  
ssx <- ss4[order(ss4$newcount, decreasing = T), ]

head(match(leu_drivers$Gene, ssx$Hugo_Symbol), 180)


ssin <- mut_ss[mut_ss$Hugo_Symbol %in% c(leu_drivers$Gene, curr_drivers), ]





mut_ss[mut_ss$Hugo_Symbol %in% tgenes[3], ]
mut_ss[mut_ss$Hugo_Symbol %in% "DNMT3A", ]
mut_ss[mut_ss$Hugo_Symbol %in% "NPM1", ]


mut_ss[mut_ss$Hugo_Symbol %in% "AIP", ]
table(mut_ss[mut_ss$Hugo_Symbol %in% "BAG1", "HGVSc"])


table(mut_ss[mut_ss$Hugo_Symbol %in% "AIP", "HGVSc"])  ## SNP!   c.682C>A  in all 298 patients

table(mut_ss[mut_ss$Hugo_Symbol %in% "DNMT3A", "HGVSc"])  ## very reasonable 

table(mut_ss[mut_ss$Hugo_Symbol %in% "ATM", "HGVSc"])
table(mut_ss[mut_ss$Hugo_Symbol %in% "BCL6", "HGVSc"])



# SIFT seems better in passenger genes, because widespread SNPs are unlikely to be very dangerous?

# so try to filter like, data-based freq > 80%?   or varied loci > 3?




dv <- list()
for(i in 1:length(samples)){
  
  temp <- read.maf(glue("{mafloc}/{samples[i]}/{samples[i]}.DeepVariant.maf"))
  temp_3 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                         "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                         "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
  colnames(temp_3)[14] <- "VEP_Impact"
  dv[[samples[i]]] <- temp_3
}
readr::write_rds(dv, "/cluster/home/yjliu_jh/projects/leu_j/data/DeepVariant.rds")

dv <- list()
for (i in 1:length(samples)){
  tryCatch({
    temp <- read.maf(glue("{mafloc}/{samples[i]}/{samples[i]}.DeepVariant.maf"))
    temp_3 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                            "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                            "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
    colnames(temp_3)[14] <- "VEP_Impact"
    dv[[samples[i]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




speedseq <- list()
for(i in 1:length(samples)){
  temp <- read.maf(glue("{mafloc}/{samples[i]}/{samples[i]}.speedseq.maf"))
  temp_3 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                          "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                          "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
  colnames(temp_3)[14] <- "VEP_Impact"
  speedseq[[samples[i]]] <- temp_3
}
readr::write_rds(speedseq, "/cluster/home/yjliu_jh/projects/leu_j/data/speedseq.rds")


hc <- list()
for(i in 1:length(samples)){
  temp <- read.maf(glue("{mafloc}/{samples[i]}/{samples[i]}.HaplotypeCalle.maf"))
  temp_3 <- temp@data[, c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele",
                          "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
                          "SIFT", "PolyPhen", "IMPACT", "vcf_qual", "gnomAD_AF", "Protein_position")]
  colnames(temp_3)[14] <- "VEP_Impact"
  hc[[samples[i]]] <- temp_3
}
readr::write_rds(hc, "/cluster/home/yjliu_jh/projects/leu_j/data/HaplotypeCaller.rds")











combined <- 






  
  
  

  



flags <- read.delim("/cluster/home/yjliu_jh/projects/leu_j/data/FLAGS.txt")








ph <- Heatmap(cf,
              name = "expression",
              #col = GenerateHeatmapColor("bluered",5),
              col = circlize::colorRamp2(c(-1.8, 0, 1.8), c("blue", "white", "red")),
              cluster_columns = T,
              cluster_rows = rev(hm_s$rowDendrogram),
              show_row_names = T,
              show_column_names = F,
              top_annotation = ha_s,
              column_title = NULL
)





pdf(file = '/cluster/home/yjliu_jh/fnttt.pdf', width = 13, height = 9)
ph
dev.off()

