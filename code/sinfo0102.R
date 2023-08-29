# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



# keep to newest update of sampleinfo
leu <- readr::read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")
sampleinfo_r <- as.data.frame(colData(leu))


info_path <- "/cluster/home/yjliu_jh/projects/leu_j/data/leukemia_anno_new.xlsx"
collected_info_sheets <- readxl::excel_sheets(info_path)
collected_info <- lapply(collected_info_sheets, function(X) readxl::read_excel(info_path, sheet = X))
names(collected_info) <- collected_info_sheets


coln_info <- list()
dtst <- character()
for (i in 1:length(collected_info)) {
  coln_info[[i]] <- colnames(collected_info[[i]])
  dtst <- c(dtst, rep(collected_info_sheets[i], ncol(collected_info[[i]])))
}

col_match <- data.frame(dataset = dtst, col = unlist(coln_info))
colns <- unique(col_match$col)


# needed column: age, subgroup, index? fusion, mutation, outcome, age_group, karyotype
# already had:  disease_type, gender, rna_library


# ============ age =================
# collect age columns and change days to years
age_related_colnames <- colns[grep("age", colns, ignore.case = T)]
age_day <- col_match$dataset[col_match$col %in% "age(days)"]
for (i in 1:length(age_day)){
  collected_info[[age_day[i]]]$"age(days)" <- as.numeric(collected_info[[age_day[i]]]$"age(days)") / 365.25
}

age_diag <- col_match$dataset[col_match$col %in% "age_at_diagnosis"]
for (i in 1:length(age_diag)){
  collected_info[[age_diag[i]]]$"age_at_diagnosis" <- as.numeric(collected_info[[age_diag[i]]]$"age_at_diagnosis") / 365.25
}

# unify colnames to "age"
arc <- age_related_colnames[c(1, 4, 6, 8, 12, 13, 14, 15, 16)]

for(i in 1:length(arc)){
  age_datasets <- col_match$dataset[col_match$col %in% arc[i]]
   for(j in 1:length(age_datasets)){
     colnames(collected_info[[age_datasets[j]]])[colnames(collected_info[[age_datasets[j]]]) %in% arc[i]] <- "age"
     collected_info[[age_datasets[j]]]$age <- as.numeric(collected_info[[age_datasets[j]]]$age)
   }
}




# ========= sub-group-related columns ===========
subgroup_related_colnames1 <- colns[grep("sub", colns, ignore.case = T)]

# unify colnames to "group"
sgrc <- c(subgroup_related_colnames1[c(2, 4, 5, 7, 9, 10, 12, 13)], "Group")

for(i in 1:length(sgrc)){
  subgroup_datasets <- col_match$dataset[col_match$col %in% sgrc[i]]
  for(j in 1:length(subgroup_datasets)){
    colnames(collected_info[[subgroup_datasets[j]]])[colnames(collected_info[[subgroup_datasets[j]]]) %in% sgrc[i]] <- "subgroup"
  }
}



# ======== fusion (annotated by original sources)  =========
fusion_related_colnames <- colns[grep("fusion", colns, ignore.case = T)]
# revise one column
tmpdset <- col_match$dataset[grep("RT-PCR for", col_match$col)]
collected_info[[tmpdset]]$fusion <- ifelse(collected_info[[tmpdset]]$"RT-PCR for \nSPI1 fusion" == "Yes", "SPI1 fusion", NA)

# unify colnames to "fusion"
frc <- fusion_related_colnames[c(1, 4, 6)]

for(i in 1:length(frc)){
  fusion_datasets <- col_match$dataset[col_match$col %in% frc[i]]
  for(j in 1:length(fusion_datasets)){
    colnames(collected_info[[fusion_datasets[j]]])[colnames(collected_info[[fusion_datasets[j]]]) %in% frc[i]] <- "fusions"
  }
}


# ======== mutation (annotated by original sources)  =========
mutation_related_colnames <- colns[grep("mutat", colns, ignore.case = T)]


# unify colnames to "mutation"
mrc <- c("mutation", "oncogenes")

for(i in 1:length(mrc)){
  mutation_datasets <- col_match$dataset[col_match$col %in% mrc[i]]
  for(j in 1:length(mutation_datasets)){
    colnames(collected_info[[mutation_datasets[j]]])[colnames(collected_info[[mutation_datasets[j]]]) %in% mrc[i]] <- "mutations"
  }
}


# ======== outcome ===========================
outcome_related_colnames <- colns[grep("outcome", colns, ignore.case = T)]

# unify colnames to "outcome"
orc <- "Outcome"

for(i in 1:length(orc)){
  outcome_datasets <- col_match$dataset[col_match$col %in% orc[i]]
  for(j in 1:length(outcome_datasets)){
    colnames(collected_info[[outcome_datasets[j]]])[colnames(collected_info[[outcome_datasets[j]]]) %in% orc[i]] <- "outcome"
  }
}

# ======== age_group =================
group_related_colnames <- colns[grep("group", colns, ignore.case = T)]

agrc <- "age Group*"

for(i in 1:length(agrc)){
  agegroup_datasets <- col_match$dataset[col_match$col %in% agrc[i]]
  for(j in 1:length(agegroup_datasets)){
    colnames(collected_info[[agegroup_datasets[j]]])[colnames(collected_info[[agegroup_datasets[j]]]) %in% agrc[i]] <- "age_group"
  }
}

# ======== karyotype =================
karyotype_related_colnames <- colns[grep("karyotype", colns, ignore.case = T)]

krc <- "Karyotype"

for(i in 1:length(krc)){
  karyotype_datasets <- col_match$dataset[col_match$col %in% krc[i]]
  for(j in 1:length(karyotype_datasets)){
    colnames(collected_info[[karyotype_datasets[j]]])[colnames(collected_info[[karyotype_datasets[j]]]) %in% krc[i]] <- "karyotype"
  }
}


# combine data

readr::write_rds(collected_info, "/cluster/home/yjliu_jh/projects/leu_j/data/collected_info_0103.rds")
new_anno_columns <- c("age", "subgroup", "fusions", "mutations", "outcome", "age_group", "karyotype")

cinfo_sub <- list()
for (i in 1:length(collected_info)) {
  cinfo_sub[[names(collected_info)[i]]] <- collected_info[[i]][colnames(collected_info[[i]]) %in% c("sample_id", "dataset", new_anno_columns)]
}

readr::write_rds(cinfo_sub, "/cluster/home/yjliu_jh/projects/leu_j/data/cinfo_sub_0103.rds")

test1 <- bind_rows(cinfo_sub)
colnames(test1)[5] <- "fusion_genes"
test1 <- unique(test1)
testxx <- as.data.frame(janitor::get_dupes(test1, sample_id))
testxx1 <- testxx[testxx$dataset %in% c("stjude", "EGAS00001001952"), ]
testxx2 <- testxx[testxx$dataset %notin% c("stjude", "EGAS00001001952"), ]
testxx3 <- testxx[testxx$sample_id %in% setdiff(unique(testxx$sample_id), unique(testxx2$sample_id)), ]
testxx3 <- testxx3[testxx3$dataset %notin% "stjude", ]

testxx2_1 <- testxx2[testxx2$sample_id %in% names(table(testxx2$sample_id)[table(testxx2$sample_id) == 1]), ]
testxx2_2 <- testxx2[testxx2$sample_id %in% names(table(testxx2$sample_id)[table(testxx2$sample_id) >= 2]), ]
keep2 <- testxx2_2[!duplicated(testxx2_2$sample_id), ] 
keep2 <- rbind(testxx2_1, keep2)
keep <- rbind(testxx3, keep2)

toadd <- anti_join(test1, testxx)
toadd <- rbind(toadd, keep[, -2])

readr::write_rds(toadd, "/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_add.rds")


sampleinfo_r[, c("karyotype", "age_group")] <- NA
sampleinfo_r$age <- as.numeric(sampleinfo_r$age)
sampleinfo_r$karyotype <- as.character(sampleinfo_r$karyotype)
sampleinfo_r$age_group <- as.character(sampleinfo_r$age_group)
sampleinfo_r$subgroup <- as.character(sampleinfo_r$subgroup)
sampleinfo_r$mutations <- as.character(sampleinfo_r$mutations)
sampleinfo_r$outcome <- as.character(sampleinfo_r$outcome)


test2 <- dplyr::rows_patch(as.data.frame(sampleinfo_r), toadd, by = "sample_id", unmatched = "ignore")

  



tempdt <- setdiff(unique(sinfo3$dataset), collected_info_sheets)








































# 1 make sure dataset and sample_id exists in sinfo

# HRA000122 match sample      
# HRA000789 already have anno, ignore
# HRA000489 match anno from S5
# phs000218 output from target_aml and 218 extra and add to excel
# phs000873 
# pnas_tall 



# ====== deleted =====
# SRR1909157 should be deleted, as it's bone marrow


# needed column: age, subgroup, index? fusion, mutation, outcome, age_group, karyotype
# already had:  disease_type, gender, rna_library



# 1 match sample_id 
# phs000218 already removed

# HRA000122

HRA000122_anno <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000122_anno.xlsx", sheet = 1)
HRA000122 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000122.xlsx", sheet = 1)
HRA000122_2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000122.xlsx", sheet = 2)
colnames(HRA000122)[4] <- "Sample"
colnames(HRA000122_2)[1] <- "sample_id"

hra122 <- left_join(HRA000122_2, HRA000122)
colnames(hra122)[7] <- "PatientID"
hra122 <- left_join(hra122, HRA000122_anno, by = "PatientID")

readr::write_csv(hra122, "/cluster/home/yjliu_jh/projects/leu_j/data/hra122_sinfo.csv")


# HRA000489

hra489_s5 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000489_S5.xlsx", skip = 1, sheet = 5)
hra489 <- collected_info[["HRA000489"]]
hra489_s5$patient_id <- paste0("BALL_", hra489_s5$`Patient No.`)
hra489 <- left_join(hra489, hra489_s5)
readr::write_csv(hra489, "/cluster/home/yjliu_jh/projects/leu_j/data/hra489_sinfo.csv")


# HRA000789
HRA000789_anno <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000789_ANNO.xlsx", skip= 1)
HRA000789 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000789.xlsx", sheet = 1)
HRA000789_2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/HRA000789.xlsx", sheet = 2)

colnames(HRA000789_2)[1] <- "sample_id"
colnames(HRA000789)[4] <- "Sample"
hra789 <- left_join(HRA000789_2, HRA000789)
colnames(hra789)[7] <- colnames(HRA000789_anno)[1]
hra789 <- left_join(hra789, HRA000789_anno, by = colnames(HRA000789_anno)[1])
readr::write_csv(hra789, "/cluster/home/yjliu_jh/projects/leu_j/data/hra789_sinfo.csv")


