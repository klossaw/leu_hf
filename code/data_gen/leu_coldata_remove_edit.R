# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# read data
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
leu <- readr::read_rds(rds_fn)
meta_leu <- metadata(leu)
tsne_leu <- as.data.frame(meta_leu$tsne)
umap_leu <- as.data.frame(meta_leu$umap)
rownames(tsne_leu) <- rownames(umap_leu)
sampleinfo_r <- colData(leu)

sampleinfo_r2 <- sampleinfo_r[match(sampleinfo_r$sample_id, rownames(tsne_leu)), ]
tsne_leu$disease_type <- sampleinfo_r2$disease_type 
#tb_disease <- table(tsne_leu$disease_type) 
#r_diseases <- names(tb_disease)[tb_disease >20]

# edit clusters
# note: works only once! at 2023.1.2 4pm
tsne_leu[tsne_leu$tSNE1 < 30 & tsne_leu$tSNE1 > 12 & tsne_leu$tSNE2 < -30, ]$disease_type <- "cll"
remove_samples <- rownames(tsne_leu[tsne_leu$tSNE1 < -5 & tsne_leu$tSNE1 > -20 & tsne_leu$tSNE2 > 50, ])

remove_samples <- c("SJASPS030009_D1",  "SJASPS030015_D1",  "SJDSRCT030010_D1",
                    "SJDSRCT030016_D1", "SJEPD003_D",       "SJEWS001303_D1",  
                    "SJGIST030018_D1",  "SJHGG001_A",       "SJHGG003_A",       
                    "SJIFS030022_D1",   "SJLGG020_D",       "SJLGG022_D",      
                    "SJLGG026_D",       "SJLGG035_D",       "SJLGG036_D",       
                    "SJLGG038_D",       "SJLGG039_D",       "SJMB002_D",       
                    "SJMB009_E",        "SJMB028_D",        "SJMB030020_D1",    
                    "SJMB030021_D1",    "SJMB067_D",        "SJMEL001004_D2",  
                    "SJMPNST030013_D1", "SJNBL030002_D1",   "SJNBL030003_D1",   
                    "SJNBL030014_D1",   "SJNBL101_D",       "SJNBL105_D",      
                    "SJNBL107_D",       "SJNBL110_D",       "SJNBL117_D",       
                    "SJOS001_M",        "SJOS013_D",        "SJRB002_D",       
                    "SJRB031_D",        "SJRB051_D",        "SJRHB007_D",       
                    "SJRHB010_D",       "SJRHB013_D",       "SJSS030012_D1",   
                    "SJSS030017_D1",    "SRR5991033",       "SRR5991034",       
                    "SRR5991035",       "SRR5991036",       "SRR5991037",      
                    "SRR5991038",       "SRR5991039",       "SRR5991040",       
                    "scmc38")

colData(leu)$remove[sampleinfo_r$sample_id %in% remove_samples] <- "yes"

readr::write_rds(leu, rds_fn)



# write data for samples to be removed (2023.1.2)
remove_data <- data.frame(sample_id = remove_samples)
remove_data <- left_join(remove_data, as.data.frame(sampleinfo_r[, c("sample_id", "dataset")]))
colnames(remove_data)[2] <- "dataset_name"
remove_data$curator_id = "yjliu"
readr::write_csv(remove_data, "/cluster/home/yjliu_jh/share/remove_samples.csv")
