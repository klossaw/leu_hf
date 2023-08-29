# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# keep newest update of sampleinfo
leu <- readr::read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")
sampleinfo_r <- as.data.frame(colData(leu))
meta_leu <- metadata(leu)
rdim1 <- meta_leu[[1]]
sampleinfo_r2 <- sampleinfo_r[match(sampleinfo_r$sample_id, rownames(rdim1)), ]
dtype_col <- brewer.pal(n = 12, name = "Paired")
names(dtype_col) <- unique(sampleinfo_r2$disease_type)
#dtype_col <- dtype_col[!is.na(names(dtype_col))]
#dtype2 <- left_join(sampleinfo_r2, data.frame(disease_type = names(dtype_col), color = dtype_col))$color

# Note that print() is needed inside a loop!!
for (rdim in 1:length(meta_leu)){
  rdimx <- as.data.frame(meta_leu[[names(meta_leu)[rdim]]])
  rdimx$disease_type <- as.factor(sampleinfo_r2$disease_type)
  pdf(glue::glue("/cluster/home/yjliu_jh/share/{names(meta_leu)[rdim]}_leu_disease_typex.pdf"), width = 7, height = 8)
  print(ggscatter(rdimx, x = colnames(rdimx)[1], y = colnames(rdimx)[2], color = "disease_type",
            palette = dtype_col, size = 0.2))
  dev.off()
}




#plot(unlist(rdimx[colnames(rdimx)[1]]), unlist(rdimx[colnames(rdimx)[2]]), col = dtype2, 
#     cex = 0.2, xlab = colnames(rdimx)[1], ylab = colnames(rdimx)[2])
#legend("bottomright",
#       legend = names(dtype_col),
#       fill = dtype_col, cex = 0.6)





dummy1 <- data.frame(x = 1:10, y = 1:10)
pdf("/cluster/home/yjliu_jh/share/testplotw.pdf", width = 7, height = 8)
ggscatter(dummy1, x = "x", y = "y")
dev.off()




