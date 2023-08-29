# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "monocle") ## note that monocle3 should not be involved
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# keep newest update of sampleinfo
leu <- readr::read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")
sampleinfo_r <- as.data.frame(colData(leu))
meta_leu <- metadata(leu)
leu_data <- assay(leu, "tpm_remove_rnatype")
leu_data <- as(leu_data, "dgCMatrix")

# ======= monocle2 ========

# construct cds in monocle2 style

fd <- data.frame(gene_short_name = row.names(leu_data), row.names = row.names(leu_data))
pd <- new('AnnotatedDataFrame', data = sampleinfo_r)
fd <- new('AnnotatedDataFrame', data = fd)

cds0 <- monocle::newCellDataSet(leu_data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = VGAM::negbinomial.size())

cds0 <- cds0 %>% estimateSizeFactors() %>% estimateDispersions()

# select DE genes among disease types
cds1 <- monocle::detectGenes(cds0, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds1), num_cells_expressed >= 0.1 * dim(leu_data)[2])) 
# if monocle3 has been loaded, detach it: detach("package:monocle3",unload = TRUE)
diff_test_res <- monocle::differentialGeneTest(cds1[expressed_genes, ],
                                               fullModelFormulaStr = "~subgroups_hjy", cores = 1)

# select 1000 genes 
deg <- subset(diff_test_res, qval < 0.001)
deg <- deg[order(deg$qval, decreasing = F), ]
ordering_genes <- row.names(deg[1:1000,])

# orderingfilter cds object for plot
cds2 <- monocle::setOrderingFilter(cds1[expressed_genes, ], ordering_genes)
cds_test2 <- monocle::reduceDimension(cds2, max_components = 2, method = 'DDRTree')




# ======= run these code when needed ========== 

#source("/cluster/home/yjliu_jh/projects/mef2d/code/order_cells.R")

# "reducedDimW<-" <- function (cds, value)
# {
#   stopifnot(is(cds, "CellDataSet"))
#   cds@reducedDimW <- value
#   validObject(cds)
#   cds
# }
# 
# 
# "reducedDimK<-" <- function (cds, value)
# {
#   stopifnot(is(cds, "CellDataSet"))
#   cds@reducedDimK <- value
#   validObject(cds)
#   cds
# }


# plot
pm1 <- monocle::plot_cell_trajectory(cds_test2, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
pm2 <- monocle::plot_cell_trajectory(cds_test2, color_by = "subgroups_hjy", size = 1, show_backbone = TRUE)

pm_all <- pm1 + pm2 + plot_layout(nrow = 1)

pdf("/cluster/home/yjliu_jh/projects/temp/leu_ddrtree_monocle2_0208.pdf", compress = F,
    width = 15, height = 12)
pm_all
dev.off()


pdf("/cluster/home/yjliu_jh/projects/temp/leu_ddrtree_monocle2_0208_sep2.pdf", compress = F,
    width = 15, height = 12)
pm2 + facet_wrap("~disease_type", nrow = 1)
dev.off()


pdf("/cluster/home/yjliu_jh/projects/temp/leu_ddrtree_monocle2_0208_sep3.pdf", compress = F,
    width = 18, height = 30)
pm2 + facet_wrap("~dataset", nrow = 8)
dev.off()
