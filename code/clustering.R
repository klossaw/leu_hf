pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load annotation data
leuloc <- "/cluster/home/jhuang/projects/leukemia/analysis/linxiangjie/human/rnaseq"
hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/hgnc_complete_set_2022-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id")]
colnames(hugo_anno)[3] <- "gene_id"


# (use cola on the small subset of leukemia)
count_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human_counts.csv"))
exp_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human.csv"))
s_fil <- exp_s[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
s_fil <- s_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
s_fil <- adjust_matrix(as.matrix(s_fil))

ssd <- transform(s_fil, SD=apply(s_fil, 1, sd, na.rm = TRUE))
s_fil2 <- s_fil[order(ssd$SD, decreasing = T), ][1:3000, ]
# stuck at CV-skmeans  500 rows are still too hard for the algorithm. 250 top rows (from 3000 rows) are at regular speed

rls <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/rls_exp.rds")
cola_report(rls, output_dir = "/cluster/home/yjliu_jh/projects/leu_j/output/")

cola_9 <- get_classes(rls, k = 9)
cola_9$sample_id <- rownames(cola_9)
sinfo2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_20220809.xlsx")
compcl <- left_join(cola_9, sinfo2[, 1:3])
compcl$j <- paste0(compcl$subgroups, "-", compcl$class)

tmp <- get_matrix(rls)

# check difference? first expression
g4 <- sinfo2[sinfo2$subgroups %in% "G4", "sample_id"]
g4$cluster <- "G4"
g5_npm1 <- sinfo2[sinfo2$characteristic %in% "NPM1" & sinfo2$subgroups %in% "G5", "sample_id"]
g5_npm1$cluster <- "G5"
col_s <- rbind(g4, g5_npm1)

s_fil <- count_s[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
s_fil <- s_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
s_fil <- as.matrix(round(s_fil))   ## as.matrix may produce non-integers? maybe another flaw of floats

mat_s <- s_fil[, colnames(s_fil) %in% col_s$sample_id]

ddss <- DESeqDataSetFromMatrix(countData = round(mat_s),
                               colData = col_s,
                               design = ~ cluster)   

keep <- rowSums(counts(ddss) >= 2) > ceiling(0.2 * ncol(mat_s))
ddss <- ddss[keep,]
keep2 <- rowSums(counts(ddss) == 0) < ceiling(0.5 * ncol(mat_s))
ddss <- ddss[keep2,]
suppressMessages(dds <- DESeq(ddss)) 


resA <- results(dds, c("cluster", "G4", "G5"))
resA$gene_name <- rownames(resA)
resA$absFC <- abs(resA$log2FoldChange)
head(resA[order(resA$padj), ], 20)
diff_genes <- as.character(na.omit(rownames(resA)[resA$padj < 0.05]))


yy <- enrichGO(diff_genes, 'org.Hs.eg.db', ont = "BP", pvalueCutoff = 0.01)
# 	head(as.data.frame(yy))






# sankey plot to compare differences
#

cl2 <- compcl %>% group_by(class, subgroups) %>% summarise(value = n()) %>% as.data.frame()
cl2$class <- as.character(cl2$class)
edges <- cl2
edges$Source <- match(cl2$class, tempnodes$name) - 1
edges$Target <- match(cl2$subgroups, tempnodes$name) - 1

tempnodes <- data.frame(
  name=c(as.character(cl2$class), 
         as.character(cl2$subgroups)) %>% unique()
)

sankeyNetwork(Links = edges, Nodes = tempnodes, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)

# seems sankey plot is not enough for visualization    just continue with no figures 


