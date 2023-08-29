cbg <- left_join(groupings, cla1)
igraph::compare(cbg$class, as.numeric(as.factor(cbg$subgroups)), method = "nmi")


cbg$class <- paste0("X", cbg$class)
cbg$comb <- paste0(cbg$subgroups, "-", cbg$class)
cbg %>% group_by(subgroups) %>% mutate(ac = n()) %>% group_by(subgroups, class) %>% 
  summarise(percentage = n() / ac) %>% unique() %>% as.data.frame()
# clusters G1 G2 G3 G4 G8 had high concordance, and these clusters are EXACTLY 
# the clusters who have distinct genomic patterns!
cbg %>% group_by(class) %>% mutate(ac = n()) %>% group_by(class, subgroups) %>% 
  summarise(percentage = n() / ac) %>% unique() %>% as.data.frame()


# test if count-based clustering is comparable with exp-based
#rls1 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/rls.rds")
rls2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/rls_exp.rds")
cla1 <- get_classes(rls2, k = 10)
cla2 <- get_classes(rls2, k = 9)
cla3 <- get_classes(rls2, k = 8) 
cla1$sample_id <- rownames(cla1)
cla2$sample_id <- rownames(cla2)
cla3$sample_id <- rownames(cla3)
# clusters using exp obviously better than count in concordance, according to tests (not shown here)
# note that it's single source data, batch correction is not needed





# now test if exp from salmon and deseq-version exp can have affects on clustering

count_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human_counts.csv"))
s_fil <- count_s[, -2] %>% left_join(hugo_anno[hugo_anno$locus_group %in% c("protein-coding gene"), ]) %>% na.omit()
s_fil <- s_fil %>% dplyr::select(-c(gene_id, locus_group)) %>% as.data.frame() %>%
  remove_rownames() %>% column_to_rownames(var = "symbol")
s_fil <- round(as.matrix(s_fil))


# construct dds object and get DESeq-normalized matrix
ddss <- DESeqDataSetFromMatrix(countData = s_fil,
                               colData = groupings,
                               design = ~ subgroups) 
vsds <- vst(ddss, blind = FALSE)
normalizedExps <- assay(vsds)
readr::write_rds(normalizedExps, "/cluster/home/yjliu_jh/projects/leu_j/data/exp_norm_small.rds")

sa_mat <- adjust_matrix(normalizedExps)
sasd <- transform(sa_mat, SD=apply(sa_mat, 1, sd, na.rm = TRUE))
sa_mat2 <- sa_mat[order(sasd$SD, decreasing = T), ][1:3000, ]

rlsa = run_all_consensus_partition_methods(s_fil2, max_k = 12,
                                           top_value_method = c("SD", "MAD", "CV", "ATC"),
                                           partition_method = c("hclust", "kmeans", "skmeans", "mclust", "pam"),
                                           cores = 20)
readr::write_rds(rlsa, "/cluster/home/yjliu_jh/projects/leu_j/data/rls_exp_adjusted.rds")

s_10 <- get_classes(rlsa, k = 10)
colnames(s_10)[1] <- "class_alter"
s_10$sample_id <- rownames(s_10)
cbg <- left_join(cbg, s_10[, c("sample_id", "class_alter")])



s_8 <- get_classes(rlsa, k = 8)
colnames(s_8)[1] <- "class_alter"
s_8$sample_id <- rownames(s_8)
cbg8 <- left_join(groupings, s_8[, c("sample_id", "class_alter")])
cl8 <- cbg8 %>% group_by(class_alter, subgroups) %>% summarise(value = n()) %>% as.data.frame()
cl8$class_alter <- paste0("C", as.character(cl8$class_alter))

nodes8 <- data.frame(name = c(cl8$class_alter, cl8$subgroups) %>% unique())
edges8 <- cl8
edges8$Source <- match(cl8$class_alter, nodes$name) - 1
edges8$Target <- match(cl8$subgroups, nodes$name) - 1

sankeyNetwork(Links = edges8, Nodes = nodes8, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)



#> igraph::compare(cbg$class_alter, as.numeric(as.factor(cbg$subgroups)), method = "nmi")
#[1] 0.6457457
#> igraph::compare(cbg$class_alter, cbg$class, method = "nmi")
#[1] 0.7973528

# get_stats(rlsa, 10)

## use heatmaps? and a sankey plot to draw the figure 

cl3 <- cbg %>% group_by(class_alter, subgroups) %>% summarise(value = n()) %>% as.data.frame()
cl3x <- cbg %>% group_by(class, class_alter) %>% summarise(value = n()) %>% as.data.frame()

cl3$class_alter <- paste0("C", as.character(cl3$class_alter))
cl3x$class_alter <- paste0("C", as.character(cl3x$class_alter))
cl3x$class <- as.character(cl3x$class)

nodes <- data.frame(name = c(cl3$class_alter, cl3$subgroups) %>% unique())
edges3 <- cl3
edges3$Source <- match(cl3$class_alter, nodes$name) - 1
edges3$Target <- match(cl3$subgroups, nodes$name) - 1

sankeyNetwork(Links = edges3, Nodes = nodes, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)

nodes3x <- data.frame(name = c(cl3x$class_alter, cl3x$class) %>% unique())
edges3x <- cl3x
edges3x$Source <- match(cl3x$class_alter, nodes3x$name) - 1
edges3x$Target <- match(cl3x$class, nodes3x$name) - 1

sankeyNetwork(Links = edges3x, Nodes = nodes3x, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)


nodesa <- data.frame(name = c(cl3$subgroups, cl3x$class_alter, cl3x$class) %>% unique())
edgesa <- rbind(cl3x, setNames(cl3, names(cl3x)))
edgesa$Source <- match(edgesa$class_alter, nodesa$name) - 1
edgesa$Target <- match(edgesa$class, nodesa$name) - 1

sankeyNetwork(Links = edgesa, Nodes = nodesa, Source = "Source",
              Target = "Target", Value = "value", NodeID = "name",
              sinksRight=FALSE)


# instead of using get_classes(), retrieve from best methods! 

stats10 <- get_stats(rlsa, 10)
stats10 <- stats10[order(stats10[, "1-PAC"], decreasing = T), ]
top5_stats10 <- rownames(stats10)[1:5]


s_10_SD_pam <- get_classes(rlsa["SD:pam"], k = 10)
temp_s_10 <- get_classes(rlsa[top5_stats10[1]], k = 10)





# next: find new characteristics for clusters



# construct and merge annotations (itd, fusion, mutation)












# match the

cbg <- left_join(cbg, sinfo2[, c("sample_id", "fusion")])
cbg$characteristic[is.na(cbg$characteristic)] <- "unknown"

quick_check2 <- function(dat, group){
  sort(table(dat[dat$class_alter %in% group, ]$fusion), decreasing = T)
}
quick_check2(cbg, "C3") 

quick_check3 <- function(dat, group){
  sort(table(dat[dat$class_alter %in% group, ]$characteristic), decreasing = T)
}
quick_check3(cbg, "C1") 


# check C1-C4 differences, and find DNA-level characteristics for C3 4 7 9 10




cola_report(rls, output_dir = "/cluster/home/yjliu_jh/projects/leu_j/output/")
cola_9 <- get_classes(rls, k = 9)
cola_9$sample_id <- rownames(cola_9)
sinfo2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_20220809.xlsx")
compcl <- left_join(cola_9, sinfo2[, 1:3])
compcl$j <- paste0(compcl$subgroups, "-", compcl$class)

tmp <- get_matrix(rls)


