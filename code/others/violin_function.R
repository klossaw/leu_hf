pkgs <- c("ggpubr", "ggthemes", "jhtools", "glue", "ggsci",
          "patchwork", "tidyverse")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# quick violin function
# required input: standard gene-exp profile, two-column group annotation file, ...

gp_violin <- function(data, gene, anno, z = 0.1, comp = NULL,
                      col = RColorBrewer::brewer.pal(n = 12, name = "Paired")){
  
  stopifnot("two columns (sample_id and groups) should be contained in anno" = c("sample_id", "groups") %in% colnames(anno))
  selected_exp <- data[data$gene_name %in% gene, ]
  long_exp <- selected_exp[, -1] %>% pivot_longer(-gene_name, names_to = "sample_id", values_to = "value")
  long_exp <- left_join(long_exp, anno)
  
  pb <- ggpubr::ggviolin(long_exp[long_exp$gene_name %in% gene, ], x = "groups", y = "value",
                         color = "groups", palette = col, add = "jitter",
                         title = gene, xlab = F, ylab = F, size = 1,
                         ggtheme = theme(panel.background = element_blank(),
                                   panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",
                                                                                                size = rel(1)), legend.key = element_blank(),
                                   strip.background = element_rect(fill = "white", colour = "black",
                                                                   size = rel(2)), complete = TRUE,
                                   text = element_text(size = 13),
                                   plot.title = element_text(size = 12, hjust = 0.02),
                                   legend.position = "none")) 
  if (!is.null(comp)){
    pb <- pb + scale_y_continuous(expand = expansion(mult = c(0, z))) +
               stat_compare_means(comparisons = comp, na.rm = T, method = "t.test", size = 3, vjust = 0.2)
  }
  return(pb)
}

# load data and construct data structure
exp_s <- read_csv(glue::glue("{leuloc}/exp/tables/linxiangjie_human.csv"))
sinfo2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_20220809.xlsx")
grouping_anno <- sinfo2[, 1:2]
colnames(grouping_anno) <- c("sample_id", "groups")
grouping_anno$groups <- factor(grouping_anno$groups, levels = paste0("G", 1:10))
interest_genes <- c("IGF2BP3", "IGF2BP2", "DDX21", "HDAC9", "ACTN1", "FLOT1", "SLC7A11")



# test plots
for (i in 1:length(interest_genes)){
  p1 <- gp_violin(exp_s, interest_genes[i], grouping_anno)
  ggsave(glue::glue("/cluster/home/yjliu_jh/projects/leu_j/output/violin_{interest_genes[i]}.pdf"),
         p1, dpi=300, width = 7, height = 4, units = "in")
}

# new annotation groups based on KMT2A fusion
#anno_f <- sinfo2[, c("sample_id", "subgroups", "fusion")]
#anno_f$groups <- "No KMT2A fusion"
#anno_f$groups[grepl("KMT2A", anno_f$fusion)] <- "KMT2A fusion"
#two_comp <- list( c("No KMT2A fusion", "KMT2A fusion") )


# KMT2A
#px <- list()
#for (i in 1:length(interest_genes)){
#  px[[interest_genes[i]]] <- gp_violin(exp_s, interest_genes[i], anno_f, comp = two_comp)
#}

# groups
px <- list()
for (i in 1:length(interest_genes)){
  px[[interest_genes[i]]] <- gp_violin(exp_s, interest_genes[i], grouping_anno)
}

# text annotation
text = paste("\n G1: PML-RARA fusion\n",
             "   G2: CBFB-MYH11 fusion\n",
             "   G3: RUNX1-RUNX1T1 fusion\n", 
             "   G4: CEBPA (bi-) mutation")
pt <- ggplot() + 
  annotate("text", x = 4, y = 25, size=8, label = text) + 
  theme_void()


#  manually merge into 1 plot
line1 <- ggarrange(px[[1]], px[[2]], px[[3]], ncol = 3, nrow = 1)
line1 <- annotate_figure(line1, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
line2 <- ggarrange(px[[4]], px[[5]], px[[6]], ncol = 3, nrow = 1)
line2 <- annotate_figure(line2, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
line3 <- ggarrange(px[[7]], pt, ncol = 3, nrow = 1)
line3 <- annotate_figure(line3, left = text_grob(bquote(log[2] ~ (FPKM)), rot = 90, vjust = 1))
merged <- ggarrange(line1, line2, line3, ncol = 1, nrow = 3)
ggsave("/cluster/home/yjliu_jh/projects/leu_j/output/violin_7_new.pdf", merged, dpi=300,
       width = 15, height = 9, units = "in")


# quick filter
anno_plt <- sinfo2[, c("sample_id", "plt")]
anno_plt$plt <- as.numeric(anno_plt$plt)
selected_exp <- exp_s[exp_s$gene_name %in% interest_genes[i], ]
long_exp <- selected_exp[, -1] %>% pivot_longer(-gene_name, names_to = "sample_id", values_to = "gene_expression")
anno_plt <- anno_plt %>% left_join(long_exp) %>% na.omit()
anno_plt <- anno_plt[anno_plt$plt < 210, ]


tmp <- ggscatter(anno_plt, x = "plt", y = "gene_expression", color = "#00AFBB",
          add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), 
          conf.int = TRUE, cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))

ggsave(glue::glue("/cluster/home/yjliu_jh/projects/leu_j/output/scatter_plt_{interest_genes[i]}.pdf"),
       tmp, dpi=300, width = 5, height = 4, units = "in")

