
# the function could be slow when number of groups increase!
# need improve: first filter based on foldchange?
comp_top <- function(x, y = "groups", z = "value", top = 5){
  x <- as.data.frame(x)
  comparisons <- combn(x[, y], 2, simplify=FALSE)
  filtered <- list()
  pvalues <- numeric()
  for (i in 1:length(comparisons)) {
    pair <- comparisons[[i]]
    comp <- paste0(pair, collapse = "_")
    ind1 <- x[, y] %in% pair[1]
    ind2 <- x[, y] %in% pair[2]
    tmp <- as.data.frame(rbind(x[ind1, ], x[ind2, ]))
    tmp[, y] <- factor(tmp[, y], levels=unique(tmp[, y]))
    res <- t.test(x[ind1, z], x[ind2, z])
    pvalues[i] <- res$p.value
    if (res$p.value < 0.05) {
      filtered[[comp]] <- pair
      pvalues <- c(pvalues, res$p.value)
    }
  }
  if(length(filtered)){
    pvalues <- sort(pvalues)
    ind <- pvalues[1:min(top, length(filtered))]
    return(filtered[ind])
  } else {
    return(NULL)
  }
}


gp_violin2 <- function(data, gene, anno, z = 0.1,
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
    comp <- comp_top(long_exp)
    pb <- pb + scale_y_continuous(expand = expansion(mult = c(0, z))) +
      stat_compare_means(comparisons = comp, na.rm = T, method = "t.test", size = 3, hjust = 0.2)

  return(pb)
}








