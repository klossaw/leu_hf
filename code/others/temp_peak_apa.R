pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "data.table")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "ebio"
chr <- "chr14"

f <- function(x){
  p <- vector("numeric", length(x)-1)
  for (i in 2:length(x)){
    p[i] = x[i] - x[i-1]
  }
  return(p)
}

# plot peaks of A80 on chr14:23117306-23117822.
X2 <- 23117306:23117821
v3 <- 23117307:23117822
dat <- data.frame(X2, v3)
sample <- "A80"
wig <- glue::glue("/cluster/home/flora_jh/projects/3UTR/{project}/03.view_peaks/intersect/{chr}/{sample}.wig.bed")
file <- read_tsv(wig, col_names = F, show_col_types = F)
test <- left_join(dat, file, by = "X2") %>% dplyr::select(X2, v3, X4)


#file$X5 <- file$X4 - c(0, file$X4[-length(file$X4)])

colnames(file) <- c("chromosome", "start", "end", "width")
file14 <- file[, 2:3]
file14$ori <- "predict"

pas14 <- pas[pas$chromosome %in% "chr14", c("start", "end")]
pas14$ori <- "database"

comp14 <- range_join(pas14, file14)


temp14 <- as.data.table(file14)[, as.data.table(reduce(IRanges(start, end)))]





for(j in 2:nrow(test)){
  if (is.na(test[j,"X4"])) test[j,"X4"] <- test[j-1,"X4"]
}
setDT(test)
n <- seq(10, 400, 10)
max <- vector("numeric", length(n))
for (i in 1:length(n)){
  avg <- test[, mean(X4), by= (seq(nrow(test)) - 1) %/% n[i]] %>% .$V1
  max[i] <- max(f(avg))
}
plot_dat <- data.frame(n, max)

pdf("A80_height_change_vs_bin.pdf")
ggplot(plot_dat, aes(x = n, y = max)) + xlab("bin size (bp)") + ylab("max height change") + geom_line() + geom_point() + theme_classic()
dev.off()
