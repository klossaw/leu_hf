
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

# after check for annotation cols
x1 <- collected_info[["GSE165656"]][, c("sample_id", "new os",  "vital.status new")]
x2 <- collected_info[["target"]][, c("sample_id", "os_days", "vital_status")]
x3 <- collected_info[["target_aml"]][, c("sample_id", "os_days", "vital_status")]
colnames(x1) <- c("sample_name", "os", "status")
colnames(x2) <- c("sample_name", "os", "status")
colnames(x3) <- c("sample_name", "os", "status")
x1$os <- 12 * x1_os

x <- rbind(x1, x2)
x <- rbind(x, x3)
x$status[x$status == "A"] <- "Alive"
x$status[x$status == "D"] <- "Dead"
x0 <- x[x$status %in% c("Alive", "Dead"), ]

itd_h <- readr::read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/leukemia/itd.csv")
itd_h$length <- itd_h$rose_end_pos - itd_h$rose_start_pos + 1
itd3 <- itd_h[itd_h$itd %in% "flt3", ]
summary(itd3$length)


# 
# itd_n <- readr::read_csv("/cluster/home/yliang_jh/projects/mRNA/leukemia/hamlet_itd.csv")
# itd_n$length <- itd_n$rose_end_pos - itd_n$rose_start_pos + 1
# itd3 <- itd_n[itd_n$itd %in% "flt3", ]



x00 <- inner_join(x0, itd3[, c("sample_name", "length", "fuzziness")])
x00 <- na.omit(x00)
x00 <- x00[order(x00$length, decreasing = T), ]

x01 <- x00[!duplicated(x00$sample_name), ]

#readr::write_rds(x01, "/cluster/home/yjliu_jh/projects/temp/tempx01.rds")

#x01 <- readr::read_rds("/cluster/home/yjliu_jh/projects/temp/tempx01.rds")

x01$len <- ifelse(x01$length > median(x01$length), "long", "short")
x01$len <- ifelse(x01$length > 45, "long", "short")
x01$sta <- ifelse(x01$status == "Alive", 0, 1)

library(survminer)
library(survival)



p_temp <- ggsurvplot(survfit(Surv(os, sta) ~ len, data = x01),
                     data = x01, palette = "jco", risk.table = T, risk.table.height = 0.3,
                     legend.labs=levels(droplevels(as.factor(x01$len))), pval = T)

       
