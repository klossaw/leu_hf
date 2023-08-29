sampleinfo2 <- readxl::read_excel("/cluster/home/flora_jh/projects/data_preprocess/0.check_fq/docs/sampleinfo2.xlsx",
                                  col_types = c("text", "text", "numeric", "text", 
                                                "text", "numeric", "text", "text", 
                                                "text", "text", "text"))

sampleinfo2$gender[sampleinfo2$gender %in% c("missing", "NA")] <- NA
tempages <- sampleinfo2$age[!is.na(sampleinfo2$age)][sampleinfo2$age[!is.na(sampleinfo2$age)] > 100] 
sampleinfo2$age[!is.na(sampleinfo2$age)][sampleinfo2$age[!is.na(sampleinfo2$age)] > 100] <- tempages / 365
sampleinfo2$rna_library[sampleinfo2$rna_library %in% c("unknown", "NA")] <- NA

# why not base R?  no  should be write.csv
write.table(sampleinfo2, "/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo2.csv", 
            na = "", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")












# what... what

his_genes1 <- c("H4C11",  "H2BC3",  "H4C4",   "H3C12",  "H2BC14", "H4C6",   "H2AC16", "H2AC11",
                "H2AC12", "H1−5",   "H3C2",   "H3C3",   "H3C11",  "H2AC17", "H2AC13", "H2BC13",
                "H2AC4",  "H2AC14", "H4C9",   "H4C1",   "H4C2",   "H4C13",  "H3C7",   "H2BC9",
                "H2BC10", "H3C8",   "H3C13",  "H2BC7",  "H2BC18", "H2BC15", "H4C14",  "H4C8",
                "H2BC17", "H3C10",  "H2BC11", "H4C12",  "H3C15",  "H2BC8",  "H2AC8",  "H2BC6",
                "H4C5",   "H2AC15", "H1−3",   "H3C4",   "H4−16",  "H2AC20", "H2AC21", "H3C1",
                "H2BC5",  "H1−4",   "H2AC7",  "H2BC4",  "H1−2",   "H2AC6")

gender_genes <- c("KDM5D", "UTY", "USP9Y", "DDX3Y", "RPS4Y1", "ZFY", "EIF1AY", "ENSG00000229807", "XIST")





