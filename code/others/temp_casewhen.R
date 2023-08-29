


quick_check2 <- function(dat, group){
  head(sort(table(dat[dat$subgroups %in% group, ]$fusion), decreasing = T), 9) / sum(sinfo2$subgroups %in% group)
}
quick_check2(ffil, "G3")  ## only 89%

# hamlet didn't find all PML-RARA fusions?


rls2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/rls_exp.rds")
cla1 <- get_classes(rls2, k = 10)
ffil2 <- left_join(ffil, cla1)
groupings <- sinfo2[, 1:3]
cbg <- left_join(groupings, cla1)

quick_check3 <- function(dat, group){
  head(sort(table(dat[dat$class == group, ]$fusion), decreasing = T), 9) / sum(cbg$class %in% group)
}
quick_check3(ffil2, 1) 


dummydata <- data.frame(col1 = letters[1:6], 
                        col2 = letters[c(1:3, 5:7)])
setdiff_v <- Vectorize(setdiff)
intersect_v <- Vectorize(intersect)
mypaste <- function(x, y){
  z <- intersect_v(x, y)
  # z <- intersect(x, y)
  # paste0(z, "z")
}
dummydata %>% mutate(xx = case_when(col1 == col2 ~ mypaste(col1, col2),
                                    col1 != col2 ~ "x",
                                    TRUE ~ "xx"))

mypaste(dummydata$col1, dummydata$col2)
mypaste(dummydata[3, 1], dummydata[3, 2])





