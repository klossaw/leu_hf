


# basics: combination of two columns 
# first based on one column, and check for each group if they are assigned to the same group2


dummydata <- data.frame(group1 = sample(letters[1:6], 100, T), group2 = sample(LETTERS[1:6], 100, T))
dummydata$join <- paste(dummydata$group1, dummydata$group2, sep = "-")
dummydata %>% group_by(group1) %>% mutate()



getMethod("get_classes","ConsensusPartitionList")





if (verbose) 
  cat("------------------------------------------------------------\n")
if (verbose) 
  cat("* adjust class labels according to the consensus classifications from all methods.\n")
reference_class = lapply(lt[[1]]@k, function(k) {
  cl_df = get_consensus_from_multiple_methods(res_list, 
                                              k)
  return(cl_df)
})
names(reference_class) = as.character(lt[[1]]@k)
if (verbose) 
  cat("  - get reference class labels from all methods, all k.\n")
rc = reference_class[[1]]$class_df$class
all_k = lt[[1]]@k
for (i in seq_along(all_k)[-1]) {
  class_df = reference_class[[i]]$class_df
  class = class_df[, "class"]
  map = relabel_class(class, rc, full_set = 1:(all_k[i]))    ## see below
  l = which((duplicated(map) | duplicated(map, fromLast = TRUE)) & 
              map != names(map))
  unmapped = setdiff(names(map), map)
  if (any(l)) {
    map[l] = unmapped
  }
  map2 = structure(names(map), names = map)
  reference_class[[i]]$class_df$class = as.numeric(map[as.character(class)])
  reference_class[[i]]$membership = reference_class[[i]]$membership[, 
                                                                    as.numeric(map2[as.character(1:all_k[i])])]
  colnames(reference_class[[i]]$membership) = paste0("p", 
                                                     1:all_k[i])
  rc = reference_class[[i]]$class_df$class
}
res_list@consensus_class = reference_class





reference_class = lapply(lt[[1]]@k, function(k) {
  cl_df = get_consensus_from_multiple_methods(res_list, k)
  return(cl_df)
})


x <- matrix(sample(1:5, 100, replace = T), nrow = 20)



cola::relabel_class(x[, 1], x[, 2])



