

library(data.table)

filename <- "/cluster/home/jhuang/projects/ngs/docs/202209/SampleSheet_2022_9_14.csv"
filename <- "/cluster/home/jhuang/projects/ngs/docs/202209/SampleSheet_2022_9_14.csv.bak"

temp <- fread(filename, skip = "I5_Index_ID")
temp$indexall <- paste0(temp$index, temp$index2)
if(sum(duplicated(temp$indexall))){
  print(paste0("第", which(duplicated(temp$indexall)), "行有序列重复"))
}

if(sum(duplicated(temp$Sample_ID))){
  print(paste0("第", which(duplicated(temp$Sample_ID)), "行有样本名重复"))
}





