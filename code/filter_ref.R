code.dir <- "~/projects/all/code"
fsetting <- paste(code.dir,"/utils/setting.R",sep="/")
futils <- paste(code.dir,"/utils/utils.R",sep="/")
source(fsetting)
source(futils)

workdir <- "~/projects/all/analysis/mpipeup"
setwd(workdir)
dat <- read.table("positions/list",header=T)
samples <- as.character(unique(dat[,1]))
for (s in samples)
{
  fout <- paste("positions/",s,"_positionlist.txt",sep="")
  fil <- dat[,1] == s
  if(file.exists(fout))
  {
    info("File existed: ",fout)
  }else
  {
    write.table(dat[fil,c(2,3)],fout,sep="\t",quote=F,col.names=F,row.names=F)
  }
}

mpileup <- function(fpos,fbam,fout)
{
  cmd=paste("samtools mpileup -Q25 -q26 -l ",fpos," -f /u2/db/hg19/hg19.fa ",fbam," > ",fout,sep="")
  info(cmd)
  system(cmd)
}

matchstr <- function (string,type)
{
  tmp <- unlist(strsplit(as.character(string),""))
  if(type=="any")
  {
    return (sum(grepl("[a-z]",tmp)) >0 | sum(grepl("[A-Z]",tmp)) >0)
  }else
  {
    return (sum(grepl("[a-z]",tmp)) >0 & sum(grepl("[A-Z]",tmp)) >0)
  }
}

percentFil <- function (string,compare,thred)
{
  tmp <- unlist(strsplit(as.character(string),""))
  total <- sum(grepl("\\.",tmp)) + sum(grepl("\\,",tmp)) + sum(grepl("[a-z]",tmp)) + sum(grepl("[A-Z]",tmp))
  percent <- (sum(grepl("[a-z]",tmp)) + sum(grepl("[A-Z]",tmp)))/total
  if(compare=="more")
  {
    return(percent>thred)
  }else
  {
    return(percent<thred)
  }
  
}

dat <- NULL
for (s in seq(1,25))
{
  if(s<10){
    s <- paste("0",s,sep="")
  }
  for (type in c("A","C"))
  {
    fn <- paste("output/L",s,type,".txt",sep="")
    fbam <- paste("~/projects/all/gatk/L",s,type,".realign.bam",sep="")
    fpos <- paste("positions/",s,"_positionlist.txt",sep="")
    if(file.exists(fn))
    {
      info(fn)
    }else
    {
      mpileup(fpos,fbam,fn)
    }
  }
  fA <- paste("output/L",s,"A.txt",sep="")
  fC <- paste("output/L",s,"C.txt",sep="")
  dat.tmpA <- read.table(fA)
  dat.tmpC <- read.table(fC)
  filA <- sapply(dat.tmpA[,5],function(x)percentFil(x,"more",0.1))
  filC <- sapply(dat.tmpC[,5],function(x)percentFil(x,"less",0.02))
  filMatchA <- sapply(dat.tmpA[,5],function(x)matchstr(x,"and"))
  fil <- filA & filC & dat.tmpA[,4] >10 & filMatchA & dat.tmpC[,4] >10
  dat.tmpA$Sample <- paste("L",s,sep="")
  if(is.null(dat))
  {
    dat <- dat.tmpA[fil,]
  }else
  {
    dat <- rbind(dat,dat.tmpA[fil,])
  }
}
fout <- "snvFilter.csv"
dat <- dat[,c(7,1,2,3,4,5)]
write.csv(dat,fout,quote=T,row.names=F)

datMerge <- read.csv("~/projects/all/mutectFilter/compare-snv.csv",row.names=1)

dat.out <- merge(dat[,c(1,2,3)],datMerge,by=c(1,2,3))
fout <- "needReCheck.csv"
write.csv(dat.out,fout,quote=T,row.names=F)