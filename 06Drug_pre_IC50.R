BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))

mydata <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUADgeneFPKMExp.txt",
header=T,sep='\t',row.names=1,check.names=F)
mydata <- log2(mydata+1)

library(parallel)
library(pRRophetic)
library(ggplot2)

result <- c(rep(0,585))
myDrugName <- list()
n <- 1
for (i in possibleDrugs2016){
  x.int <- try(a <- pRRopheticPredict(testMatrix=as.matrix(mydata),
    drug=i,tissueType = "lung",batchCorrect = "eb",
    selection=1,dataset = "cgp2016"),silent=TRUE)
  if ('try-error' %in% class(x.int)){
    next
  }else{
    a <- pRRopheticPredict(testMatrix=as.matrix(mydata),
      drug=i,tissueType = "lung",batchCorrect = "eb",
      selection=1,dataset = "cgp2016")
    result <- data.frame(result,a)
    myDrugName[n] <- i
    n <- n+1
  }
}
result1 <- result[,-1]

setwd("D:\\5.Research\\02Lactic\\16Drug")
write.table(result1,"DrugIC50_1.txt",quote=F,sep='\t')
write.table(myDrugName,"DrugName.txt",quote=F,sep='\t')