library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)

setwd("D:\\5.Research\\02Lactic\\18estimate")
mydata <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUADgeneFPKMExp.txt")
mydata1 <- log2(mydata+1)
filterCommonGenes(input.f="D:\\5.Research\\02Lactic\\01Data\\TCGA-LUADgeneFPKMExp.txt", 
 output.f="TCGA_LUADgenes.gct", 
 id="GeneSymbol")

estimateScore(input.ds = "TCGA_LUADgenes.gct",
 output.ds="TCGALUAD_estimate_score.gct", 
 platform="affymetrix")

plotPurity(scores="OV_estimate_score.gct", samples="s516", 
 platform="affymetrix")

scores=read.table("TCGALUAD_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
write.table(scores,"estimate_score.txt",quote=F,sep='\t')


