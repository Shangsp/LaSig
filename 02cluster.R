setwd("D:\\1.Research\\02Lactic\\01Data")
mydata <- read.table("TCGA-LUAD.htseq_fpkm.tsv",header=T,sep='\t',row.names=1,check.names=F) #读取表达文件,USCS Xena下载的文件为log(data+1)转换过
mydata1 <- 2^mydata-1 #转换文件表达值

myfile <- read.table("Gencodev22Type.txt",header=F,sep='\t',row.names=1,check.names=F)

###比对ensembleID到基因名
data1 <- data.frame(genename=myfile[rownames(mydata1),1],mydata1)

###相同基因取均值
group=data1$genename
myMeanFun<-function(x){
  tapply(as.double(x),group,mean)
}

data2=data1[,-1]
mean=apply(data2,2,myMeanFun)

Result <- data.frame(genename=rownames(mean),mean)
colnames(Result)[2:dim(Result)[2]] <- colnames(mydata)
write.table(Result,"TCGA-LUADgeneFPKMExp.txt",quote=F,sep='\t',row.names=F)


		   
setwd("D:\\1.Research\\02Lactic\\02cluster")
colnames(mean) <- colnames(mydata)
genename <- c("ACTN3","HAGH","HIF1A","LDHA","LDHAL6B","LDHB",
"LDHC","LDHD","MIR210","PARK7","PER2","PFKFB2","PNKD","SLC25A12","C12orf5","TP53")

cancerSample <- colnames(mean)[which(substr(colnames(mean),14,15)<10)]
samplegroup <- ifelse(substr(colnames(mean),14,15)<10,"Cancer","Normal")
cancerData <- mean[genename,cancerSample]

library(pheatmap)
pheatmap(cancerData,)
mygroup <- data.frame(group=clustergroup[sampleName,])

pdf("pheatmap.pdf",width = 8, height = 7)
pheatmap(cancerData, scale = "row",color = colorRampPalette(c("green","white","red"))(256),
annotation_col = mygroup,show_rownames=T,show_colnames=F,cluster_cols=T)
dev.off()

mygroup <- data.frame(group=samplegroup)
rownames(mygroup) <- colnames(mean)
LaData <- mean[genename,c(rownames(mygroup)[which(mygroup[,1]=="Normal")],
rownames(mygroup)[which(mygroup[,1]=="Cancer")])]
pheatmap(LaData, scale = "row",color = colorRampPalette(c("#0000b8","blue","white","red","#8d0000"))(256),
annotation_col = mygroup,show_rownames=T,show_colnames=F,cluster_cols=F)


a <- pheatmap(cancerData,scale="row",color = colorRampPalette(c("#0000b8","blue","white","red","#8d0000"))(256),
annotation_col = mygroup,show_rownames=T,show_colnames=F,cluster_cols=T)

colnames(cancerData)[a$tree_col$order]
LaData <- mean[genename,c(rownames(mygroup)[which(mygroup[,1]=="Normal")],
colnames(cancerData)[a$tree_col$order])]
pdf("lactic_pheatmap.pdf",height=7,width=7)
pheatmap(LaData, scale = "row",color = colorRampPalette(c("#0000b8","blue","white","red","#8d0000"))(256),
annotation_col = mygroup,show_rownames=T,show_colnames=F,cluster_cols=F)
dev.off()

cancerData <-data.matrix(cancerData)
library(ConsensusClusterPlus)

res <- ConsensusClusterPlus(cancerData, maxK = 10, reps = 1000, 
                            pItem = 0.8, distance = "euclidean",
                            pFeature = 1, clusterAlg = "pam",
                            seed=123456, title = "cluster_pam10",
                            writeTable=T,plot = "pdf")

resICL = calcICL(res,title="cluster_pam10",writeTable=T,plot = "pdf")

clustergroup <- read.table("D:\\1.Research\\02Lactic\\02cluster\\cluster_pam10\\cluster_pam10.k=2.consensusClass.csv",sep=',',row.names=1)
mySur <- read.table("D:\\1.Research\\02Lactic\\01Data\\TCGA-LUAD.survival.tsv",header=T,sep='\t',row.names=1,check.names=F)
sampleName <- intersect(rownames(clustergroup),rownames(mySur))
data1 <- data.frame(mySur[sampleName,c(3,1)],group=clustergroup[sampleName,])

library(survminer)
library(survival)
library(ggplot2)

gene1 <- survfit(Surv(OS.time, OS)~group, data=data1) #OS.time是生存时间的列名,OS是生存结局的列名,group是分的亚型

pdf("cluster_survival.pdf",width=7,height=7,onefile=F)
ggsurvplot(gene1, conf.int=F, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("cluster1","cluster2"),legend.title="",  #High Risk,Low Risk是group里面分的类型
           palette=c("#E8B34B", "#5F5CA5"), 
           title="", 
           risk.table.height=0.3)
dev.off()