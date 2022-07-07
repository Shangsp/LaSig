setwd("D:\\5.Research\\02Lactic\\01Data")
mydata <-read.table("TCGA-LUADgeneFPKMExp.txt",header=T,sep='\t',row.names=1,check.names=F)

LAC <- read.table("D:\\5.Research\\02Lactic\\01Data\\LACTATE_gene.txt",sep='\t',check.names=F)
genename <- as.character(LAC[,1])
LACExp <- log2(mydata[as.character(genename),]+1)

group <- ifelse(as.numeric(substr(colnames(LACExp),14,15))<10,"Cancer","Normal")
mergedata <- data.frame(t(LACExp),group)

#difgene
result <- matrix(,25,1)
for (i in 1:25){
  result[i,1] <- wilcox.test(mergedata[,i]~group,mergedata)$p.value
}
rownames(result) <- colnames(mergedata)[1:25]

setwd("D:\\5.Research\\02Lactic\\07difgene")
write.table(result,"LAC_difgene.txt",quote=F,sep='\t',col.names=F)

difgene <- rownames(result)[which(result[,1]<0.05)]

#ggplot2
difgeneExp <- mergedata[,c(difgene,"group")]
library(reshape2)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)

alldata <- melt(difgeneExp,id.vars="group",variable.name="gene",value.name="Expression")
alldata[,1]<-factor(alldata[,1])

pdf("difgene.pdf",width=7,height=4)
p <- ggplot(alldata,aes(x=gene,y=Expression,fill=group))+ 
  geom_boxplot()+
  scale_fill_manual(values=c("#EE4000","#4292C6"))+
  geom_boxplot(outlier.size=0.1,outlier.alpha=0.1)+
  theme_classic()+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size = 10, 
                                   color = "black", face = "plain", 
                                   vjust = 1, hjust = 1,
                                   angle = 60))+
  theme(legend.spacing.x=unit(1.2,"cm")) + theme(axis.ticks.x = element_blank())+
  stat_compare_means(method="wilcox.test",label="p.signif")
print(p)
dev.off()


#cor

library(Hmisc)
cor_res <- rcorr(t(LACExp))

res <- matrix(,625,4)
m=1
for (i in 1:25){
  for (j in 1:25){
    if (j >i){
      res[m,1] <- rownames(cor_res$r)[i]
      res[m,2] <- colnames(cor_res$r)[j]
      res[m,3] <- cor_res$r[i,j]
      res[m,4] <- cor_res$P[i,j]
      m <- m+1
    }
  }
}

colnames(res) <- c("gene1","gene2","cor","pvalue")
res <- na.omit(res)
write.table(res,"D:\\5.Research\\02Lactic\\07difgene\\LAC_cor_Result.txt",quote=F,sep='\t',row.names=F)

mySur <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUAD.survival.tsv",header=T,sep='\t',row.names=1,check.names=F)
name1 <- intersect(rownames(mySur),colnames(LACExp))

Surdata <- data.frame(t(LACExp)[name1,],mySur[name1,c(3,1)])
result <- matrix(,25,6)
library(survival)
for (i in 1:25){
c <- coxph(Surv(OS.time,OS)~Surdata[,i],Surdata)
result[i,2] <- summary(c)[[7]][[1]][1]
result[i,3] <- summary(c)[[7]][[2]][1]
result[i,4] <- summary(c)[[8]][3]
result[i,5] <- summary(c)[[8]][4]
result[i,6] <- summary(c)[[7]][[5]][1]
}
result[,1] <- colnames(Surdata)[1:25]
colnames(result) <- c("Gene","coef","HR","lower.95","upper.95","p")
write.table(result,"D:\\5.Research\\02Lactic\\07difgene\\LACcoxResult.txt",quote=F,sep='\t',row.names=F)

###cluster
setwd("D:\\5.Research\\02Lactic\\02cluster")
mydata <- read.table("D:\\5.Research\\02Lactic\\01Data\\LAClogExp.txt",header=T,sep='\t',row.names=1,check.names=F)
cancerSample <- colnames(mydata)[which(substr(colnames(mydata),14,15)<10)]
cancerData <- mydata[,cancerSample]

cancerData <-data.matrix(cancerData)
library(ConsensusClusterPlus)

res <- ConsensusClusterPlus(cancerData, maxK = 8, reps = 1000, 
                            pItem = 0.8, distance = "euclidean",
                            pFeature = 1, clusterAlg = "pam",
                            seed=123456, title = "LAC_cluster_pam8",
                            writeTable=T,plot = "pdf")

resICL = calcICL(res,title="LAC_cluster_pam8",writeTable=T,plot = "pdf")


clustergroup <- read.table("D:\\5.Research\\02Lactic\\02cluster\\LAC_cluster_pam8\\LAC_cluster_pam8.k=2.consensusClass.csv",sep=',',row.names=1)
mySur <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUAD.survival.tsv",header=T,sep='\t',row.names=1,check.names=F)
Surdata <- mySur
Surdata[which(mySur$OS.time>=3650),1]=0
Surdata[which(mySur$OS.time<3650 & mySur$OS==1),1]=1
Surdata[which(mySur$OS.time>=3650),3]=3650
sampleName <- intersect(rownames(clustergroup),rownames(Surdata))
data1 <- data.frame(Surdata[sampleName,c(3,1)],group=clustergroup[sampleName,])

library(survminer)
library(survival)
library(ggplot2)

gene1 <- survfit(Surv(OS.time, OS)~group, data=data1) #OS.time是生存时间的列名,OS是生存结局的列名,group是分的亚型

pdf("LAC_cluster_survival10.pdf",width=7,height=7,onefile=F)
ggsurvplot(gene1, conf.int=F, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("cluster1","cluster2"),legend.title="",  #High Risk,Low Risk是group里面分的类型
           palette=c("#E8B34B", "#5F5CA5"), 
           title="", 
           risk.table.height=0.3)
dev.off()