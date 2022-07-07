library('gplots')
library('limma')
setwd("D:\\5.Research\\02Lactic\\12DEgene")
mydata <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUADgeneFPKMExp.txt",
header=T,sep='\t',row.names=1,check.names=F)
mydata <- log2(mydata+1)

mycount0 <- function(x){
return(length(which(x!=0)))}
result <- apply(mydata,1,mycount0)
genename <- which(result>585*0.2)


group <- read.table("D:\\5.Research\\02Lactic\\10cluster\\LAC_cluster_pam8\\LAC_cluster_pam8.k=2.consensusClass.csv",
sep=',',row.names=1)
samplename <- intersect(colnames(mydata),rownames(group))
datamatrix <- mydata[genename,samplename]
group1 <- group[samplename,]
mygroup <- paste("cluster",group1,sep='')
mygroup <- as.matrix(mygroup)
rownames(mygroup) <- samplename
design <- model.matrix(~0+mygroup)
colnames(design)=levels(factor(mygroup))
rownames(design)=rownames(mygroup)
contrast.matrix<-makeContrasts(paste(c("cluster1","cluster2"),collapse = "-"),levels = design)

##step1

fit <- lmFit(datamatrix,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
length(which(abs(nrDEG[,1])>1 & nrDEG[,5]<0.01))

write.table(nrDEG,"LUADFPKMlimmaResult.txt",sep='\t',quote=F)

###
DEgene <- rownames(nrDEG)[which(abs(nrDEG[,1])>1 & nrDEG[,5]<0.01)]
library("clusterProfiler")

x <- as.character(DEgene)
eg <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
genelist <- unique(eg$ENTREZID)
go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05,keyType = 'ENTREZID')
write.table(go,"difgeneGO.txt",quote=F,sep='\t')
barplot(go,showCategory=10,drop=T)
dotplot(go,showCategory=10)

###筛选其中个别通路
ego2=filter(as.data.frame(go),ID=="GO:0002283" | ID=="GO:0002429"|ID=="GO:0045088"|ID=="GO:0045089"|
ID=="GO:0016064"|ID=="GO:0042110"|ID=="GO:0050852"|ID=="GO:0042098"|ID=="GO:0048002"|ID=="GO:0051054")

y <- new("enrichResult",
      result=ego2)
pdf("GO_dotplot.pdf",width=6,height=7)
dotplot(y)
dev.off()

#KEGG

kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
write.table(kegg,"difgeneKegg.txt",quote=F,sep='\t')
barplot(kegg,showCategory=20,drop=T)
dotplot(kegg,showCategory=10)

##
ego1=filter(as.data.frame(kegg),ID=="hsa05132" | ID=="hsa05166"|ID=="hsa05222"|ID=="hsa05210"|
ID=="hsa04145"|ID=="hsa04144"|ID=="hsa00010"|ID=="hsa05219"|ID=="hsa04142"|ID=="hsa01232")

y <- new("enrichResult",
      result=ego1)
pdf("KEGG_dotplot.pdf",width=6,height=7)
dotplot(y)
dev.off()
###
setwd("D:\\5.Research\\02Lactic\\13RiskScore")
library(caret)
library(lasso)
library(glmnet)
idx2 <- createDataPartition(data2$OS, p = 0.67, list = FALSE)

traindata <- data2[idx2,]
testdata <- data2[-idx2,]
write.table(traindata,"traindata.txt",quote=F,sep='\t')
write.table(testdata,"testdata.txt",quote=F,sep='\t')

coxgene <- result[which(as.numeric(result[,6])<0.05),1]
coxgene <- gsub("-",".",coxgene)
x <- as.matrix(traindata[,coxgene])
y <- Surv(traindata$OS.time,traindata$OS)


cv.fit = cv.glmnet(x,y, family = "cox",nfolds=10,alpha=1)
pdf("cv.fit.pdf",height=5,width=5)
plot(cv.fit)
dev.off()

fit <- glmnet(x,y, family = "cox",alpha=1)
pdf("fit.pdf",height=5,width=5)
plot(fit)
dev.off()
coefficients<-coef(fit,s=cv.fit$lambda.min) #0.09373289
Active.Index<-which(as.numeric(coefficients)!=0)
Active.coefficients<-coefficients[Active.Index]

write.table(as.matrix(coefficients)[Active.Index,],"coef.txt",quote=F,sep='\t')


###DEgene

DEdata <- mydata[DEgene,]
mySur <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUAD.survival.tsv",header=T,sep='\t',row.names=1,check.names=F)
name1 <- intersect(colnames(DEdata),rownames(mySur))

data2 <- data.frame(t(DEdata[,name1]),mySur[name1,c(3,1)])
result <- matrix(,4318,6)
library(survival)
for (i in 1:4318){
c <- coxph(Surv(OS.time,OS)~data2[,i],data2)
result[i,2] <- summary(c)[[7]][[1]][1]
result[i,3] <- summary(c)[[7]][[2]][1]
result[i,4] <- summary(c)[[8]][3]
result[i,5] <- summary(c)[[8]][4]
result[i,6] <- summary(c)[[7]][[5]][1]
}
result[,1] <- rownames(DEdata)
colnames(result) <- c("Gene","coef","HR","lower.95","upper.95","p")
write.table(result,"DEgene_single_cox.txt",quote=F,sep='\t',row.names=F)

setwd("D:\\5.Research\\02Lactic\\13RiskScore")
library(caret)
library(lasso)
library(glmnet)
idx2 <- createDataPartition(data2$OS, p = 0.67, list = FALSE)

traindata <- data2[idx2,]
testdata <- data2[-idx2,]
write.table(traindata,"traindata.txt",quote=F,sep='\t')
write.table(testdata,"testdata.txt",quote=F,sep='\t')

coxgene <- result[which(as.numeric(result[,6])<0.05),1]
coxgene <- gsub("-",".",coxgene)
x <- as.matrix(traindata[,coxgene])
y <- Surv(traindata$OS.time,traindata$OS)

cv.fit = cv.glmnet(x,y, family = "cox",nfolds=10,alpha=1)
fit <- glmnet(x,y, family = "cox",alpha=1)
coefficients<-coef(fit,s=cv.fit$lambda.min) #0.09373289
Active.Index<-which(as.numeric(coefficients)!=0)
Active.coefficients<-coefficients[Active.Index]

write.table(as.matrix(coefficients)[Active.Index,],"coef.txt",quote=F,sep='\t')



coef <- read.table("coef.txt",row.names=1,header=T,check.names=F,sep='\t')
trainScore <- colSums(t(traindata[,rownames(coef)])*coef[,1])
trainmerge <- data.frame(score=trainScore,traindata[,c("OS.time","OS")])

trainmerge$OS.time <- trainmerge$OS.time/30
trainmerge$group <- ifelse(trainmerge$score>median(trainmerge$score),"HighRisk","LowRisk") #0.150205

write.table(trainmerge,"trainmerge.txt",quote=F,sep='\t')
library(survminer)
library(survival)
library(ggplot2)
gene1 <- survfit(Surv(OS.time, OS)~group, data=trainmerge) #OS.time是生存时间的列名,OS是生存结局的列名,group是分的亚型
pdf("train_survival.pdf",height=5,width=5,onefile=F)
p <- ggsurvplot(gene1, conf.int=F, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("High Risk","Low Risk"), legend.title="",  #High Risk,Low Risk是group里面分的类型
           palette=c("dodgerblue2", "orchid2"), 
           title="", 
           risk.table.height=0.25)
print(p)
dev.off()


testScore <- colSums(t(testdata[,rownames(coef)])*coef[,1])
testmerge <- data.frame(score=testScore,testdata[,c("OS.time","OS")])

testmerge$OS.time <- testmerge$OS.time/30
testmerge$group <- ifelse(testmerge$score>median(testmerge$score),"HighRisk","LowRisk") #0.1558242

write.table(testmerge,"testmerge.txt",quote=F,sep='\t')
gene2 <- survfit(Surv(OS.time, OS)~group, data=testmerge) #OS.time是生存时间的列名,OS是生存结局的列名,group是分的亚型
pdf("test_survival.pdf",height=5,width=5,onefile=F)
p <- ggsurvplot(gene2, conf.int=F, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("High Risk","Low Risk"), legend.title="",  #High Risk,Low Risk是group里面分的类型
           palette=c("dodgerblue2", "orchid2"), 
           title="", 
           risk.table.height=0.25)
print(p)
dev.off()




###GEO验证集
setwd("D:\\5.Research\\02Lactic\\13RiskScore")
mydata <- read.table("D:\\5.Research\\02Lactic\\01Data\\GEO\\GSE19188\\GSE19188matrix.txt",
header=T,sep='\t',row.names=1,check.names=F)

coef <- read.table("coef.txt",row.names=1,header=T,check.names=F,sep='\t')

mygroup <- read.table("D:\\5.Research\\02Lactic\\01Data\\GEO\\GSE19188\\GSE19188group.txt",
header=T,sep='\t',row.names=1,check.names=F)

data1 <- mydata[,rownames(mygroup)]

GPL570 <- read.table("D:\\5.Research\\02Lactic\\01Data\\GEO\\GPL570ID2gene.txt",
header=T,sep='\t',row.names=1,check.names=F)

data2 <- data.frame(gene=GPL570[rownames(data1),1],data1)

data3 <- data2[which(data2[,1] %in% rownames(coef)),]
group=data3$gene
myMeanFun<-function(x){
  tapply(as.double(x),group,mean)  
}
data4=data3[,-1]
mean=apply(data4,2,myMeanFun)


riskScore <- colSums(mean[rownames(coef),]*coef[,1])
GEOmerge <- data.frame(score=riskScore,mygroup[names(riskScore),c("OS.time","OS")])
GEOmerge$group <- ifelse(GEOmerge$score>median(GEOmerge$score),"HighRisk","LowRisk")#0.02116663
write.table(GEOmerge,"GSE19188data.txt",quote=F,sep='\t')

library(survminer)
library(survival)
library(ggplot2)

gene1 <- survfit(Surv(OS.time, OS)~group, data=GEOmerge) #OS.time是生存时间的列名,OS是生存结局的列名,group是分的亚型
pdf("GSE19188_survival.pdf",height=5,width=5,onefile=F)
p <- ggsurvplot(gene1, conf.int=F, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("High Risk","Low Risk"), legend.title="",  #High Risk,Low Risk是group里面分的类型
           palette=c("dodgerblue2", "orchid2"), 
           title="", 
           risk.table.height=0.25)
print(p)
dev.off()

###risk factor图
setwd("D:\\5.Research\\02Lactic\\13RiskScore")
mydata <- read.table("RiskGeneCacnerData.txt",
header=T,sep='\t',row.names=1,check.names=F)
mySur <- read.table("D:\\5.Research\\02Lactic\\01Data\\TCGA-LUAD.survival.tsv",header=T,sep='\t',row.names=1,check.names=F)
coef <- read.table("D:\\5.Research\\02Lactic\\13RiskScore\\coef.txt",
row.names=1,header=T,check.names=F,sep='\t')
mygroup <- read.table("trainmerge.txt",header=T,sep='\t',row.names=1,check.names=F)
#write.table(mygroup,"trainmerge.txt",quote=F,sep='\t')
name1 <- intersect(rownames(mygroup),rownames(mySur))
mygroup <- mygroup[name1,]
mygroup[,c("OS.time","OS")] <- mySur[name1,c("OS.time","OS")]
samplename <- intersect(colnames(mydata),rownames(mygroup))


library(pheatmap)
library(ggplot2)
library(ggplotify)
library(cowplot)

mygroup <- mygroup[order(as.numeric(mygroup$score)),]
group <- data.frame(mygroup[,4]) #分组信息
colnames(group) <- "group"
rownames(group) <- rownames(mygroup)
data1 <- mydata[rownames(coef),rownames(group)]
p1=pheatmap(data1, 
            show_colnames = F ,
            scale="row",
            fontsize_row=10, 
            color =colorRampPalette(c("blue", "white","red"))(256),
            legend = T,
            show_rownames = T,
		cluster_cols = F,
            cluster_rows = T,
		annotation_col = group)


p1=ggplotify::as.ggplot(p1)

str(mygroup$OS) 
mygroup$OS=factor(mygroup$OS,labels = c("Alive","Dead"))
p2<-ggplot(data=mygroup) + geom_point(aes(x=seq(0:338),y=OS.time/365,color=OS))+scale_x_continuous(breaks=seq(0, 338, 100))+
scale_y_continuous(breaks=seq(0, 21, 3)) #x轴为样本数，样本总数为248

p3=p2+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
axis.line = element_line(colour = "black"))+labs(x="",y="Follow up years")+ 
theme(axis.line = element_line(size=1, colour = "black"))+ geom_vline(aes(xintercept=169), 
colour="#BB0000", linetype="dashed")

biomarker_risk<-mygroup[order(as.numeric(mygroup$score) ),]

biomarker_risk$group<-ifelse(biomarker_risk$score>=median(biomarker_risk$score),"high_risk","low_risk") #按照中值等高低风险组
table(biomarker_risk$group)
p4<-ggplot(data=biomarker_risk) + geom_point(aes(x=seq(0:338),y=score,color=group) )+scale_x_continuous(breaks=seq(0, 338, 100))+
scale_y_continuous(breaks=seq(0, 0.4, 0.1))+ geom_vline(aes(xintercept=169), colour="#BB0000", linetype="dashed") #244为高低风险组时样本的阈值

p5=p4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
axis.line = element_line(colour = "black"))+ labs(x="",y="Risk score")+ theme(axis.line = element_line(size=1, colour = "black"))


library(cowplot)
pdf("risk.pdf",width=7,height=7)
#plot_grid( plot_grid(p3, p5, ncol = 1,axis = 'lrtb',align = 'hv',labels = c("A", 'B')),p1,ncol = 1,axis = 'lrtb',labels = c("", 'C'))

plot_grid( p3, p5, p1,ncol = 1,axis = 'l',align = 'v',labels = c("A", 'B','C'))
dev.off()

