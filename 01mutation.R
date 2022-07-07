setwd("D:\\1.Research\\02Lactic\\01Data")
mydata <- read.table("TCGA-LUADgeneFPKMExp.txt",header=T,sep='\t',row.names=1,check.names=F)

genename <- c("ACTN3","HAGH","HIF1A","LDHA","LDHAL6B","LDHB",
"LDHC","LDHD","MIR210","PARK7","PER2","PFKFB2","PNKD","SLC25A12","C12orf5","TP53")

samplegroup <- ifelse(substr(colnames(mydata),14,15)<10,"Cancer","Normal")
LaData <- log2(mydata[genename,]+1)
mergeData <- data.frame(t(LaData),group=samplegroup)
library(reshape2)
mergeData <- melt(mergeData,id.vars="group",variable.name="gene",value.name="value")

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
pdf("boxplot_group.pdf",height=7,width=10)
ggplot(mergeData,aes(x=gene,y=value,fill=group))+
  geom_boxplot(outlier.size=0.1,outlier.alpha=0.1)+
  scale_fill_manual(values=c("#BA2353","#126BAF"))+
  theme_classic() +
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size = 10, 
                                   color = "black", face = "plain", 
                                   vjust = 1, hjust = 1,
                                   angle = 60))+
  theme(legend.spacing.x=unit(1.2,"cm")) + theme(axis.ticks.x = element_blank())+
  stat_compare_means(aes(group=group), label = "p.signif",method = "anova",
  hide.ns = TRUE)
dev.off()

pdf("cluster_survival.pdf",width=7,height=7,onefile=F)
ggsurvplot(gene1, conf.int=F, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("cluster1","cluster2"),legend.title="",  #High Risk,Low Risk是group里面分的类型
           palette=c("#E8B34B", "#5F5CA5"), 
           title="", 
           risk.table.height=0.3)
dev.off()


##mutation
library(maftools)
setwd("D:\\1.Research\\02Lactic\\03mutation")
maf<-read.maf('D:\\1.Research\\02Lactic\\01Data\\TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf')

#col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
#names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
#               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
col = c("Missense_Mutation" = "#2E58A4", "Frame_Shift_Del" = "#9A58A4","Multi_Hit"="#9D9DF1","In_Frame_Ins"="#F39D71","Frame_Shift_Ins"="#459D71","In_Frame_Del"="#B69D71","Nonsense_Mutation"="brown","Splice_Site"="yellow")

genename <- c("ACTN3","HAGH","HIF1A","LDHA","LDHAL6B","LDHB",
"LDHC","LDHD","MIR210","PARK7","PER2","PFKFB2","PNKD","SLC25A12","C12orf5","TP53")
pdf("mutation.pdf",height=7,width=7)
oncoplot(maf = maf, removeNonMutated = T, genes=genename,
         colors = col,writeMatrix = T,
         annotationFontSize = 0.1)
dev.off()
###CNV
setwd("D:\\1.Research\\02Lactic\\04GISTIC")
mydata <- read.table("all_data_by_genes.txt",header=T,sep='\t',row.names=1,check.names=F)

genename <- c("ACTN3","HAGH","HIF1A","LDHA","LDHAL6B","LDHB",
"LDHC","LDHD","MIR210","PARK7","PER2","PFKFB2","PNKD","SLC25A12","C12orf5","TP53")

myname <- t(as.data.frame(strsplit(rownames(mydata),split="[|]")))
data1 <- mydata[which(myname %in% genename),]
data2 <- data1[1:16,c(-1,-2)]

result <- matrix(,16,4)
for (i in 1:16){
  result[i,1] <- length(which(data2[i,]>0.3))
  result[i,2] <- length(which(data2[i,]<(0-0.3)))
}
result[,3] <- result[,1]/(dim(data2)[2])
result[,4] <- result[,2]/(dim(data2)[2])
rownames(result) <- rownames(data2)
colnames(result) <- c("AmpCount","DelCount","AmpProb","DelProb")

library(ggplot2)
setwd("D:\\1.Research\\02Lactic\\04GISTIC")
CNV <- read.table("LUAD_LA_CNV_picture.txt",header=T,sep='\t',check.names=F)

CNV$gene <- factor(CNV$gene,level=c("PFKFB2","SLC25A12","ACTN3","PNKD",
"HAGH","PER2","LDHB","HIF1A","PARK7","C12orf5","LDHC","LDHA","LDHD","MIR210","LDHAL6B","TP53"))
pdf("CNV_boxplot.pdf",width=7,height=7)

ggplot(data=CNV,aes(gene,prob,fill=type))+
geom_bar(stat="identity",position="stack", color="black", width=0.7,size=0.25)+
scale_fill_manual(values = c("#DA4453","#4A89DC")) +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))
dev.off()


###corplot
setwd("D:\\1.Research\\02Lactic\\05CIBERSORT")
mydata <- read.table("TCGA.Kallisto.fullIDs.cibersort.relative.tsv",header=T,sep='\t',check.names=F)
LUAD <- mydata[which(mydata[,2]=="LUAD"),]
LUAD[,1] <- substr(LUAD[,1],1,16)
LUAD[,1] <- gsub("[.]","-",LUAD[,1])

LUAD <- LUAD[!duplicated(LUAD$SampleID),]
rownames(LUAD) <- LUAD[,1]
LUAD <- LUAD[which(LUAD[,25]<0.05),3:24]



LUADfpkm <- read.table("D:\\1.Research\\02Lactic\\01Data\\TCGA-LUADgeneFPKMExp.txt",header=T,sep='\t',row.names=1,check.names=F)
genename <- c("PFKFB2","SLC25A12","ACTN3","PNKD",
"HAGH","PER2","LDHB","HIF1A","PARK7","C12orf5","LDHC","LDHA","LDHD","MIR210","LDHAL6B","TP53")

LAfpkm <- log2(LUADfpkm[genename,]+1)
samplename <- intersect(rownames(LUAD),colnames(LAfpkm))

geneVScellCor <- matrix(,22,16)
geneVScellp <- matrix(,22,16)
for (i in 1:22){
  for (j in 1:16){
    geneVScellCor[i,j] <- cor.test(as.numeric(t(LAfpkm)[samplename,j]),as.numeric(LUAD[samplename,i]))$estimate
    geneVScellp[i,j] <- cor.test(as.numeric(t(LAfpkm)[samplename,j]),as.numeric(LUAD[samplename,i]))$p.value
  }
}
rownames(geneVScellCor) <- colnames(LUAD)
rownames(geneVScellp) <- colnames(LUAD)
colnames(geneVScellCor) <- rownames(LAfpkm) 
colnames(geneVScellp) <- rownames(LAfpkm)

write.table(geneVScellCor,"geneVScellCor.txt",quote=F,sep='\t')
write.table(geneVScellp,"geneVScellp.txt",quote=F,sep='\t')

library(ggcorrplot)
pdf("corplot.pdf",height=7,width=7)
ggcorrplot(geneVScellCor,hc.order = F,  #分等级聚类重排矩阵
           ggtheme = ggplot2::theme_void(base_size = 15), #主题修改
           colors = c("#4A89DC","white","#DA4453"), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
           lab = T,lab_size = 2,    #相关系数文本字体大小
           tl.cex = 12,             #坐标轴字体大小
           p.mat = geneVScellp,         #添加显著性信息
           sig.level = 0.05,        #显著性水平
           pch = 4,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
           pch.cex = 6)
dev.off()


