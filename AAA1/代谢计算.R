setwd("E:\\AAA_scRNA")
library(phangorn)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(AUCell)
library(GSEABase)
library(GSVA)
library(Seurat)
load("aaa_hg_mye.Rdata")
memory.limit(size = 10000000)
score = 1:85
pbmc = sc.metabolism.Seurat(pbmc,method = "AUCell", imputation =F, ncores = 8, metabolism.type = "KEGG")
save(pbmc,file= "aaa_hg_mye.Rdata")
for( i in  unique(pbmc$cell)){
  part = subset(pbmc, orig.ident == i)
  part = sc.metabolism.Seurat(obj = part, method = "AUCell", imputation =F, ncores = 8, metabolism.type = "KEGG")
  gc()
  part = part@assays$METABOLISM$score
  score = cbind(score,part)
}
score = score[,-1]
pbmc@assays$METABOLISM = score
save(pbmc,file="8AAA.Rdata")

load("6Normal.Rdata")
con = pbmc@assays$METABOLISM
load("8AAA.Rdata")
aaa =pbmc@assays$METABOLISM
metadata = rownames(aaa)
pvalue = NULL
mean =NULL
for(i in 1:85){
  re=t.test(as.numeric(aaa[i,]),as.numeric(con[i,]))
  m = mean(as.numeric(aaa[i,]))-mean(as.numeric(con[i,]))
      pvalue=c(pvalue,re$p.value)
      mean = c(mean,m)
}
metadata=cbind(metadata,pvalue,mean)
colnames(metadata) = c("Metabolism","Pvalue","Differences")  
metadata= as.data.frame(metadata)
metadata$Pvalue =as.numeric(metadata$Pvalue)
metadata$Differences = as.numeric(metadata$Differences)
metadata = metadata[metadata$Pvalue < 0.05,]
ggplot(data, aes(x=Celltype, y=Expression,fill=sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                               method = "t.test")+
    geom_violin(trim=FALSE) + geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线???
    scale_fill_manual(values = c("#1F77B4FF","#D62728FF"))+ 
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度???15度，并向下调???1(hjust = 1)，字体簇为Times大小???20
          axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属???
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显???(size=1)
          legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属???
                                   size=12),
          legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属???
                                    size=12),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+  #不显示网格线
    xlab("")
ggplot(data, aes(x=Celltype, y=Expression,fill=sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                             method = "t.test")+
  geom_violin(trim=FALSE) + geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线???
  scale_fill_manual(values = c("#1F77B4FF","#D62728FF"))+ 
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度???15度，并向下调???1(hjust = 1)，字体簇为Times大小???20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属???
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显???(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属???
                                 size=12),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属???
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  xlab("")