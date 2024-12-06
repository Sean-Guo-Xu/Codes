library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
setwd("D:\\AAA\\scRNA")
load("M2 seurat.Rdata")
AAA = pbmc
load("5Normal.Rdata")
nor = pbmc
nor = subset(nor,cell =="Macrophage"  )
pbmc=merge(AAA,nor)
pbmc$sample = c(rep("AAA",length(AAA$sample)),rep("Normal",length(nor$sample))) 
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)    
pbmc.markers=FindMarkers(pbmc,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
write.table(pbmc.markers,"diff.txt",sep="\t",quote=F,row.names =T ,col.names = T)
df_id<-bitr(rownames(pbmc.markers), #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
marker=cbind(rownames(pbmc.markers),pbmc.markers)
marker=marker[which(marker$`rownames(pbmc.markers)` %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker$avg_log2FC
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
save(gene_fc,file="genelist.Rdata")
setwd("D:\\AAA\\scRNA")
load("genelist.Rdata")
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
KEGG <- gseKEGG(gene_fc, organism = "hsa")
load()

library("KEGGREST")
path <- keggGet("hsa00600")
gene<- path[[1]]$GENE
gene <- unlist(lapply(gene,function(x) strsplit(x,";")))
gene <- gene[1:length(gene)%%3 == 2]
sp = pbmc.markers[rownames(pbmc.markers) %in% gene,]

data1 =cbind(pbmc$sample,pbmc$sphingolipid,rep("Sphingolipid",length(pbmc$sample)))
data2 =cbind(pbmc$sample,pbmc$apoptosis,rep("Apoptosis",length(pbmc$sample)))
data= rbind(data1,data2)
colnames(data) = c("Sample","AUC","Pathway")
data =as.data.frame(data)
data$AUC = as.numeric(data$AUC)
library(ggplot2)
library(ggpubr)
library(ggsci)
ggplot(data, aes(x=Pathway, y=AUC,fill=Sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                     method = "t.test")+
  geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_npg()+ #设置填充的颜�?
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显�?(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属�?
                                 size=12),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属�?
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("AUC")+xlab("")
ggsave("AUC_AAA_Nor.tiff",width = 4,height = 4)
data = cbind(pbmc$sphingolipid,pbmc$apoptosis,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Sphingolipid","Apoptosis","Sample")
data$Sphingolipid = as.numeric(data$Sphingolipid)
data$Apoptosis = as.numeric(data$Apoptosis)
ggplot(data, aes(x=Sphingolipid, y=Apoptosis, color=Sample)) + scale_color_manual(values = c("brown3","cyan3"))+
  geom_point()+theme_bw()+ stat_smooth(method=lm,formula=y~x,color="black")+stat_cor(color="black",aes(x =Sphingolipid, y =Apoptosis))
ggsave("col_AAA_Nor.tiff",width = 6,height = 5)
