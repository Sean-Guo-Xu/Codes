setwd("D:\\bigslice\\GCF cluster")
library(ggplot2)
library(Seurat)
load("seurat gcfmatrix.Rdata")
class=read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
genus=read.table("D:\\bigslice\\gcf-rank\\Genus.txt",header = T,sep="\t")
genus=genus[-6,]
c=c(genus$Name[1:10])
for (i in c) {
  
  mc = class[which(class$V10 %in% i),]
  gcf=unique(mc$V4)
  ml =cbind( pbmc$ano,rownames(pbmc@meta.data))
  ml[,2]=gsub("_f","",ml[,2])
  ml[which((ml[,2] %in% gcf) & (ml[,1] %in% "1")),1]="2"
  pbmc$ml = ml[,1]
  
  p=DimPlot(pbmc,reduction="tsne",group.by = "ml",repel=T,cols =c("0"="white","1"="gray","2"="blue"))+ggtitle (i)+ theme(legend.position = "none") #TSNE??????
  ggsave(paste(i,".tiff"),p,height = 5,width = 5)
}
pbmc=subset(pbmc,deleteano != "0")
p=DimPlot(pbmc,reduction="tsne",group.by = "seurat_clusters",repel=T,label = F)+ggtitle ("Unknown")+ theme(legend.position = "none") #TSNE??????
ggsave(paste("all",".tiff"),p,height = 5,width = 5)
