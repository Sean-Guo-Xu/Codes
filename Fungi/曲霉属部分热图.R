library(ggplot2)
library(dplyr)
library(scales)
setwd("D:\\bigslice\\gcfscreen")
genus="Aspergillus fumigatus"
class = read.table("D:\\bigslice\\化合物\\bgc_withclass.txt",sep="\t",header=F,fill=T,stringsAsFactors = F)
ano=read.table("D:\\bigslice\\化合物\\GCFannotation.xls",header=T,sep="\t",stringsAsFactors = F)
iano=class[-which(class$V4 %in% ano$gcf_id),]
ano=class[which(class$V4 %in% ano$gcf_id),]

uniclass=class[1,]
for(i in genus){
  class2=class[which(class$V11 %in% i),]
  class1=class[-which(class$V1 %in% class2$V1),]
  gcf=c()
  for(j in unique(class2$V4)){
    if(j %in% class1$V4){
      
    }
    else{
      gcf=c(gcf,j)
    }
  }
  class2=class2[which(class2$V4 %in% gcf),]
  uniclass=rbind(uniclass,class2)
}
uniclass=uniclass[-1,]


all=uniclass

mat=matrix(data=0,nrow = length(unique(all$V4)),ncol=length(unique(all$V1)))

allclass=all[!duplicated(all$V1),]
allclass$V11=as.character(allclass$V11)
allclass=allclass[order(allclass$V11),]

rownames(mat)=unique(all$V4)
colnames(mat)=allclass$V1
mat=t(mat)
for(i in 1:length(uniano$V4)){
  r=as.character(uniano[i,1])
  c=as.character(uniano[i,4])
  mat[rownames(mat) %in% r, colnames(mat) %in% c]=1
  
}
for(i in 1:length(uni$V4)){
  r=as.character(uni[i,1])
  c=as.character(uni[i,4])
  mat[rownames(mat) %in% r, colnames(mat) %in% c]=2
  
}
for(i in 1:length(ano$V4)){
  r=as.character(ano[i,1])
  c=as.character(ano[i,4])
  mat[rownames(mat) %in% r, colnames(mat) %in% c]=3
  
}
for(i in 1:length(class$V4)){
  r=as.character(class[i,1])
  c=as.character(class[i,4])
  mat[which(rownames(mat) %in% r),which( colnames(mat) %in% c)]=4
  
}

pdf("fumigatusheatmap.pdf",height=8,width=15)  
pheatmap::pheatmap(mat,fontsize_row=1,show_colnames=T,show_rownames =F,cluster_cols = F,cluster_rows = F, color=c("steelblue4","pink"))
dev.off()
#写代码是一项乐趣

