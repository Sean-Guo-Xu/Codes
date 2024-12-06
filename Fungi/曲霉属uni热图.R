library(ggplot2)
library(dplyr)
library(scales)
setwd("D:\\bigslice\\gcfscreen")
genus="Saccharomycetes"
class = read.table("D:\\bigslice\\化合物\\bgc_withclass.txt",sep="\t",header=F,fill=T,stringsAsFactors = F)
ano=read.table("D:\\bigslice\\化合物\\GCFannotation.xls",header=T,sep="\t",stringsAsFactors = F)
iano=class[-which(class$V4 %in% ano$gcf_id),]
ano=class[which(class$V4 %in% ano$gcf_id),]

uniclass=class[1,]
for(i in genus){
  class2=class[which(class$V7 %in% i),]
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
#gcf 类型判断
gcf = read.table("D:\\bigslice\\化合物\\GCF class.txt",header=T)
gcfclass=c()
for(i in 1:length(gcf$GCF)){
  
  col=as.numeric(which(gcf[i,3:12] == max(gcf[i,3:12]), arr.ind=TRUE )[2])
  gcfclass=c(gcfclass,colnames(gcf)[col+2])
}
gcf=cbind(gcf,gcfclass)
gcf=gcf[order(gcf$gcfclass),]
unigcf=gcf[which(gcf$GCF %in% uniclass$V4),]
rownames(mat)=unique(unigcf$GCF)
colnames(mat)=allclass$V1
mat=t(mat)
for(i in 1:length(uniclass$V4)){
  r=as.character(uniclass[i,1])
  c=as.character(uniclass[i,4])
  mat[rownames(mat) %in% r, colnames(mat) %in% c]=1
  
}

col=as.character(unigcf$gcfclass)
col=as.data.frame(col)
rownames(col)=unigcf$GCF
row=as.data.frame(allclass$V11 )
col=as.data.frame(col,stringsAsFactors = F)
pdf("酵母纲heatmap.pdf",height=20,width=8)  
pheatmap::pheatmap(mat,annotation_col = col,fontsize_row=1,fontsize_col = 6,show_colnames=T,labels_row=allclass$V11 ,show_rownames = T,cluster_cols = F,cluster_rows = F, color=c("steelblue4","pink"))
dev.off()
#写代码是一项乐趣
gcfclass=read.table("D:\\bigslice\\化合物\\GCF class.txt",header=T,sep="\t")
gcfclass=gcfclass[gcfclass$GCF %in% uni$V4,]
mat = rbind(unigcf$gcfclass,mat)
write.table(mat,"D:\\bigslice\\gcfscreen\\yeast unigcf.txt",col.names = T,row.names = T,sep = "\t",quote = F)


