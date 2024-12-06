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

class=class[which(class$V7 %in% genus),]
class=class[-which(class$V2 %in% uniclass$V2),]

uniano=uniclass[which(uniclass$V4 %in% ano$V4),]
uni=uniclass
ano=class[which(class$V4 %in% ano$V4),]
class=class[-which(class$V4 %in% ano$V4),]
pie=c("a",1)
pie=rbind(pie,c("uni&ano",length(unique(uniano$V4))))
pie=rbind(pie,c("uni",length(unique(uni$V4))))
pie=rbind(pie,c("ano",length(unique(ano$V4))))
pie=rbind(pie,c("gcf",length(unique(class$V4))))
pie=pie[-1,]
pie=as.data.frame(pie)
pie[,2]=as.character(pie[,2])
pie[,1]=as.character(pie[,1])
pie[,2]=as.numeric(pie[,2])
pie[,2]=pie[,2]/sum(pie[,2])*100
p=ggplot(pie, aes(x="", y=pie$V2, fill=pie$V1))+ theme_classic()+theme(axis.text.x=element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) + scale_fill_brewer(palette="Accent") +
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+scale_fill_manual(values=c(uniano="pink",uni="yellow",gcf="brown3",ano="green"))
ggsave("D:\\bigslice\\pie.tiff",p,height = 5,width = 6)
all=rbind(uniano,uni,ano,class)
mat=matrix(data=0,nrow = length(unique(all$V4)),ncol=length(unique(all$V1)))

allclass=all[!duplicated(all$V1),]
allclass$V11=as.character(allclass$V11)
allclass=allclass[order(allclass$V11),]
#gcf 类型判断
gcf = read.table("D:\\bigslice\\化合物\\GCF class.txt",header=T)
gcfclass=c()
for(i in gcf$GCF){
  
  col=as.numeric(which(gcf[i,3:12] == max(gcf[i,3:12]), arr.ind=TRUE )[2])
  gcfclass=c(gcfclass,colnames(gcf)[col+2])
}
gcf=cbind(gcf,gcfclass)
gcf=gcf[order(gcf$gcfclass),]

unianogcf=gcf[which(gcf$GCF %in% uniano$V4),]
unigcf=gcf[which(gcf$GCF %in% uni$V4),]
anogcf=gcf[which(gcf$GCF %in% ano$V4),]
gcf=gcf[which(gcf$GCF %in% class$V4),]
l=c(length(rownames(unianogcf)),length(rownames(unigcf))+length(rownames(unianogcf)),length(rownames(unigcf))+length(rownames(unianogcf))+length(rownames(anogcf)))
gcf=rbind(unianogcf,unigcf,anogcf,gcf)

rownames(mat)=unique(gcf$GCF)
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
col=as.character(gcf$gcfclass)
col=as.data.frame(col,stringsAsFactors = F)
rownames(col)=gcf$GCF
row=as.data.frame(allclass$V7,stringsAsFactors = F )
pdf("yeastheatmap.pdf",height=25,width=15)  
pheatmap::pheatmap(mat,annotation_col = col,fontsize_row=1,show_colnames=F,labels_row=allclass$V11 ,show_rownames = T,gaps_col = l,cluster_cols = F,cluster_rows = F, color=c("steelblue4","yellow","green","brown3"))
dev.off()
#写代码是一项乐趣
gcfclass=read.table("D:\\bigslice\\化合物\\GCF class.txt",header=T,sep="\t")
gcfclass=gcfclass[gcfclass$GCF %in% uni$V4,]
write.table(gcfclass,"D:\\bigslice\\gcfscreen\\Aspergillus unigcf.txt",col.names = T,row.names = F,sep = "\t",quote = F)


