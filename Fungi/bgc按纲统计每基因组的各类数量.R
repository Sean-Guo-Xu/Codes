library(ggplot2)
all=read.table("D:\\bigslice\\new11608.txt",sep="\t")
withclass = read.table("D:\\bigslice\\化合物\\bgc_forheatmap4.txt",sep="\t",header=F,fill=T)
classname=unique(withclass$V3)
name=unique(withclass$V7)
genomename=unique(withclass$V1)
genome=unique(all$V1)
m=matrix(0,nrow = length(genome),ncol=length(classname))
rownames(m)=genome
colnames(m)=classname
for (i in genome)
{

  if (i %in% genomename){
    data=withclass[which(withclass$V1 %in% i),3]
    for(j in data)
    {
      m[i,j]=m[i,j]+1
    }
  }
}
num=c(0,0,0,0,0,0,0,0,0)
for (i in name){
  data=all[which(all$V4 %in% i),1]
  data=unique(data)
  n=length(data)
  data=m[which(rownames(m) %in% data),]
  data=rbind(data,c(0,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0))
  data=apply(data,2,sum)/n
  num=rbind(num,data)
}
num=num[-1,]
rownames(num)=name
colnames(num)=colnames(m)
num=round(num,2)
write.table(num,"D:\\bigslice\\化合物\\基因组化合物分类统计.txt",sep="\t",quote=F)
