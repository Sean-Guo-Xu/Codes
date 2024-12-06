setwd("D:\\bigslice")
withclass = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
spe=read.table(paste("D:\\bigslice\\gcf-rank\\Genus.txt",sep=""),header=T,sep = "\t")
m=matrix(0,nrow=1006,ncol = 26825)
colnames(m)=c(1:26825)
rownames(m)=spe$Name
for (i in 1:293926){
  y=withclass[i,10]
  x=as.numeric(withclass[i,4])
  x=as.character(x)
  m[y,x]=1
}
t=apply(m,2,sum)
t=sort(t,decreasing = T)
c=m[,names(t)]
c=c[,1:500]
write.table(c,"gcf-属分布矩阵.txt",row.names = T,col.names = T,sep = "\t",quote = F)


withclass = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
genus=withclass[which(withclass$V10 %in% "Aspergillus"),]
spe=read.table(paste("D:\\bigslice\\gcf-rank\\Species.txt",sep=""),header=T,sep = "\t")
spe=spe$Name
spe=spe[which(spe %in% genus$V11)]

withclass=withclass[which(withclass$V11 %in% spe),]
withclass=withclass[order(withclass[,4],decreasing = F),]
m=matrix(0,nrow=length(spe),ncol = nrow(withclass[!duplicated(withclass$V4),]))
colnames(m)=withclass[!duplicated(withclass$V4),4]
rownames(m)=spe
for (i in withclass[!duplicated(withclass$V4),4]){
  y=withclass[i,11]
  x=as.numeric(withclass[i,4])
  x=as.character(x)
  m[y,x]=1
}
write.table(m,"Aspergillus-部分分布矩阵.txt",row.names = T,col.names = T,sep = "\t",quote = F)


withclass = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
spe=read.table(paste("D:\\bigslice\\gcf-rank\\Species.txt",sep=""),header=T,sep = "\t")
need=c()
spe=spe$Name
for (i in spe){
  gcf=withclass[which(withclass$V11 %in% i),4]
  if (length(unique(gcf))>10){
    need=c(need,i)
  }
}
withclass=withclass[which(withclass$V11 %in% need),]
withclass=withclass[order(withclass[,4],decreasing = F),]
m=matrix(0,nrow=2013,ncol = nrow(withclass[!duplicated(withclass$V4),]))
colnames(m)=withclass[!duplicated(withclass$V4),4]
rownames(m)=need
for (i in 1:26){
  y=withclass[i,11]
  x=as.numeric(withclass[i,4])
  x=as.character(x)
  m[y,x]=1
}
write.table(m,"gcf-部分分布矩阵.txt",row.names = T,col.names = T,sep = "\t",quote = F)
