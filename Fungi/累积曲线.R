library(ggplot2)
setwd("D:\\bigslice")
class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
genome=class[!duplicated(class$V1),]
data=c(seq(1,11051,50),11097)
for(j in 1:10){
  t=c(rep("1",11097))
  ctt=sample(11097)
  for(i in 1:11097){
    t[i]=genome[ctt[i],1]
  }
  withclass=class
  time=withclass[which(withclass$V1 %in% t[1]),]
  time=time[!duplicated(time$V4),]
  num=c(nrow(time))
  withclass=withclass[-which(withclass$V4 %in% time$V4),]
  for (i in 1:221)
  {
    time=withclass[which(withclass$V1 %in% t[(2+50*(i-1)):(1+50*i)]),]
    nrow(time[!duplicated(time$V1),])
    time=time[!duplicated(time$V4),]
    num=c(num,nrow(time))
    if (nrow(time)==0)
    {
      withclass=withclass
    }else{
      withclass=withclass[-which(withclass$V4 %in% time$V4),]
    }
  }
  time=withclass[which(withclass$V1 %in% t[11052:11097]),]
  time=time[!duplicated(time$V4),]
  num=c(num,nrow(time))
  withclass=withclass[-which(withclass$V4 %in% time$V4),]
  data=cbind(data,num)
}
data=data[,-1]
means=c(0)
for (i in 1:length(data[,1])) {
  means=c(means,mean(data[i,]))
}
means=means[-1]
sumgcf=c(0)
for (i in 1:length(means)){
  sumgcf=c(sumgcf,sum(means[1:i]))
}
sumgcf=sumgcf[-1]
#速度
data=cbind(c(seq(1,11051,50),11097),means)
data[,1]=as.numeric(data[,1])
data[,2]=as.numeric(data[,2])
data=as.data.frame(data)
colnames(data)=c("GenomesNumber","IncreasingSpeed")
p1=ggplot(data,aes(GenomesNumber,IncreasingSpeed))+geom_point()+geom_smooth(method="loess")

ggsave("speed.tiff",p1)
#总和
data=cbind(c(seq(1,11051,50),11097),sumgcf)
data[,1]=as.numeric(data[,1])
data[,2]=as.numeric(data[,2])
data=as.data.frame(data)
colnames(data)=c("GenomesNumber","GCFNumber")
p2=ggplot(data,aes(GenomesNumber,GCFNumber))+geom_point()+geom_smooth(method="loess")
ggsave("numbers.tiff",p2)