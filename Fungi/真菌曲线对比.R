library(ggplot2)
setwd("D:\\bigslice\\化合物稀释曲线")
number=20304
class = read.table("D:\\bigslice\\化合物稀释曲线\\真菌聚类.txt",sep="\t",header=F,fill=T)
genome=class[!duplicated(class$V1),]
data=c(seq(1,20301,50),number)
for(j in 1:10){
  t=c(rep("1",number))
  ctt=sample(number)
  for(i in 1:number){
    t[i]=genome[ctt[i],1]
  }
  withclass=class
  time=withclass[which(withclass$V1 %in% t[1]),]
  time=time[!duplicated(time$V2),]
  num=c(nrow(time))
  withclass=withclass[-which(withclass$V2 %in% time$V2),]
  for (i in 1:406)
  {
    time=withclass[which(withclass$V1 %in% t[(2+50*(i-1)):(1+50*i)]),]
    nrow(time[!duplicated(time$V1),])
    time=time[!duplicated(time$V2),]
    num=c(num,nrow(time))
    if (nrow(time)==0)
    {
      withclass=withclass
    }else{
      withclass=withclass[-which(withclass$V2 %in% time$V2),]
    }
  }
  time=withclass[which(withclass$V1 %in% t[20302:number]),]
  time=time[!duplicated(time$V2),]
  num=c(num,nrow(time))
  withclass=withclass[-which(withclass$V2 %in% time$V2),]
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
#閫熷害
data=cbind(c(seq(1,20301,50),number),means)
data[,1]=as.numeric(data[,1])
data[,2]=as.numeric(data[,2])
data=as.data.frame(data)
colnames(data)=c("CompoundsNumber","IncreasingSpeed")
comspeed=data

#鎬诲拰
data=cbind(c(seq(1,20301,50),number),sumgcf)
data[,1]=as.numeric(data[,1])
data[,2]=as.numeric(data[,2])
data=as.data.frame(data)
colnames(data)=c("CompoundsNumber","GCFNumber")
comnum=data

#Bacteria compounds
number=13068
class = read.table("D:\\bigslice\\化合物稀释曲线\\细菌聚类.txt",sep="\t",header=F,fill=T)
genome=class[!duplicated(class$V1),]
data=c(seq(1,13051,50),number)
for(j in 1:10){
  t=c(rep("1",number))
  ctt=sample(number)
  for(i in 1:number){
    t[i]=genome[ctt[i],1]
  }
  withclass=class
  time=withclass[which(withclass$V1 %in% t[1]),]
  time=time[!duplicated(time$V2),]
  num=c(nrow(time))
  withclass=withclass[-which(withclass$V2 %in% time$V2),]
  for (i in 1:261)
  {
    time=withclass[which(withclass$V1 %in% t[(2+50*(i-1)):(1+50*i)]),]
    nrow(time[!duplicated(time$V1),])
    time=time[!duplicated(time$V2),]
    num=c(num,nrow(time))
    if (nrow(time)==0)
    {
      withclass=withclass
    }else{
      withclass=withclass[-which(withclass$V2 %in% time$V2),]
    }
  }
  time=withclass[which(withclass$V1 %in% t[13052:number]),]
  time=time[!duplicated(time$V2),]
  num=c(num,nrow(time))
  withclass=withclass[-which(withclass$V2 %in% time$V2),]
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
data=cbind(c(seq(1,13051,50),number),sumgcf)
data[,1]=as.numeric(data[,1])
data[,2]=as.numeric(data[,2])
data=as.data.frame(data)
colnames(data)=c("CompoundsNumber","GCFNumber")
baccomnum=data


#bgc part
number=20304
bgc=read.table("bgc-gcf.txt")
bgcre=c(rep(0,408))
for(sii in 1:10){
modom=sample(293926,number)
class = bgc[which(bgc$V1 %in% modom),]
genome=class[!duplicated(class$V1),]
data=c(seq(1,20301,50),number)
for(j in 1:10){
  t=c(rep("1",number))
  ctt=sample(number)
  for(i in 1:number){
    t[i]=genome[ctt[i],1]
  }
  withclass=class
  time=withclass[which(withclass$V1 %in% t[1]),]
  time=time[!duplicated(time$V2),]
  num=c(nrow(time))
  withclass=withclass[-which(withclass$V2 %in% time$V2),]
  for (i in 1:406)
  {
    time=withclass[which(withclass$V1 %in% t[(2+50*(i-1)):(1+50*i)]),]
    nrow(time[!duplicated(time$V1),])
    time=time[!duplicated(time$V2),]
    num=c(num,nrow(time))
    if (nrow(time)==0)
    {
      withclass=withclass
    }else{
      withclass=withclass[-which(withclass$V2 %in% time$V2),]
    }
  }
  time=withclass[which(withclass$V1 %in% t[20302:number]),]
  time=time[!duplicated(time$V2),]
  num=c(num,nrow(time))
  withclass=withclass[-which(withclass$V2 %in% time$V2),]
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
bgcre=cbind(bgcre,sumgcf)
}


bgcre=bgcre[,-1]
data=bgcre
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

data=cbind(c(seq(1,20301,50),number),means)
data[,1]=as.numeric(data[,1])
data[,2]=as.numeric(data[,2])
data=as.data.frame(data)
colnames(data)=c("BGCNumber","GCFNumber")
bgcnum=data


ggplot()+geom_line(data=bgcnum,aes(BGCNumber,GCFNumber),color="red",lwd=1)+geom_line(data=comnum,aes(CompoundsNumber,GCFNumber),lwd=1)+geom_line(data=baccomnum,aes(CompoundsNumber,GCFNumber),color="blue",lwd=1)+theme_bw()
