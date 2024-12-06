library("iNEXT")
library("ggplot2")
setwd("D:\\bigslice\\化合物稀释曲线")
class = read.table("细菌聚类.txt",sep="\t",header=F,fill=T)
name = class[!duplicated(class$V2),2]
data1=c(nrow(class))
for(i in name){
  gcf=class[which(class$V2 == i),]
  gcf=gcf[!duplicated(gcf$V1),]
  data1=c(data1,nrow(gcf))
}
memory.limit(size=80000)
class = read.table("真菌聚类.txt",sep="\t",header=F,fill=T)
name = class[!duplicated(class$V2),2]
data2=c(nrow(class))
for(i in name){
  gcf=class[which(class$V2 == i),]
  gcf=gcf[!duplicated(gcf$V1),]
  data2=c(data2,nrow(gcf))
}


class=read.table("bgc-gcf.txt")
bgc=sample(293925,20304)
bgc=class[bgc,]
name = bgc[!duplicated(bgc$V2),2]
data3=c(20304)
for(i in name){
  gcf=bgc[which(bgc$V2 == i),]
  gcf=gcf[!duplicated(gcf$V1),]
  data3=c(data3,nrow(gcf))
}
y3=iNEXT(data3, datatype="incidence_freq", size=round(seq(1,300000,500)), se=T)
y3=y3$iNextEst
y3=y3[,c(1,2,4)]
for (i in 1:99){
     bgc=sample(293925,20304)
     bgc=class[bgc,]
     name = bgc[!duplicated(bgc$V2),2]
     data3=c(20304)
     for(i in name){
       gcf=bgc[which(bgc$V2 == i),]
       data3=c(data3,nrow(gcf))
     }
     y=iNEXT(data3, datatype="incidence_freq", size=round(seq(1,300000,500)), se=T)
     y=y$iNextEst
     y3=cbind(y3,y$qD)}
write.table(y3,"100times.txt",col.names = F,row.names = F,quote = F,sep = "\t")
test=y3[,3:102]
m=apply(test,1,mean)
s=apply(test, 1, sd)
bgcdata=cbind(y3[,1],m,m+s,m-s)
colnames(bgcdata)=c("t","qD","up","down")


data=list()
data[["Bacteria"]]=data1
data[["Fungi"]]=data2
y=iNEXT(data, datatype="incidence_freq", size=round(seq(1,300000,500)), se=T)

bac=y$iNextEst$Bacteria
fun=y$iNextEst$Fungi
bac=bac[,c(1,4)]
fun=fun[,c(1,4)]
bac=cbind(bac,rep("bac",603))
fun=cbind(fun,rep("fun",603))
colnames(bac)[3]="Name"
colnames(fun)[3]="Name"
bb=cbind(bgcdata[,1:2],rep("bgc",603))
colnames(bb)[3]="Name"
data=rbind(bb,bac[,1:3],fun[,1:3])
data$t=as.numeric(data$t)
data$qD=as.double(data$qD)
data=as.data.frame(data)
data1=data[c(1:43,604:613,1207:1249),]
data2=data[-c(1:43,604:613,1207:1249),]
bgcdata=as.data.frame(bgcdata)
p=ggplot() +
  geom_line(data=data1,aes(x=data1$t/1000, y=qD/1000,color=data1$Name),lwd=1) +  scale_color_manual(values = c("grey","cyan3","brown3"))+
  geom_line(data=data2,aes(x=data2$t/1000, y=data2$qD/1000,color=data2$Name),linetype="dashed" , lwd=1)+geom_ribbon(data=bgcdata,aes(ymin=bgcdata$up/1000, ymax=bgcdata$down/1000, x=bgcdata$t/1000), fill = "cyan3", alpha = 0.3) +
  theme_bw()+xlab("BGC/CompoundNumber")+ylab("ClusterNumber")
ggsave("D:\\bigslice\\com&bgc.tiff",p,width = 6,height = 6)
