library(iNEXT)
library(ggplot2)
library(ggsci)
library(ggthemes)
class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
class = class[!(class$V10 %in% "Others") ,]
genus=class[!duplicated(class$V10),10]
curve=list()
listlenght=1
for (j in 1:length(genus)){
 
  genusdata=class[which(class$V10 %in% genus[j]),]
  gcf=c(length(genusdata[!duplicated(genusdata$V1),1]))
  for (i in genusdata[!duplicated(genusdata$V4),4]){
    t=genusdata[which(genusdata$V4 %in% i),]
  t=t[!duplicated(t$V1),]
  gcf=c(gcf,as.numeric(nrow(t)))
  }
  if (length(gcf)>=3){
  curve[[genus[j]]]=gcf
  listlenght=listlenght+1}  

}
curve[[481]]<-NULL

y=iNEXT(curve, datatype="incidence_freq", size=round(seq(1,10000,10)), se=FALSE)

name= names(y$iNextEst)
data=c(0,"method",0,"name")

for (i in 1:993){
  gen=y$iNextEst[[i]]
  gen=gen[,c(1,2,4)]
  gen=cbind(gen,rep(name[i],nrow(gen))) 
  colnames(gen)[4]="name"
  data=rbind(data,gen)
}
data=data[-1,]
data=as.data.frame(data)
data=data[!(data$name %in% "fungal"),]
data$qD=as.numeric(data$qD)
alldata = data[order(data$qD,decreasing = T),]
alldata = alldata[!duplicated(alldata$name),] 
data[,1]=as.numeric(data[,1])
data[,3]=as.numeric(data[,3])
data1=data[which(data$method %in% "interpolated"),]
data2=data[which(data$method %in% "extrapolated"),]
save(data,file="D:\\bigslice\\data\\all_curves_data.Rdata")
genus=read.table("D:\\bigslice\\gcf-rank\\genus.txt",header=T)
genus=genus[genus$Name %in% alldata$name[1:10],]
col=pal_d3()(10)
names(col) = genus$Name
name= names(y$iNextEst)
data=c(0,"method",0,"name")
for (i in genus$Name){
  gen=y$iNextEst[[i]]
  gen=gen[,c(1,2,4)]
  gen=cbind(gen,rep(i,nrow(gen))) 
  colnames(gen)[4]="name"
  data=rbind(data,gen)
}
data=data[-1,]
data=as.data.frame(data)
data[,1]=as.numeric(data[,1])
data[,3]=as.numeric(data[,3])
data3=data[which(data$method %in% "interpolated"),]
data4=data[which(data$method %in% "extrapolated"),]
data3$name = factor(data3$name,levels = alldata$name[1:10])
data4$name = factor(data4$name,levels = alldata$name[1:10])
data3 = data3[order(data3$name),]
data4 = data4[order(data4$name),]
data1=data1[!(data1$name %in% data3$name),]
data2=data2[!(data2$name %in% data3$name),]
pdf("D:\\bigslice\\top10genuscurve.pdf",width = 6,height = 8)
ggplot() +
  geom_line(data=data3,aes(x=data3$t, y=qD,color=data3$name),lwd=1) +
  geom_line(data=data4,aes(x=data4$t, y=data4$qD,color=data4$name),linetype="dashed" , lwd=1)+
  scale_color_d3(name="")+
  geom_line(data=data1,aes(x=data1$t, y=qD,fill=data1$name),lwd=0.8) +
  geom_line(data=data2,aes(x=data2$t, y=data2$qD,fill=data2$name),linetype="dashed" , lwd=0.8)+
  theme_bw()+theme(legend.position = "right")+ylab("Count of GCFs")+xlab("Count of Genomes")+
geom_rect(aes(xmin=0,xmax=10000,ymin=0,ymax=1000),
          fill="green",alpha=0.3)+
  geom_rect(aes(xmin=0,xmax=10000,ymin=1000,ymax=2000),
            fill="yellow",alpha=0.2)+
  geom_rect(aes(xmin=0,xmax=10000,ymin=2000,ymax=Inf),
            fill="red",alpha=0.2)+scale_y_continuous(breaks=c(1000,2000,4000,6000,8000),expand = c(0,0),limits = c(0,8000))+
  scale_x_continuous(expand = c(0,0))
dev.off()
colnames(genus)[1]="name"
out = merge(genus,alldata[1:10,],by="name")
colnames(out)[6] = "predicted"
write.table(out,"D:\\bigslice\\top10data.txt",col.names = T,row.names = F,sep="\t",quote = F)
