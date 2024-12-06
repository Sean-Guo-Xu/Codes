library(ggplot2)
library(scales)
library(grid)
all=read.table("D:\\bigslice\\new11608.txt",sep="\t")
withclass = read.table("D:\\bigslice\\化合物\\bgc_withclass.txt",sep="\t",header=F,fill=T)
genome=unique(withclass$V1)
num=c()
for (i in genome){
  num=c(num,nrow(withclass[which(withclass$V1 %in% i),]))
}
bgc=cbind(genome,num)
allbgc=all[-which(all$V1 %in% genome),]
nobgc=cbind(allbgc$V1,rep(0,length(allbgc$V1)))
allbgc=rbind(bgc,nobgc)
class=unique(all$V4)
data=c("class",1)
for (i in class)
{
  n=all[which(all$V4 %in% i),1]
  for (j in n){
    data=rbind(data,c(i,allbgc[which(allbgc[,1] %in% j),2]))
  }
}
data=data[-1,]
data=as.data.frame(data)
data$num=as.numeric(data$num)

ord=read.table("D:\\bigslice\\化合物\\纲，从上到下或从右到左.txt")
ord$V1=rev(ord$V1)
data=data[data$V1 %in% ord$V1,]
data$V1=factor(data$V1, levels = ord$V1)
ggplot(data, aes(x=data$V1, y=data$num,color="white",fill="yellow" ))+scale_x_discrete("")+
  theme_classic() +  geom_boxplot(aes(fill="yellow"))+
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),panel.background=element_rect(fill="dodgerblue4",color="grey50"), panel.grid=element_line(color="grey50",size=2),axis.text.x = element_text(angle = 270, hjust = 0.5, vjust = 0.5,size=12),legend.position = "none")+ylab("BGCs/genome")
ggsave("D:\\bigslice\\bgcgeno箱图.tiff",width = 2,height =6 )
class[-which(class %in% ord$V1)]
