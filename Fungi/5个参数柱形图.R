library(ggplot2)
library(ggsci)
library(ggthemes)
num=c(350,450,550,650,750)

data=c("Genus",1,"group")
for (i in num){
  number=read.table(paste("D:\\bigslice\\10¸ö-5²ÎÊý\\",i,"genus.txt",sep=""),header=T)
  number=cbind(number,rep(i,10))
  colnames(number)[3]="group"
  data=rbind(data,number)
}
data=data[-1,]
data=as.data.frame(data)
data[,2]=as.numeric(data[,2])
genus=read.table("D:\\bigslice\\gcf-rank\\genus.txt",header=T)
genus = genus[!(genus$Name %in% "Others"),]
genus=genus[1:10,]
col=pal_d3()(10)
names(col) = genus$Name
ggplot(data,aes(x=group,y=GCFNumber,fill=factor(data$Genus,levels = data[1:10,1])))+
  geom_bar(stat = 'identity', position = 'dodge', 
           width = 0.8,color='black')+        
  theme_stata()+  scale_fill_manual(values = col,name="")+ylab("Count of GCFs")+xlab("T")+
  theme(axis.text.x = element_text(angle = 315, hjust = 0.5, vjust = 0.5,size=12))
ggsave("D:\\bigslice\\5barfamily.tiff",width = 8,height = 5)
