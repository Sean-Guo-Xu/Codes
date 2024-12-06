library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)
library(PupillometryR)
setwd("D:\\bigslice")
withclass = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
class = c("Phylum","Class","Order","Family","Genus","Species")
c1=c("Name")
c2=c("GCFNumber")
c3=("GenomeNumber")
c4=("Class")
standard=c(1)
allbox=cbind(c1,c2,c3,c4)
for (i in 6:11)
{c1=c("Name")
 c2=c("GCFNumber")
 c3=("GenomeNumber")
 c4=("Class")
 
 t = withclass[!duplicated(withclass[i]),i]
 t = t[!(t =="Other")]
  for (j in 1:length(t)){
   name=withclass[withclass[i]==t[j],]
   name1=name[!duplicated(name$V4),]
   name2=name[!duplicated(name$V1),]
   c1<-append(c1,t[j])
   c2<-append(c2,nrow(name1))
   c3<-append(c3,nrow(name2))
   }
  count=cbind(c1,c2,c3)
  colnames(count)=count[1,]
  count = count[-1,]
  count=as.data.frame(count)
  count[,2]=as.numeric(count[,2])
  count[,3]=as.numeric(count[,3])
  count = count[order(count[,2],decreasing = T),]
  write.table(count,paste(class[i-5],".txt",sep = ""),quote = F,row.names = F,sep = "\t")
   
  addboxdata=cbind(count,c(rep(class[i-5],nrow(count))))
  colnames(addboxdata)=c("c1","c2","c3","c4")
  allbox=rbind(allbox,addboxdata)
  if (1 == 2){c1=count[1:10,1] 
  c2=count[1:10,2]
  c3=count[1:10,3]
  rt = rbind(cbind(c1,c2)[-1,],cbind(c1,c3)[-1,])
  rt = cbind(rt,c(rep("GCFNumber",length(c1[-1])),rep("GenomeNumber",length(c1[-1]))))
  rt=as.matrix(rt)
  colnames(rt)=c("Name","Value","Group") 
  rt=as.data.frame(rt)
  rt[,2]=as.numeric(rt[,2])
  rt = rt[order(rt[,2],decreasing = T),]
  rt1=rt[which(rt$Group %in% "GCFNumber"),]
  rt2=rt[which(rt$Group %in% "GenomeNumber"),]
   ggplot(rt, aes(x=rt$Name,y= rt$Value, fill = Group)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.7) +

    scale_fill_manual(values = c("#3C5488B2","#00A087B2")) +
    labs(title = NULL, x = NULL, y = 'Relative abundance (%)', fill = NULL) +
    theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0))}
  
  standard=append(standard,1/2*(mean(addboxdata$c2)+sqrt(mean(addboxdata$c2)^2+var(addboxdata$c2))))
}
colnames(allbox)=allbox[1,]
allbox=allbox[-1,]
allbox = as.data.frame(allbox)
allbox[,2]=as.numeric(allbox[,2])
allbox$Class <- factor(allbox$Class,levels=c("Phylum","Class","Order","Family","Genus","Species"))
ggplot(allbox, aes(x=Class, y=log(GCFNumber,10)) )+scale_x_discrete("")+
  theme_classic()+
  geom_boxplot(aes(fill=Class))+theme(axis.text.x = element_text(angle = 315, hjust = 0.5, vjust = 0.5,size=12))+geom_jitter(aes(fill=Class),width =0.2,shape = 21,size=1)
