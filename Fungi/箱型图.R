library(ggplot2)
setwd("D:\\bigslice\\gcf-rank")
class = c("Class","Order","Family","Genus","Species","Organism")
withclass = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
class1=c("Phylum","Class","Order","Family","Genus","Species")
out = cbind("Name","Variance","Rank")
for (i in 1:6 ){
  
  up = read.table(paste(class[i],".txt",sep=""),header=T,sep="\t")
  t = withclass[!duplicated(withclass[i+5]),i+5]
  outtab = cbind("Name","Variance")
  for (j in 1:length(t)){
    data=withclass[which(withclass[,i+5] %in% t[j]),i+6]
    data=data[!duplicated(data)]
    value=up[which(up$Name %in% data),2]
    if (length(value) != length(data)){
      print(paste(i,j,"wrong"))
    }
    value = as.numeric(value)
    outtab=rbind(outtab,cbind(t[j],as.character(sd(na.omit(value),na.rm = T))))

  }
  colnames(outtab) = outtab[1,]
  outtab=outtab[-1,]
  outtab=cbind(outtab,c(rep(class1[i],length(outtab[,1])))) 
  out= rbind(out,outtab)
}
colnames(out) = out[1,]
out=out[-1,]
out=as.data.frame(out)
out[,2]=as.numeric(out[,2])
out=na.omit(out)
write.table(out,"variance.txt",row.names = F,quote = F,sep = "\t")
out = out[!(out$Rank %in% "Phylum"),]
out$Rank <- factor(out$Rank,levels=c("Class","Order","Family","Genus","Species"))
col=c("#925E9FFF","#FDAF91FF","#42B540FF","#0099B4FF","#00468BFF")
names(col)=c("Class","Order","Family","Genus","Species")
library(ggsci)
pdf("boxplot.pdf",height = 5,width=6)
ggplot(out, aes(x=out$Rank, y=log(out$Variance)) )+scale_x_discrete("")+ geom_violin(trim=FALSE,aes(fill=out$Rank)) +scale_fill_manual(values=col,name="")+ylab("Variance of biosynthetic diversity
")+stat_summary(fun=mean, geom="line",color="red", aes(group=1)) +   geom_boxplot(width=0.2,position=position_dodge(0.9))+theme(axis.text.x = element_text(angle = 315, hjust = 0.5, vjust = 0.5,size=12))+geom_jitter(aes(fill=Rank),width =0.2,shape = 21,size=1)+
theme_classic()
dev.off()
ggsave("boxplot.tiff",height = 5,width=6)
