setwd("C:\\Users\\56988\\Documents\\Tencent Files\\569886166\\FileRecv")
drug<-read.table("D:\\cescdrug.txt")
merge<-read.table("C:\\Users\\56988\\Documents\\Tencent Files\\569886166\\FileRecv\\hl.txt")
drug<-rbind(merge[1,],drug)
write.table (drug, file ="cescdrug.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)
drug<-read.table("cescdrug.txt",header=T)
risk<-read.table("geneRisk.txt",header=T)
drug<-drug[!duplicated(drug$id),]
rownames(drug)<-drug[[1]]
drug<-drug[,-1]
c<-1
h<-risk[which(risk$risk %in% "high"),]
l<-risk[which(risk$risk %in% "low"),]
colnames(drug)<-gsub("[.]","-",colnames(drug))
ich<-drug[,which(colnames(drug) %in% h$id)]
icl<-drug[,which(colnames(drug) %in% l$id)]
for(i in 1:length(colnames(drug))){
  t<-wilcox.test(as.numeric(ich[i,]),as.numeric(icl[i,]))
  print(t$p.value)
  c<-append(c,t$p.value,after = i)
  
}
c<-c[-1]
c<-as.numeric(c)
p<-cbind(rownames(drug),c)
p<-as.data.frame(p)
p[,2]<-as.numeric(p[[2]])
p<-p[p$c<0.01,]


h<-ich
l<-icl
drug<-cbind(h,l)
hm<-rowMeans(h,na.rm = T)
lm<-rowMeans(l,na.rm=T)
for ( i in 1:length(h[,1])) {
  if ((hm[[i]]>lm[[i]]) ){
    drug<-drug[-which(rownames(drug) %in% rownames(h)[i]),]
  }
  
}
drug<-drug[which(rownames(drug) %in% p$V1),]
library(ggsignif)
library(ggplot2)
exp<-""
color<-""
drugname<-""
for (i in 1:length(drug[,1])){
  c1<-drug[i,1:((length(ich[1,])))]
  c2<-drug[i,(length(ich[1,])+1):(length(ich[1,])+length(icl[1,]))]
 e<- c(c1,c2)
 exp<-c(exp,e)
 c<-c(rep("high",length(ich[1,])),rep("low",length(icl[1,])))
 color<-c(color,c)
 d<-c(rep(rownames(drug)[i],length(drug[1,])))
 drugname<-c(drugname,d)
}
all<-cbind(as.numeric(exp),color)
all<-cbind(all,drugname)
all<-as.data.frame(all)
all<-all[-1,]
all[,1]<-as.character(all[[1]])
all[,1]<-as.numeric(all[[1]])
colnames(all)[1]<-"IC50"
colnames(all)[2]<-"Risk"
tiff(file="alldrug.tiff",width = 16, height =8,units ="cm",compression="lzw",bg="white",res=600)
ggplot(all, aes(x=drugname, y=IC50, fill=Risk)) +scale_x_discrete("")+
  theme_classic()+
    geom_boxplot()+theme(axis.text.x = element_text(angle = 315, hjust = 0.5, vjust = 0.5))
dev.off()

for (t in 16:24) {
  i<-rownames(drug)[t]
  
  c1<-rep("high",152)
  c2<-rep("low",152)
  df<-t(drug[which(rownames(drug) %in% i),])
  df1<-df[which(rownames(df) %in% colnames(ich)),]
  df2<-df[which(rownames(df) %in% colnames(icl)),]
  df1<-cbind(df1,c1)
  df2<-cbind(df2,c2)
  ddf<-rbind(df1,df2)
  ddf<-as.data.frame(ddf)
  path<-paste("D:\\",i,".tiff",sep="")
  compared_list = list(c("high", "low"))
  
  IC50<-ddf[[1]]
  Risk<-ddf[[2]]
  ddf<-cbind(IC50,Risk)
  ddf<-as.data.frame(ddf)
  ddf[,1]<-as.numeric(ddf[,1])
  p<-ggplot(ddf ,aes(Risk, IC50))+labs(y=i)+geom_boxplot(aes(fill =Risk)) + geom_signif(comparisons = compared_list, test = wilcox.test, step_increase = 0.2)+theme_bw() 
  ggsave(p,file=path)

}
dev.off()