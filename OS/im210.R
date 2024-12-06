library(IMvigor210CoreBiologies)
library(plyr)
library(stringr)
library(ggplot2)
library(cowplot)
setwd("C:\\Users\\56988\\Documents\\Tencent Files\\569886166\\FileRecv")
data(cds)
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))
rownames(eset)<-eff_length2$symbol
eset<-log(eset+1,2)
c<-c("HK2","RASSF5","SLC22A3","SNX10","URB2")
coxexp<-eset[which(rownames(eset)%in% c ),]
p<-cds@phenoData
p<-p@data
time<-p[,21:22]
c<-c("futime","fustat")
colnames(time)<-c
coxexp<-t(coxexp)
coxexp<-cbind(time,coxexp)
coxexp$futime<-coxexp$futime/12
write.table (coxexp, file ="coxexp.txt", sep ="\t", row.names =T, col.names =T, quote =FALSE)



#有无反应
riskOut<-read.table("210Risk.txt",header=T)
riskOut<-riskOut[order(riskOut$id),]
p<-p[order(rownames(p)),]

rt<-cbind(riskOut[,length(colnames(riskOut))],p[,2])
rt<-na.omit(rt)
colnames(rt)<-c("risk","BOR_binary")
rt<-as.matrix(rt)

Freq<-c(length(rt[which(rt[,1]==1 & rt[,2]==1) ,][,1]),
length(rt[which(rt[,1]==1 & rt[,2]==2) ,][,1]),
length(rt[which(rt[,1]==2 & rt[,2]==1) ,][,1]),
length(rt[which(rt[,1]==2 & rt[,2]==2) ,][,1]))
Var1<-c("high","high","low","low")
Var2<-c("R","NR","R","NR")
a<-cbind(Var1,Var2,Freq)
a<-as.data.frame(a)
a[,3]<-as.character(a[,3])
a[,3]<-as.numeric(a[,3])
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
c<-a[,3]

pvalue <- chisq.test(c(c,ncol=2))$p.value #卡方检验
library(plyr)
pdf("heatmap1.pdf",height=5,width=3) 
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c("pink","cyan"),label=c("SD/PD","CR/PR"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Risk",y="Percent Weidght",
       fill="")+
 annotate(geom = "text",
           cex=4,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()
write.table (a, file ="h1.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)

#IClevel
rt<-cbind(riskOut[,length(colnames(riskOut))],p[,4])
rt<-na.omit(rt)
colnames(rt)<-c("risk","BOR_binary")
rt<-as.matrix(rt)

Freq<-c(length(rt[which(rt[,1]==1 & rt[,2]==1) ,][,1]),
        length(rt[which(rt[,1]==1 & rt[,2]==2) ,][,1]),
        length(rt[which(rt[,1]==1 & rt[,2]==3) ,][,1]),
        length(rt[which(rt[,1]==2 & rt[,2]==1) ,][,1]),
        length(rt[which(rt[,1]==2 & rt[,2]==2) ,][,1]),
        length(rt[which(rt[,1]==2 & rt[,2]==3) ,][,1]))
Var1<-c("high","high","high","low","low","low")
Var2<-c("IC0","IC1","IC2+","IC0","IC1","IC2+")
a<-cbind(Var1,Var2,Freq)
a<-as.data.frame(a)
a[,3]<-as.character(a[,3])
a[,3]<-as.numeric(a[,3])
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
c<-a[,3]
write.table (a, file ="h2.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)

pvalue <- chisq.test(c(c,ncol=3))$p.value #卡方检验
library(plyr)
pdf("heatmap2.pdf",height=5,width=3) 
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = "stack")+
  scale_fill_manual(values = c("pink","cyan","red"),label=c("IC0","IC1","IC2+"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Risk",y="Percent Weidght",
       fill="")+
  annotate(geom = "text",
           cex=4,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()
#TClevel
rt<-cbind(riskOut[,length(colnames(riskOut))],p[,5])
rt<-na.omit(rt)
colnames(rt)<-c("risk","BOR_binary")
rt<-as.matrix(rt)

Freq<-c(length(rt[which(rt[,1]==1 & rt[,2]==1) ,][,1]),
        length(rt[which(rt[,1]==1 & rt[,2]==2) ,][,1]),
        length(rt[which(rt[,1]==1 & rt[,2]==3) ,][,1]),
        length(rt[which(rt[,1]==2 & rt[,2]==1) ,][,1]),
        length(rt[which(rt[,1]==2 & rt[,2]==2) ,][,1]),
        length(rt[which(rt[,1]==2 & rt[,2]==3) ,][,1]))
Var1<-c("high","high","high","low","low","low")
Var2<-c("TC0","TC1","TC2+","TC0","TC1","TC2+")
a<-cbind(Var1,Var2,Freq)
a<-as.data.frame(a)
a[,3]<-as.character(a[,3])
a[,3]<-as.numeric(a[,3])
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
c<-a[,3]
write.table (a, file ="h3.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)

pvalue <- chisq.test(c(c,ncol=3))$p.value #卡方检验
library(plyr)
pdf("heatmap3.pdf",height=5,width=3)
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = "stack")+
  scale_fill_manual(values = c("pink","cyan","red"),label=c("TC0","TC1","TC2+"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Risk",y="Percent Weidght",
       fill="")+
 annotate(geom = "text",
           cex=4,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()

