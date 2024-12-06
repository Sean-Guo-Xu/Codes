library(Hmisc)
library(ggplot2)
library(cowplot)
library(ggpubr)
setwd("D:\\生信\\target数据分析\\药学")
logfc<-read.table("lhdiff.xls")
diff<-read.table("diffexp.txt",header=T)
risk<-read.table("geneRisk.txt",header=T)
rownames(diff)<-diff[,1]
diff<-diff[,-1]
diff<-log(diff+1,2)
diff<-t(diff)
  c<-diff[,1:length(colnames(diff))]
 df<-cbind(risk$futime,c)
 colnames(df)[1]<-"time"
 res2<-rcorr(as.matrix(df))
 flattenCorrMatrix <- function(cormat, pmat) {
   ut <- upper.tri(cormat)
   data.frame(
     row = rownames(cormat)[row(cormat)[ut]],
     column = rownames(cormat)[col(cormat)[ut]],
     cor  =(cormat)[ut],
     p = pmat[ut]
   )
 }
 res<-flattenCorrMatrix(res2$r,res2$P)
 res<-res[which((res$row %in% "time")),]
 res<-res[res$p<0.05,]
 diff<-diff[,which(colnames(diff) %in% res$column)]
 plist<-list()
 for (i in 1:length(colnames(diff))) {
   
   
   V1<-risk$futime
   V2<-diff[,i]
   df1<-cbind(V1,V2)
   
   df1<-cbind(df1,risk$risk)
   df1<-as.data.frame(df1)
   df1[,1]<-as.numeric(df1[[1]])
   df1[,2]<-as.numeric(df1[[2]])
   p<-ggplot(df1, aes(V1, V2)) + 
     xlab("futime")+ylab(colnames(diff)[i])+
     geom_point(aes(colour=V3))+ 
     geom_smooth(method=lm) + 
     theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
     stat_cor(aes(x =V1, y =V2))+guides(size="none",colour = "none")+scale_colour_gradient(low="pink",high = "cyan")
   plist[[i]]<-p
 }
 plot_grid(plist[[1]],plist[[2]],plist[[3]],plist[[4]],plist[[5]],plist[[6]],plist[[7]],plist[[8]],plist[[9]],plist[[10]],plist[[11]],plist[[12]],plist[[13]],plist[[14]],plist[[15]],plist[[16]],plist[[17]],plist[[18]],plist[[19]])
 resup<-res[res$cor>0,]
 resdown<-res[res$cor<0,]
 logfc<-logfc[which(rownames(logfc) %in% res$column),]
       write.table (logfc, file ="diffcor.xls", sep ="\t", row.names =T, col.names =T, quote =FALSE)
 write.table (resup$column, file ="down.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)
 write.table (resdown$column, file ="up.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)
 