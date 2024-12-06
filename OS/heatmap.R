library(pheatmap)
setwd("C:\\Users\\56988\\Documents\\Tencent Files\\569886166\\FileRecv")
drug<-read.table("top20.txt",sep = "\t",header=T)
res<-read.table("SUR.txt",header=T,sep = "\t")
res<-res[which(res$name %in% drug$Name),]
res<-res[!duplicated(res$name),]
for (i in 1:length(res[,1])){
  for (j in 1:length(res[1,])) {
    
  
  res[i,j]<-gsub("NaN",0,res[i,j])}
}
res<-res[order(res$name),]
drug<-drug[order(drug$Name),]
res<-cbind(res,drug$Score)
rownames(res)<-res[,1]
res<-res[,-1]
res<-apply(res, 2, as.numeric)
rownames(res)<-drug$Name
res<-t(res)
ress<-as.data.frame(res[10,])
res<-res[-10,]

ress<-t(ress)
rownames(ress)[1]<-"Comprehensive"
pdf("heatmap1.pdf",height=6,width=20)    

#保存输出结果
#rt <- rt[apply(rt, 1, function(x) sd(x)!=0),]

pheatmap(res,
         
         color = colorRampPalette(c("cyan", "white", "pink"))(50),
         show_rownames = T,
         show_colnames = F,
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         scale="none", 
         border_color=F,
         fontsize_number=12,
         fontsize_row=14,
         fontsize_col=14,
         legend_labels = "Score",
         #gaps_row = 9:10,
         angle_col = 315)
dev.off()
pdf("heatmap2.pdf",height=3,width=20)    

#保存输出结果
#rt <- rt[apply(rt, 1, function(x) sd(x)!=0),]

pheatmap(ress,
         
         color = colorRampPalette(c("cyan", "white", "pink"))(50),
         show_rownames = T,
         show_colnames = T,
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         display_numbers=T,
         legend = F,
         border_color=F,
         
         fontsize_number=12,
         fontsize_row=14,
         fontsize_col=14,
         #gaps_row = 9:10,
         angle_col = 315)
dev.off()
