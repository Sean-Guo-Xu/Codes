
#install.packages("pheatmap")

setwd("D:\\????\\70GTEx\\12.heatmap")      #???ù???Ŀ¼
data=read.table("circdiff.txt",sep="\t",header=T,row.names=1,check.names=F)
Rank <- ceiling(rank(data, ties.method="first")/ceiling(nrow(data)*ncol(data)/10))
dimnames1=list(rownames(data),colnames(data))
rt=matrix(as.numeric(as.matrix(Rank)),nrow=nrow(data),dimnames=dimnames1)

library(pheatmap)
Type=c(rep("N",3),rep("T",3))    #?޸Ķ??պʹ???????Ʒ??Ŀ
names(Type)=colnames(rt)
Type=as.data.frame(Type)
anncolor=list(Type=c(N="pink",T="cyan"))
pdf("heatmap.pdf",height=10,width=10)     #????????????
rt <- rt[apply(rt, 1, function(x) sd(x)!=0),]
pheatmap(rt,
         
         annotation=Type, 
         show_colnames=F,
         show_rownames = F,
         color = c(colorRampPalette(colors = c("blue","white"))(length(c(seq(-10,10,by=0.01)))/2),colorRampPalette(colors = c("white","red"))(length(c(seq(-10,10,by=0.01)))/2)),
         legend_breaks=seq(-10,10,2),
         cluster_cols =F,
         fontsize = 10,
         scale="row", 
         border_color=F,
         breaks=c(seq(-10,10,by=0.01)),
         fontsize_row=1,
         annotation_colors=anncolor,
         fontsize_col=3)
dev.off()

