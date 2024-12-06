setwd("E:\\AAA_scRNA")
library(clusterProfiler)
library(org.Hs.eg.db)
marker=read.table("EP300logFC.txt")
df_id<-bitr(marker$V1, #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
marker=marker[which(marker$V1 %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker[,2]
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
GO <- gseGO(gene_fc,ont = "ALL",OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
re = GO@result
View(re)
write.table(re,"EP300_GO_result.txt",sep = "\t",quote=F)