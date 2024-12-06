setwd("D:\\AAA\\omics")
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
marker=read.table("EP300logFC.txt")
df_id<-bitr(marker$V1, #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Mm.eg.db")
marker=marker[which(marker$V1 %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker[,2]
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
GO <- gseGO(gene_fc,ont = "ALL",OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
re = GO@result
View(re)
write.table(re,"EP300_GO_result.txt",sep = "\t",quote=F)

############GO,KEGG#######
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
df_id<-bitr(genes$HGNC.symbol, #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
R.utils::setOption ("clusterProfiler.download.method",'auto')
kegg=enrichKEGG(gene = df_id$ENTREZID, 
           organism = "hsa",keyType = "kegg"
)
go <- enrichGO(gene = df_id$ENTREZID, # Entrez ID б 
               OrgDb = org.Hs.eg.db, # ָ         ݿ 
               keyType = "ENTREZID", # ָ                
               ont = "ALL", #   ѡ,BP(    ѧ    )/CC(ϸ     )/MF(   ӹ   )/ALL(ͬʱָ  )
)
