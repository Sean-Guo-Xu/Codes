library("FGNet")
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(org.Sc.sgd.db)
library(ReactomePA)
library("org.Hs.eg.db")      
library(patchwork)#引用包
rt=read.table("D:\\生信\\target数据分析\\Cmap\\pinkgenes.txt",sep="\t",check.names=F,header=T)    #读取文件
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
geneList=out$logFC
names(geneList)=out$entrezID
gene <- names(geneList)[abs(geneList) > 2]
res1 <- enrichPathway(gene)
res1 <- setReadable(res1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

res1 <- enrichplot::pairwise_termsim(res1)
pdf(file="D:\\生信\\target数据分析\\Cmap\\pinkheatmap.pdf",width = 10,height = 3.8)
heatplot(res1, foldChange=geneList)
dev.off()