library("limma")
setwd("D:\\����\\target���ݷ���\\��mi-m")
merge = read.table("merge.txt",sep="\t",header=T)
risk = read.table("geneRisk.txt",header=T)
rownames(merge) = merge[,1]
merge = merge[,-1]
merge = merge[,-(1:396)]
colnames(merge) = substr(colnames(merge),1,16)
colnames(merge) = gsub("[.]","-",colnames(merge))
merge = merge[,!duplicated(colnames(merge))]
h = risk[risk$risk %in% "high",]
l = risk[risk$risk %in%"low",]
h = megre[,colnames(merge) %in% h$id]
l = megre[,colnames(merge) %in% l$id]

fdrFilter=0.05                                                    #fdr�ٽ�ֵ
logFCfilter=2                                                     #logFC�ٽ�ֵ
conNum=48                                                       #normal����Ʒ��Ŀ
treatNum=47
exp = cbind(l,h)
grade=c(rep(1,conNum),rep(2,treatNum))
exp=as.matrix(rt)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.0,]

outTab=data.frame()
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=treatGeneMeans-conGeneMeans
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
####################GSEA##################333
outTab$logFC = as.numeric(outTab$logFC)
library(org.Hs.eg.db)

R.utils::setOption("clusterProfiler.download.method",'auto') 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
df_id<-bitr(outTab$gene, #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
marker=outTab
marker=marker[which(marker$gene %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker$logFC
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
KEGG <- gseKEGG(gene_fc, organism = "hsa")

kk_gse <- KEGG
kk_gse_entrez <- KEGG_kk_entrez

###条件筛�? 
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj�?<0.25的通路是显著富集的
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

#选择展现NES前几个通路 
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]
pdf(file=paste(unique(i$cell),"gsea.pdf"),width=10,height=12)
gseaplot2(kk_gse,
          up_gsea$ID,#富集的ID编号
          title = "UP_GSEA_all",#标题
          color = "red",#GSEA线条颜色
          base_size = 20,#基础字体大小
          rel_heights = c(1.5, 0.5, 1),#副图的相对高�?
          subplots = 1:3, #要显示哪些副�? 如subplots=c(1,3) #只要第一和第三个�?
          ES_geom = "line",#enrichment score用线还是用点"dot"
) #显示pvalue等信
gseaplot2(kk_gse,
          down_gsea$ID,#富集的ID编号
          title = "UP_GSEA_all",#标题
          color = "red",#GSEA线条颜色
          base_size = 20,#基础字体大小
          rel_heights = c(1.5, 0.5, 1),#副图的相对高�?
          subplots = 1:3, #要显示哪些副�? 如subplots=c(1,3) #只要第一和第三个�?
          ES_geom = "line",#enrichment score用线还是用点"dot"
) #显示pvalue等信
dev.off()

re = KEGG@result
write.table(re,paste(unique(i$cell),"gsea.txt"),sep = "\t",quote = F,row.names = F)