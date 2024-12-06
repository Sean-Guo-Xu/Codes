library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
R.utils::setOption("clusterProfiler.download.method",'auto') 
setwd("D:\\AAA\\scRNA")
load("all.Rdata")
pbmc = SplitObject (testAB.integrated, split.by = "cell")
for (i in pbmc){
  pbmc.markers=FindMarkers(i,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
  df_id<-bitr(rownames(pbmc.markers), #??????????df??????????SYMBOL??
              fromType = "SYMBOL",#????????ID????
              toType = "ENTREZID",#????????ID????
              OrgDb = "org.Hs.eg.db")
  marker=cbind(rownames(pbmc.markers),pbmc.markers)
  marker=marker[which(marker$`rownames(pbmc.markers)` %in% df_id$SYMBOL),]
  df_id=df_id[!duplicated(df_id$SYMBOL),]
  colnames(marker)[1]="SYMBOL"
  marker=merge(marker, df_id,by="SYMBOL")
  gene_fc=marker$avg_log2FC
  names(gene_fc)=marker$ENTREZID
  gene_fc=gene_fc[order(gene_fc,decreasing = T)]
  KEGG <- gseKEGG(gene_fc, organism = "hsa")
  
  kk_gse <- KEGG
  kk_gse_entrez <- KEGG_kk_entrez
  
  ###条件筛选 
  #一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
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
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line",#enrichment score用线还是用点"dot"
                     ) #显示pvalue等信
  gseaplot2(kk_gse,
                      down_gsea$ID,#富集的ID编号
                      title = "UP_GSEA_all",#标题
                      color = "red",#GSEA线条颜色
                      base_size = 20,#基础字体大小
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line",#enrichment score用线还是用点"dot"
  ) #显示pvalue等信
  dev.off()
  
  re = KEGG@result
  write.table(re,paste(unique(i$cell),"gsea.txt"),sep = "\t",quote = F,row.names = F)}