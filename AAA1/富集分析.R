library("clusterProfiler")
library("org.Hs.eg.db")
library(org.Mm.eg.db)
library("enrichplot")
library("ggplot2")
library(Seurat)
library(stringr)
setwd("E:\\AAA_scRNA\\enrich")
load("E:\\AAA_scRNA\\human_aaa_myeloid.Rdata")
pbmc@active.ident = as.factor(pbmc$cell)
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.25)
for (i in unique(pbmc$cell)) {
part = marker[marker$cluster %in% i,]  
ids = bitr(part$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db" )
go.res<- enrichGO(gene = ids$ENTREZID,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
go.res@result$Description = str_to_title(go.res@result$Description)
go.res <- data.frame(go.res) 
write.table(go.res,file=gsub("/","",paste(i,"human GO.txt")),sep="\t",quote=F,row.names = F)
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:5,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:5,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:5,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图
 go.bar=ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ # 图边距
  scale_fill_manual(values = c("#90d1a4","#fbe18c","#FB9A99"))+xlab("")
ggsave(gsub("/","",paste(i,"human GO.tiff")),go.bar,width = 8,height = 5)

}
load("E:\\AAA_scRNA\\mouse_aaa_myeloid.Rdata")
pbmc@active.ident = as.factor(pbmc$cell)
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.25)
for (i in unique(pbmc$cell)) {
  part = marker[marker$cluster %in% i,]  
  ids = bitr(part$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db" )
  go.res<- enrichGO(gene = ids$ENTREZID,
                    OrgDb = org.Mm.eg.db, 
                    pvalueCutoff =0.05, 
                    qvalueCutoff = 0.05,
                    ont="all",
                    readable =T)
  go.res@result$Description = str_to_title(go.res@result$Description)
  go.res <- data.frame(go.res) 
  write.table(go.res,file=gsub("/","",paste(i,"mouse GO.txt")),sep="\t",quote=F,row.names = F)
  goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:5,]
  goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:5,]
  goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:5,]
  go.df <- rbind(goBP,goCC,goMF)
  # 使画出的GO term的顺序与输入一致
  go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
  # 绘图
  go.bar=ggplot(data = go.df, # 绘图使用的数据
                aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
    geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
    coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
    labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
    theme(axis.title = element_text(size = 13), # 坐标轴标题大小
          axis.text = element_text(size = 11), # 坐标轴标签大小
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
          legend.title = element_text(size = 13), # 图例标题大小
          legend.text = element_text(size = 11), # 图例标签大小
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ # 图边距
    scale_fill_manual(values = c("#90d1a4","#fbe18c","#FB9A99"))+xlab("")
  ggsave(gsub("/","",paste(i,"mosue GO.tiff")),go.bar,width = 8,height = 5)
  
}
