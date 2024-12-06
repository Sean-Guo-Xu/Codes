rm(list=ls())
setwd("D:/学习/大创/AAA_AD/转录组学/AAA/GSE7084")

load(file='uni_inter.Rdata')
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegig)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) +
    geom_bar(stat="identity") +
    scale_fill_gradient(low="blue",high="red",guide = FALSE) +
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment")
}

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
# DEG <- DEG[which(DEG$change!='NOT' ,arr.ind = T),]

# df <- bitr(rownames(DEG), fromType = "SYMBOL",
#            toType = c("ENTREZID"),
#            OrgDb = org.Hs.eg.db)

# 用的uni，因为inter基因太少了，腹肌效果很差
DEG <- uni_deg

df <- bitr(rownames(DEG), fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db)

head(df)
head(DEG)
DEG$SYMBOL = rownames(DEG)
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)

gene_up= DEG[DEG$change == 'UP','ENTREZID'] 
gene_down=DEG[DEG$change == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)

## KEGG pathway analysis

###   over-representation test
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    universe     = gene_all,
                    pvalueCutoff = 0.01,
                    qvalueCutoff =0.05)
head(kk.up)[,1:6]
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.05)
head(kk.down)[,1:6]
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]
kk <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    universe     = gene_all,
                    pvalueCutoff = 0.01,
                    qvalueCutoff =0.05)

kegg_diff_dt <- as.data.frame(kk.diff)
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)

ggsave(g_kegg,filename = 'kegg_up_down.png')

save(up_kegg, down_kegg, file = 'kegg_result.Rdata')

# gene.KEGG <- enrichKEGG(gene = df$ENTREZID,
#                         organism = "hsa",
#                         keyType = "kegg",
#                         pvalueCutoff = 0.01,
#                         qvalueCutoff = 0.05 )
# library(enrichplot)
# pic <- barplot(gene.KEGG, showCategory=30)
# print(pic)
# ggsave(pic,filename = 'kegg.png', width = 17, height = 13)

###  GSEA 
# kk_gse <- gseKEGG(geneList     = geneList,
#                   organism     = 'hsa',
#                   nPerm        = 1000,
#                   minGSSize    = 120,
#                   pvalueCutoff = 0.01,
#                   verbose      = FALSE)
# head(kk_gse)[,1:6]
# gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
# 
# down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
# up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
# 
# g_kegg=kegg_plot(up_kegg,down_kegg)
# print(g_kegg)
# ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')


### GO database analysis 

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)


go_enrich_results <- lapply( g_list , function(gene) {
  lapply( c('BP','MF','CC') , function(ont) {
    cat(paste('Now process ',ont ))
    ego <- enrichGO(gene          = gene,
                    universe      = gene_all,
                    OrgDb         = org.Hs.eg.db,
                    ont           = ont ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    print( head(ego) )
    return(ego)
  })
})
save(go_enrich_results,file = 'go_enrich_results.Rdata')

load(file = 'go_enrich_results.Rdata')

# TODO:
# ~~~~~~
setwd("D:/学习/大创/AAA_AD/转录组学/AAA/GSE7084")
load(file = 'uni_inter.Rdata')
write.csv(rownames(inter_deg), file = "intersected_DEGs.txt", row.names = F, quote = F)

uni_up <- rownames(uni_deg[uni_deg$change=="UP",])
uni_down <- rownames(uni_deg[uni_deg$change=="DOWN",])

setwd("D:/学习/大创/AAA_AD/AAA R/富集分析/新_富集分析")
write.csv(uni_up, file = "uni_DEGs_up.txt", row.names = F, quote = F)
write.csv(uni_down, file = "uni_DEGs_down.txt", row.names = F, quote = F)

intersect(kegg$Description,AAA_kegg$X)
