library(Seurat)
setwd("D:\\bigslice\\pictures")
load("D:\\bigslice\\GCF cluster\\seurat gcfmatrix.Rdata")
marker = FindAllMarkers(pbmc,
                        only.pos = FALSE,
                        min.pct = 0.25,
                        logfc.threshold = 0.25
)
domain = read.table("hmm.tsv",sep="\t",header = T)
data = merge(marker,domain,by.x="gene",by.y = "id")
pbmc$sample[!(pbmc$ano %in% "0")] = "withoutano"
marker = FindMarkers(pbmc,ident.1 = "withoutano",ident.2 = "known",group.by = "sample")
marker = cbind(rownames(marker),marker)
colnames(marker)[1] = "id"
data = merge(marker,domain,by="id")
data = data[order(data$avg_log2FC,decreasing = T),]
FeaturePlot(pbmc,features = "1610", reduction = 'tsne')
FeaturePlot(pbmc,features = "6830", reduction = 'tsne')

write.table(data,"marker domain.xls",sep = "\t",row.names = F,col.names = T)


library(paletteer)
library(ggplot2)
library(scico)
library(RColorBrewer)
scico_palette_show()
library(ggsci)
col<-colorRampPalette(brewer.pal(8,'Set2'))(36)
pbmc$seurat_clusters=pbmc$seurat_clusters-1

ggsave("all_cluster.tiff",width = 6 ,height = 5)
