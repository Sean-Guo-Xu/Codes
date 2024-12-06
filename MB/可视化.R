color=c("#b71f48","#e45d52","#f57245","#feac60","#fbe18c","#e6f49a","#90d1a4","#4fa5b2","#4f98c6","#5371b3")
COLORTHREE=c("#FDBF6F","#FB9A99" , "#E31A1C")
color = c("#FDBF6F","#e6f49a","#90d1a4","#4fa5b2","#4f98c6","#5371b3", "#E31A1C","#fbe18c")
color14=c("#E41A1C", "#377EB8","#E78AC3", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",  "#A6D854")

library(Seurat)
library(ggplot2)
setwd("E:\\brain")
load("4type.Rdata")
load("meta.Rdata")
score=t(score)
score=as.data.frame(score)
pbmc$`Glycolysis / Gluconeogenesis` = score$`Glycolysis / Gluconeogenesis`
FeaturePlot(pbmc,features = c("Glycolysis / Gluconeogenesis"),cols = c("#5371b3", "#E31A1C"))
ggsave("meta_feature_plot.png",width = 5,height = 5)
unique(pscoreunique(pbmc$cell)
names(color) = unique(pbmc$cell)
ggsave("dimplot_4type.png",width = 6.5,height = 5)
markergene = c("KIT", "LHX5" ,"PAX2","NXPH1","SLC1A3","SOX2","NES","GAP43","GRIK2","MGP","ATOH1",
         "RORA","PCP4","OTX2","RSPO3","EOMES","RFC3")
DotPlot(pbmc,features = markergene,group.by = "cell")+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient(low = "#5371b3",high= "#E31A1C")
ggsave("dotplot_4type.png",width =8,height = 5)
pbmc@active.ident = as.factor(pbmc$cell)
pbmc = ScaleData(pbmc)
marker = FindAllMarkers(pbmc, min.pct = 0.2,logfc.threshold = 1)
top10 = NULL
for (i in unique(pbmc$cell)) {
  part = marker[marker$cluster %in% i,]
  part = part[1:10,]
  top10 = rbind(top10,part)
}
DoHeatmap(pbmc,features = top10$gene, group.colors =color,size = 2,angle = 90)+scale_fill_gradient(low="#6BB7CA",high="#E07B54")
ggsave("heatmap_type4.png",width = 8, height = 12)

load("meta.Rdata")
DoHeatmap(pbmc,features = top10$gene, group.colors =color,size = 2,angle = 90)+scale_fill_gradient(low="#6BB7CA",high="#E07B54")
##############################
##############################
load("immune_cancer.Rdata")
DimPlot(pbmc,group.by = "cell",cols = color14)+xlab("UMAP")+ylab("UMAP")
ggsave("im_cancer_dimplot.png",height = 5,width = 6.5)
markergene = c("CD3D", "CD3E","CD8A", "CD8B","IL7R","CD27","FGFBP2", "NKG7","FOXP3","MKI67","TOP2A","CD19", "CD79A",
               "CD68", "CD163" ,"CD86", "FCGR3A","MRC1","P2RY12","TMEM119", "CX3CR1","CCL3","TNF",
               "FCN1","THBS1","CD1C", "FCER1A","CLEC10A", "LYZ", "CD1D","SOX4", "SOX11")
DotPlot(pbmc,features = markergene,group.by = "cell")+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient(low = "#5371b3",high= "#E31A1C")
ggsave("dotplot_im_cancer.png",width =10,height = 5)
pbmc@active.ident = as.factor(pbmc$cell)
pbmc = ScaleData(pbmc)
marker = FindAllMarkers(pbmc, min.pct = 0.2,logfc.threshold = 1)
top10 = NULL
for (i in unique(pbmc$cell)) {
  part = marker[marker$cluster %in% i,]
  part = part[1:10,]
  top10 = rbind(top10,part)
}
DoHeatmap(pbmc,features = top10$gene, group.colors =color,size = 2,angle = 90)+scale_fill_gradient(low="#6BB7CA",high="#E07B54")
ggsave("heatmap_im_cancer.png",width = 8, height = 12)

##############################
##############################
load("immune_normal.Rdata")
DimPlot(pbmc,group.by = "cell",cols = color14[c(1:8,10)])+xlab("UMAP")+ylab("UMAP")
ggsave("im_normal_dimplot.png",height = 5,width = 6.5)
markergene = c("CD3D", "CD3E",
               "CD68", "CD163" ,"CD86", "FCGR3A","MRC1","P2RY12","TMEM119", "CX3CR1","S100A1", "CD14",
               "FCN1","THBS1","SERPINE2", "FABP7", "BCAN","MAP2","SOX4", "SOX11")
DotPlot(pbmc,features = markergene,group.by = "cell")+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient(low = "#5371b3",high= "#E31A1C")
ggsave("dotplot_normal_cancer.png",width =8,height = 5)
pbmc@active.ident = as.factor(pbmc$cell)
pbmc = ScaleData(pbmc)
marker = FindAllMarkers(pbmc, min.pct = 0.2,logfc.threshold = 1)
top10 = NULL
for (i in unique(pbmc$cell)) {
  part = marker[marker$cluster %in% i,]
  part = part[1:10,]
  top10 = rbind(top10,part)
}
DoHeatmap(pbmc,features = top10$gene, group.colors =color14,size = 2,angle = 90)+scale_fill_gradient(low="#6BB7CA",high="#E07B54")
ggsave("heatmap_normal_cancer.png",width = 8, height = 12)
