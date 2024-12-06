setwd("E://AAA_rupture")
library(Seurat)
load("Aorta.Rdata")
unique(pbmc$subcell)
pbmc = subset(pbmc,subcell %in% c("SMC","NK&T cell"),invert = T)
pbmc$subcell[pbmc$subcell %in% c("B cell")] = "Memory B cell"
pbmc$cell[pbmc$subcell %in% c("Naive B cell","Memory B cell")] = "B cell"
pbmc$cell[pbmc$subcell %in% c( "Treg", "Naive CD4+T", "Memory CD4+T", "Effect CD8+T", "eMemory CD8+T","Th2" )] = "T cell"
pbmc$cell[pbmc$subcell %in% c( "Inf Macrophage",  "Resident Macrophage", "Proinf Macrophage" , "Foam Macrophage")] = "Macrophage"
pbmc$cell[pbmc$subcell %in% c( "Myeloid DC" ,"cDC","Plasmacytoid DC")] = "DC"
pbmc$cell[pbmc$subcell %in% c( "Fibroblast")] = "Fibroblast"
pbmc$cell[pbmc$subcell %in% c( "EC")] = "EC"
pbmc$cell[pbmc$subcell %in% c( "NK cell")] = "NK cell"
pbmc$cell[pbmc$subcell %in% c( "Contractile SMC","Synthetic SMC","Lipid related SMC","Inflammatory SMC")] = "SMC"
pbmc$sample = factor(pbmc$sample,levels = c("rAAA","AAA","Normal"))
DimPlot(pbmc,label = T,group.by = "cell")

load("Blood.Rdata")
unique(pbmc$subcell)
pbmc$cell[pbmc$subcell %in% c( "Effect CD8+T"  , "Naive CD8+T" , "Th2" ,"Naive CD4+T","γδT","Memory CD8+T","Memory CD4+T","Treg","Th17" )] = "T cell"
pbmc$cell[pbmc$subcell %in% c( "CD14 Monocyte","C1QA Monocyte","CD16 CD14 Monocyte")] = "Monocyte"
pbmc$cell[pbmc$subcell %in% c( "NK cell")] = "NK cell"
pbmc$cell[pbmc$subcell %in% c( "Plasma" )] = "Plasma" 
pbmc$cell[pbmc$subcell %in% c( "Neutrophil")] ="Neutrophil"
pbmc$cell[pbmc$subcell %in% c( "Megakaryocyte")] ="Megakaryocyte"
pbmc$cell[pbmc$subcell %in% c( "Pro Myeloid")] ="Pro Myeloid"
pbmc$cell[pbmc$subcell %in% c(    "Cycling cell" )] =   "Cycling cell" 
pbmc$cell[pbmc$subcell %in% c(    "Mast cell")] = "Mast cell"
pbmc$cell[pbmc$subcell %in% c(    "cDC")] = "cDC"
pbmc$cell[pbmc$subcell %in% c(    "pDC")] = "pDC"
pbmc$cell[pbmc$subcell %in% c(    "MDSC")] = "MDSC"
pbmc$sample = factor(pbmc$sample,levels = c("rAAA","AAA","Normal"))
save(pbmc,file = "Blood.Rdata")

color = c("#D95D39", "#FF8C42", "#4A90E2", "#7B6C95", "#6ABE4A", "#D62828", "#F77F00", "#3D348B", "#FFD60A", "#D9BF77", "#023047", "#8D99AE", "#2EC4B6", "#FFBF69")
names(color) = unique(pbmc$cell)
library(ggplot2)
DimPlot(pbmc,label = F,group.by = "cell",reduction = "umap",split.by = "sample",cols = color)+ggtitle("")+xlab("")+ylab("")+theme(     axis.text = element_blank(), # 文本
                                                                                                                                      axis.ticks = element_blank(),axis.line = element_blank())

ggsave("Aorta_umap.png",width = 14,height = 5) 
library(dplyr)
library(ggalluvial)
Ratio <- pbmc@meta.data %>%
  group_by(sample, cell) %>% # 鍒嗙粍
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
ggplot(Ratio, aes(x =sample, y= relative_freq, fill = cell,
                  stratum=cell, alluvium=cell)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 鍙傛暟knot.pos璁剧疆涓?0.5浣胯繛鎺ヤ负鏇茬嚎闈㈢Н锛屽氨鍍忓父瑙佺殑妗戝熀鍥?
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = color)+theme(axis.line.y = element_blank(),plot.margin = margin(5.5, 10, 5.5, 5.5),axis.ticks.y = element_blank() ,legend.position = "none")+scale_y_continuous(expand = c(0, 0))+ylab("")+xlab("")
ggsave("Aorta_proportion.png",width = 7,height = 3) 

pbmc@active.ident = factor(pbmc$cell,levels = unique(pbmc$cell)[order(unique(pbmc$cell))])
#marker = FindAllMarkers(pbmc, min.pct = 0.2,logfc.threshold = 1)
markergene =c("CD79A","CD19","CD3D","CD3E","KLRK1","NKG7","CD68","CD163","CD1C","FCER1A","KIT","PECAM1","CDH5","PDGFRA","COL1A1","ACTA2","MYH11","MKI67")
markergene =c("CD79A","CD19","CD3D","CD3E","KLRK1","NKG7","CD68","CD163","CD1C","FCER1A","KIT","PECAM1","CDH5","PDGFRA","COL1A1","ACTA2","MYH11","MKI67")
pbmc$newcell = factor(pbmc$cell, levels = rev(c("B cell","T cell","NK cell","Macrophage","DC","Mast cell","EC","Fibroblast","SMC","Proliferation")))
library(ggplot2)
DotPlot(pbmc,features = markergene,group.by = "newcell" )+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low =  "#4A90E2",mid = "white",high= "#D95D39")
ggsave("Aorta_dotplot.png",width = 7,height = 4)
library(pheatmap)
anno = cbind(as.character(pbmc$sample),as.character(pbmc$cell))
rownames(anno) = colnames(pbmc)
score = t(pbmc@assays$Hallmark)
score = score[order(factor(anno[,2],levels =c("B cell","T cell","NK cell","Macrophage","DC","Mast cell","EC","Fibroblast","SMC","Proliferation"))),]
anno = anno[order(factor(anno[,2],levels =c("B cell","T cell","NK cell","Macrophage","DC","Mast cell","EC","Fibroblast","SMC","Proliferation"))),]
score = score[order(factor(anno[,1],levels = c("rAAA","AAA","Normal"))),]
anno = anno[order(factor(anno[,1],levels = c("rAAA","AAA","Normal"))),]
anno = as.data.frame(anno)
annotation_colors <- list(
  V1 = c(rAAA="brown3"  ,AAA = "orange1" ,Normal ="cyan3" ),
  V2 = color[1:length(unique(pbmc$cell))]
)
colnames(score) = gsub("HALLMARK_","",colnames(score))
pdf("Aorta_heatmap.pdf",width = 11,height =6)
pheatmap(t(score),annotation_col = anno,
         annotation_names_col = FALSE,
         cluster_rows = T,cluster_cols = F,show_colnames = F,fontsize_row=6,color = colorRampPalette(c("white", "#D95D39"))(50),annotation_colors = annotation_colors)
dev.off()

score=t(pbmc@assays$Hallmark)
allmean = apply(score,2,mean)
allmean = allmean[order(allmean,decreasing = T)]
padata =NULL
for (j in names(allmean)[1:9]) {
  partdata =cbind(as.character(pbmc$cell),as.character(pbmc$sample),score[,j],rep(j,length(pbmc$cell)))
  padata = rbind(padata,partdata)
}
padata=as.data.frame(padata)
padata$V3=as.numeric(padata$V3)
colnames(padata) = c("cell","sample","score","pathway")
padata$pathway = gsub("HALLMARK_","",padata$pathway)
padata$sample = factor(padata$sample,levels = c("rAAA","AAA","Normal"))
library(ggpubr)
library(ggplot2)
ggviolin(padata,x="sample",y="score",fill = "sample",facet.by = "pathway",add = "mean") +
  stat_compare_means(comparisons = list(c("rAAA", "AAA"),
                                        c("rAAA", "Normal"),
                                        c("AAA", "Normal")), y_position = c(8.0, 8.5, 7.5), 
                     tip_length = c(0.1, 0.05, 0.1,0.05,0.1,0.05), label = "p.signif")+
  scale_fill_manual(values = c("#D95D39", "#FF8C42", "#4A90E2"))+theme(legend.position = "None")+xlab("")+ylab("")
ggsave(paste("Aorta","pathway_vln.png"),height = 12,width = 10)

#################PBMC
load("Blood.Rdata")
color = c("#D95D39", "#FF8C42", "#4A90E2", "#7B6C95", "#6ABE4A", "#D62828", "cyan", "#3D348B", "#FFD60A", "#D9BF77", "brown4", "#8D99AE", "#2EC4B6", "#FFBF69")
names(color) = unique(pbmc$cell)
library(ggplot2)
DimPlot(pbmc,label = F,group.by = "cell",reduction = "umap",split.by = "sample",cols = color)+ggtitle("")+xlab("")+ylab("")+theme(     axis.text = element_blank(), # 文本
                                                                                                                                       axis.ticks = element_blank(),axis.line = element_blank())
ggsave("PBMC_umap.png",width = 14,height = 5) 
library(dplyr)
library(ggalluvial)
Ratio <- pbmc@meta.data %>%
  group_by(sample, cell) %>% # 鍒嗙粍
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
ggplot(Ratio, aes(x =sample, y= relative_freq, fill = cell,color = NA,
                  stratum=cell, alluvium=cell)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 鍙傛暟knot.pos璁剧疆涓?0.5浣胯繛鎺ヤ负鏇茬嚎闈㈢Н锛屽氨鍍忓父瑙佺殑妗戝熀鍥?
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = color)+theme(axis.line.y = element_blank(),plot.margin = margin(5.5, 10, 5.5, 5.5),axis.ticks.y = element_blank() ,legend.position = "none")+scale_y_continuous(expand = c(0, 0))+ylab("")+xlab("")
ggsave("PBMC_proportion.png",width = 7,height = 3) 

pbmc@active.ident = factor(pbmc$cell,levels = unique(pbmc$cell)[order(unique(pbmc$cell))])
marker = FindAllMarkers(pbmc, min.pct = 0.2,logfc.threshold = 1)
markergene =c("CD79A","CD19","CD38","XBP1","CD3D","CD3E","KLRK1","NKG7","CD16","FCN1","CD1C","FCER1A","CLEC4C","PF4","ARG1","S100A8","FCGR3B", "CSF3R","CPA3","ENPP3","KIT","CD34","MKI67")
pbmc$newcell = factor(pbmc$cell, levels = rev(c("B cell","Plasma","T cell","NK cell","Monocyte","cDC","pDC","Megakaryocyte","MDSC","Neutrophil","Mast cell","Pro Myeloid","Cycling cell")))
library(ggplot2)
library(Seurat)
cellrank  =c("B cell","Plasma","T cell","NK cell","Monocyte","cDC","pDC","Megakaryocyte","MDSC","Neutrophil","Mast cell","Pro Myeloid","Cycling cell")
pbmc@assays$Hallmark = NULL
DotPlot(pbmc,features = markergene,group.by = "newcell" )+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low =  "#4A90E2",mid = "white",high= "#D95D39")
ggsave("PBMC_dotplot.png",width = 8,height = 4)
load("Blood.Rdata")
library(pheatmap)
anno = cbind(as.character(pbmc$sample),as.character(pbmc$cell))
rownames(anno) = colnames(pbmc)
score = t(pbmc@assays$Hallmark)
score = score[order(factor(anno[,2],levels =cellrank)),]
anno = anno[order(factor(anno[,2],levels =cellrank)),]
score = score[order(factor(anno[,1],levels = c("rAAA","AAA","Normal"))),]
anno = anno[order(factor(anno[,1],levels = c("rAAA","AAA","Normal"))),]
anno = as.data.frame(anno)
annotation_colors <- list(
  V1 = c(rAAA="brown3"  ,AAA = "orange1" ,Normal ="cyan3" ),
  V2 = color[1:length(unique(pbmc$cell))]
)
colnames(score) = gsub("HALLMARK_","",colnames(score))
pdf("PBMC_heatmap.pdf",width = 11,height =6)
pheatmap(t(score),annotation_col = anno,
         annotation_names_col = FALSE,
         cluster_rows = T,cluster_cols = F,show_colnames = F,fontsize_row=6,color = colorRampPalette(c("white", "#D95D39"))(50),annotation_colors = annotation_colors)
dev.off()

score=t(pbmc@assays$Hallmark)
allmean = apply(score,2,mean)
allmean = allmean[order(allmean,decreasing = T)]
padata =NULL
for (j in names(allmean)[1:9]) {
  partdata =cbind(as.character(pbmc$cell),as.character(pbmc$sample),score[,j],rep(j,length(pbmc$cell)))
  padata = rbind(padata,partdata)
}
padata=as.data.frame(padata)
padata$V3=as.numeric(padata$V3)
colnames(padata) = c("cell","sample","score","pathway")
padata$pathway = gsub("HALLMARK_","",padata$pathway)
padata$sample = factor(padata$sample,levels = c("rAAA","AAA","Normal"))
library(ggpubr)
library(ggplot2)
ggviolin(padata,x="sample",y="score",fill = "sample",facet.by = "pathway",add = "mean") +
  stat_compare_means(comparisons = list(c("rAAA", "AAA"),
                                        c("rAAA", "Normal"),
                                        c("AAA", "Normal")), y_position = c(8.0, 8.5, 7.5), 
                     tip_length = c(0.1, 0.05, 0.1,0.05,0.1,0.05), label = "p.signif")+
  scale_fill_manual(values = c("#D95D39", "#FF8C42", "#4A90E2"))+theme(legend.position = "None")+xlab("")+ylab("")
ggsave(paste("Blood","pathway_vln.png"),height = 12,width = 10)

