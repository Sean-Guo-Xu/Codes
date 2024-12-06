library(Seurat)
setwd("E:\\brain")
load("immune_cancer.Rdata")
options(timeout=10000)
pbmc = im
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1.2, cluster.name = "unintegrated_clusters")
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
library(harmony)
pbmc<- RunHarmony(pbmc,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1.5)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:30,reduction.name = "umap")

DimPlot(pbmc,label=T,group.by = c("seurat_clusters","cell"),reduction = "umap")

immarker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc$cell = as.character(pbmc$seurat_clusters)
pbmc$cell[pbmc$seurat_clusters %in% c(4)] = "CD8+T cell" #	CD3D, CD3E,CD8A CD8B
pbmc$cell[pbmc$seurat_clusters %in% c(13)] = "B cell" #CD19, CD79A
pbmc$cell[pbmc$seurat_clusters %in% c(8)] = "NK cell" #	FGFBP2, NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(17,11)] = "Microglial" #P2RY12、TMEM119
pbmc$cell[pbmc$seurat_clusters %in% c(3,23)] = "Myeloid DC" #"FCN1","VCAN",S100A12
pbmc$cell[pbmc$seurat_clusters %in% c(9)] = "cDC" #CD1C CLEC9A
pbmc$cell[pbmc$seurat_clusters %in% c(24)] = "pDC" #GZMB CLIC3
pbmc$cell[pbmc$seurat_clusters %in% c(15)] = "Inf DC" #CD1D LYZ 
pbmc$cell[pbmc$seurat_clusters %in% c(1,2,16,22,10)] = "M2" #MRC1
pbmc$cell[pbmc$seurat_clusters %in% c(18,19)] = "Cycling cells" #MKI67,TOP2A
pbmc$cell[pbmc$seurat_clusters %in% c(14)] = "Neuronal progenitor cell" #MAP2 SOX4 SOX11
pbmc$cell[pbmc$seurat_clusters %in% c(0,12,6,20)] = "M1" #CD86 IL1B CXCL3
pbmc$cell[pbmc$seurat_clusters %in% c(21)] = "Treg" #FOXP3 	
pbmc$cell[pbmc$seurat_clusters %in% c(5)] = "Naive CD4+T cell" #IL7R CCR7 TCF7
pbmc$cell[pbmc$seurat_clusters %in% c(7)] = "Memory CD4+T cell" #IL7R

DimPlot(pbmc,group.by = "cell",label = T,reduction = "umap")
save(pbmc,file="immune_cancer.Rdata")
load("immmune_marker.Rdata")
genes = c("CD3D","CD3E","IL7R","CCR7","TCF7","CD8A","CD8B","FOXP3","MKI67","TOP2A","CD19","CD79A","FGFBP2","NKG7","P2RY12","TMEM119","S100A8","S100A9","VCAN","CD1C","FCER1A","GZMB","CD1D","CD1E","IL1B","CXCL2", "CLIC3","MRC1","CD86","FCGR3A","MAP2","SOX4","SOX11")
DotPlot(pbmc,features = genes)
DotPlot(pbmc,features = genes,group.by = "cell")
save(pbmc,file="im cell.Rdata")

load("immune_normal.Rdata")
pbmc <- NormalizeData(object = pbmc, normalizmakation.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1.2, cluster.name = "unintegrated_clusters")
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
library(harmony)
pbmc<- RunHarmony(pbmc,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:30,reduction.name = "umap")
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1.2)

DimPlot(pbmc,label = T,reduction = "umap",group.by = c("seurat_clusters","orig.ident"))
load("immune_marker.Rdata")
save(pbmc,file="im cell.Rdata")
cellselect = CellSelector(p)
pbmc$index = colnames(pbmc)
pbmc = subset(pbmc,index %in% cellselect,invert = T)
p=DimPlot(pbmc,label = T)
cellselect = CellSelector(p)
pbmc = subset(pbmc,index %in% cellselect,invert = T)
DimPlot(pbmc,label = T,group.by = c("seurat_clusters","SingleR"),reduction = "umap")
write.csv(marker,"immune_marker.csv")
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc), ref = ref, labels = ref$label.main)
pbmc$SingleR = cellpred$labels
DimPlot(pbmc,label = T,group.by = c("orig.ident","SingleR"))
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
save(pbmc,file="immune_normal.Rdata")

pbmc$cell = "Microglial"
pbmc$cell[pbmc$seurat_clusters %in% c(0,5,8,16,7,2,13,10,11,9) ] = "Resting Microglial" #P2RY12、TMEM119
pbmc$cell[pbmc$seurat_clusters %in% c(20,17,23) ] = "M2" #MRC1
pbmc$cell[pbmc$seurat_clusters %in% c(12) ] = "M1" #CD86
pbmc$cell[pbmc$seurat_clusters %in% c(15,18) ] = "Monocyte" #HLA-DRA FCN1
pbmc$cell[pbmc$seurat_clusters %in% c(24) ] = "T cell" #CD3E,CD3D
pbmc$cell[pbmc$seurat_clusters %in% c(14,6,21,22,4,3)] = "Stem cells" #MAP2 SOX4 SOX11
pbmc$cell[pbmc$seurat_clusters %in% c(19)] = "Astrocyte" #SERPINE2 FABP7 BCAN
pbmc$cell[pbmc$seurat_clusters %in% c(1)] = "Pro-inf Macrophage" #S100A1 CD14 CD68
genes = c("P2RY12","TMEM119","CX3CR1","MRC1","CD86","FCN1","HLA-DRA","CD3D","CD3E","SERPINE2","BCAN","TREM2","SOX4","S100A1","CD14","CD68","SOX11")
DotPlot(pbmc,features = genes)
DimPlot(pbmc,label = T)
save(pbmc,file="immune_normal.Rdata")
write.csv(marker,"immune_marker.csv")

##################胶质#######
load("neuron cell.Rdata")
pbmc = neuron
pbmc <- NormalizeData(object = pbmc, normalizmakation.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.8, cluster.name = "unintegrated_clusters")
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

library(harmony)
pbmc<- RunHarmony(pbmc,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "harmony", reduction.name = "umap")

DimPlot(pbmc,label = T,group.by = c("seurat_clusters","orig.ident"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)

pbmc$cell[pbmc$seurat_clusters %in% c(4) ] = "Neuron" #PCP4
pbmc$cell[pbmc$seurat_clusters %in% c(1,0,10,8,7,5)] = "Oligodendrocyte" #OLIG2 CSPG4 MBP MOG
pbmc$cell[pbmc$seurat_clusters %in% c(2,3,6,9)] = "Astrocyte" #GFAP S100B ALDH1L1http://127.0.0.1:20109/graphics/plot_zoom_png?width=1692&height=676
save(pbmc,file="Neuron cell.Rdata")
DotPlot(pbmc,features = c("CSPG4","OLIG2","MOG","GFAP","PCP4"),group.by = "cell")
write.csv(marker[marker$cluster %in% 6,],"cluster6.csv")


###########marker合并
library(Seurat)
setwd("E:\\brain")
load("immune_cancer.Rdata")
im = pbmc
load("Neuron cell.Rdata")
neuron =pbmc
load("cancer.Rdata")
pbmc$index = colnames(pbmc)
pbmc$cell[colnames(pbmc) %in% colnames(im)] = as.character(im$cell)
pbmc$cell[colnames(pbmc) %in% colnames(neuron)] = as.character(neuron$cell)
pbmc$cell[pbmc$cell %in% "Mono-like DC"] = "Myeloid DC"


save(pbmc,file = "cancer.Rdata")
unique(pbmc$cell)
DimPlot(pbmc,label = T,group.by = c("cell","sample"))

library(ggplot2)
library(ggsci)
color23 = c("#E69F00", "#56B4E9", "#34495E", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#8E44AD", "#1ABC9C", "#F39C12", "#2ECC71", "#E74C3C", "#3498DB", "#9B59B6", "#16A085", "#F1C40F", "#56B4E9", "#E67E22", "#C0392B","#2980B9", "#27AE60")
imcellname=c("Naive CD4+T cell","Memory CD4+T cell","CD8+T cell","Treg","Cycling cells","B cell","NK cell","Microglial","Myeloid DC","cDC","Inf DC","pDC","M2","M1","Neuron", "Oligodendrocyte" , "Astrocyte" ,"G3_MB","G4_MB", "SHH","WNT","Endothelial","Pericytes")
names(color23) = imcellname

load("immune_cancer.Rdata")
pbmc$sample[pbmc$sample %in% "GP4"] ="G4_MB"
pbmc$sample[pbmc$sample %in% "GP3"] ="G3_MB"
pbmc$sample[pbmc$sample %in% "GP3/4"] ="G3_MB"
pbmc = subset(pbmc,cell != "G3_MB")
color = color23[names(color23) %in% unique(pbmc$cell)]
im = pbmc

sample_table <- as.data.frame(table(im$sample,im$cell))
sample_table$Var1 =as.character(sample_table$Var1)
sample_table$Var2 =as.character(sample_table$Var2)
names(sample_table) <- c("Samples","celltype","CellNumber")
sample_table$celltype = factor(sample_table$celltype,levels = names(color))
ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
ggsave("immune/cellproportion.png",height = 7,width = 5)



DimPlot(pbmc,group.by = "cell",cols = color,reduction = "umap")
ggsave("immune/umap.png",height = 5,width = 6.5)
markergene = c("CD3D","CD3E","IL7R","CCR7","TCF7","CD8A","CD8B","FOXP3","MKI67","TOP2A","CD19","CD79A","FGFBP2","NKG7","P2RY12","TMEM119","CD1D","S100A12","FCN1","VCAN","FCER1A","CD1C","ITGAX","FLT3","LYZ","GZMB","CLIC3","CD68","MRC1","CXCL3","IL1B")
DotPlot(pbmc,features = markergene,group.by = "cell")+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low = "#5371b3",mid = "white",high= "#E31A1C")
ggsave("immune/dotplot_im_cancer.png",width =12,height = 6)


pbmc = ScaleData(pbmc)
pbmc@active.ident = pbmc$cell
marker = FindAllMarkers(pbmc, min.pct = 0.2,logfc.threshold = 1)
top10 = NULL
for (i in unique(pbmc$cell)) {
  part = marker[marker$cluster %in% i,]
  part = part[1:10,]
  top10 = rbind(top10,part)
}
library(ggplot2)
DoHeatmap(pbmc,group.by = "cell",features = top10$gene, group.colors =color,size = 4,angle = 90)+scale_fill_gradient2(low = "#5371b3",mid = "#F0E442",high= "#E31A1C")
ggsave("immune/heatmap_im_cancer.png",width = 14, height = 16)
unique(pbmc$cell)
