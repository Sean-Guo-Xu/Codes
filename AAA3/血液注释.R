library(Seurat)
setwd("E:\\AAA_rupture")
load("Blood.Rdata")
pbmc = subset(pbmc,cell %in% "Myeloid")
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc = IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution =1)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "integrated.cca")

library(SingleR)
library(Seurat)
library(celldex)  # Contains pre-built reference datasets
library(SummarizedExperiment)
ref <- celldex::HumanPrimaryCellAtlasData()
singleR_results <- SingleR(test =GetAssayData(pbmc, slot = "data"), ref = ref, labels = ref$label.main)
pbmc$singleR <- singleR_results$labels
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","singleR"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)

library(ggplot2)
pbmc$cell = "1"
pbmc$cell[pbmc$seurat_clusters %in% c(14)] = "cDC" #CD1C
pbmc$cell[pbmc$seurat_clusters %in% c(16)] = "pDC" #LILRA4
pbmc$cell[pbmc$seurat_clusters %in% c(6,11,8,4,3,7)] = "CD14 Monocyte" #CD14
pbmc$cell[pbmc$seurat_clusters %in% c(13)] = "CD16 CD14 Monocyte" #CD16
pbmc$cell[pbmc$seurat_clusters %in% c(12)] = "C1QA Monocyte" #C1QA
pbmc$cell[pbmc$seurat_clusters %in% c(1)] = "MDSC" #S100A8 S100A9
pbmc$cell[pbmc$seurat_clusters %in% c(2,0)] = "Neutrophil" #FCGR3B, CSF3R
pbmc$cell[pbmc$seurat_clusters %in% c(9)] = "Megakaryocyte" #PPBP PF4
pbmc$cell[pbmc$seurat_clusters %in% c(15)] = "Mast cell" #CPA3
pbmc$cell[pbmc$seurat_clusters %in% c(17)] = "Pro Myeloid" #CD34, KIT
pbmc$cell[pbmc$seurat_clusters %in% c(5,10)] = "NK&T cell" #
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","cell"),reduction = "umap")
DimPlot(subset(pbmc,cell %in% c("NK&T cell","B cell"),invert = T),label = T,group.by = "cell",reduction = "umap")
save(pbmc,file="Blood_Myeloid.Rdata")
load("Blood_Myeloid.Rdata")
pbmc = subset(pbmc,cell != "NK&T cell")
markergene =c("CD14","FCGR3A","C1QA","S100A8","S100A9","CD1C","LILRA4","PPBP","PF4","CPA3","CD34","KIT","FCGR3B","CSF3R")
pbmc$cell = factor(pbmc$cell,levels = c("Neutrophil","Pro Myeloid","Mast cell","Megakaryocyte","pDC","cDC","MDSC","C1QA Monocyte","CD16 CD14 Monocyte", "CD14 Monocyte"))
library(ggplot2)
DotPlot(pbmc,features = markergene,group.by = "cell" )+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low = "#4A90E2",mid = "white",high= "#D95D39")

#########################3########################
library(Seurat)
setwd("E:\\AAA_rupture")
load("Blood.Rdata")
pbmc = subset(pbmc,cell %in% "NK&T cell")
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc = IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution =1)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "integrated.cca")

library(SingleR)
library(Seurat)
library(celldex)  # Contains pre-built reference datasets
library(SummarizedExperiment)
ref <- celldex::HumanPrimaryCellAtlasData()
singleR_results <- SingleR(test =GetAssayData(pbmc, slot = "data"), ref = ref, labels = ref$label.main)
pbmc$singleR <- singleR_results$labels
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","singleR"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)

library(ggplot2)
pbmc$cell="1"
pbmc$cell[pbmc$seurat_clusters %in% c(15)] = "γδT" #TRDV2, TRGV9
pbmc$cell[pbmc$seurat_clusters %in% c(11)] = "Treg" #FOXP3, IL2RA
pbmc$cell[pbmc$seurat_clusters %in% c(2,4,16,7)] = "NK cell" #KLRF1, NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(3)] = "Naive CD8+T" #CD8A CD8B CCR7
pbmc$cell[pbmc$seurat_clusters %in% c(6)] = "Memory CD8+T" #GPR183
pbmc$cell[pbmc$seurat_clusters %in% c(1)] = "Effect CD8+T" #GZMH
pbmc$cell[pbmc$seurat_clusters %in% c(0,17)] = "Naive CD4+T" #CD40LG CCR7 TCF7
pbmc$cell[pbmc$seurat_clusters %in% c(5,18,8)] = "Memory CD4+T" #GPR183 
pbmc$cell[pbmc$seurat_clusters %in% c(12)] = "Th2" #GATA3
pbmc$cell[pbmc$seurat_clusters %in% c(10)] = "Th17" #RORC

pbmc$cell[pbmc$seurat_clusters %in% c(13)] = "Megakaryocyte" 
pbmc$cell[pbmc$seurat_clusters %in% c(19)] = "pDC" 
pbmc$cell[pbmc$seurat_clusters %in% c(9,14)] = "MDSC" 
save(pbmc,file="Blood_NK&T.Rdata")

DimPlot(pbmc,label=T,group.by = c("seurat_clusters","cell"),reduction = "umap")
DimPlot(subset(pbmc,cell %in% c("NK&T cell","B cell"),invert = T),label = T,group.by = "cell",reduction = "umap")

load("Blood_NK&T.Rdata")
pbmc = subset(pbmc,cell %in% c("pDC","MDSC","Megakaryocyte"),invert=T)
markergene =c("CD40LG","CCR7","TCF7","GPR183","GATA3","RORC","TRGV9","TRDV2","FOXP3","IL2RA","CD8A","CD8B","GZMH","KLRF1","NKG7")
pbmc$cell = factor(pbmc$cell,levels =rev( c("Naive CD4+T","Memory CD4+T","Th2","Th17","γδT","Treg","Naive CD8+T","Memory CD8+T","Effect CD8+T","NK cell")))
library(ggplot2)
DotPlot(pbmc,features = markergene,group.by = "cell" )+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low = "#4A90E2",mid = "white",high= "#D95D39")

 
