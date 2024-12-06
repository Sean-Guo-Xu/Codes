library(Seurat)
setwd("E:\\AAA_rupture")
load("Aorta.Rdata")
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
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","cell"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)


load("Aorta_Myeloid.Rdata")
TRM=c("LYVE1","MRC1","FOLR2","F13A1","CD163")
Inflammatorymacrophages=c("IL1A","TLR2","NFKBIA","TNF")
Foamymacrophages=c("TREM2","SPP1","APOC1","APOE","CTSB","FABP5")
Proinf = c("CCL3","CCL4","EGR1","EGR2")
cdc =c("CD1C","CD1E")
myedc = c("S100A12","LYZ","FCN1")
gene = c(TRM,Inflammatorymacrophages,Proinf,Foamymacrophages,cdc,myedc)
pbmc$cell="1"
pbmc$cell[pbmc$seurat_clusters %in% c(11,1)] = "Inf Macrophage" #TRDV2, TRGV9
pbmc$cell[pbmc$seurat_clusters %in% c(2,3,14,9)] = "Foam Macrophage" #FOXP3, IL2RA
pbmc$cell[pbmc$seurat_clusters %in% c(5,6)] = "Resident Macrophage" #KLRF1, NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(7,0,10)] = "Proinf Macrophage" #KLRF1, NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(8)] = "cDC" #CD8A CD8B CCR7
pbmc$cell[pbmc$seurat_clusters %in% c(12)] = "Myeloid DC" #GPR183

pbmc$cell[pbmc$seurat_clusters %in% c(13)] = "B cell" #GZMH
pbmc$cell[pbmc$seurat_clusters %in% c(4)] = "SMC" #CD40LG CCR7 TCF7
pbmc=subset(pbmc,cell %in% c("B cell","SMC"),invert=T)
save(pbmc,file="Aorta_Myeloid.Rdata")
library(ggplot2)
pbmc$cell = factor(pbmc$cell, levels = c("Resident Macrophage","Inf Macrophage","Proinf Macrophage","Foam Macrophage", "cDC", "Myeloid DC"  ))
DotPlot(subset(pbmc,cell %in% c("B cell","SMC"),invert = T),features = gene,group.by = "cell")+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low = "#4A90E2",mid = "white",high= "#D95D39")

library(Seurat)
setwd("E:\\AAA_rupture")
load("Aorta.Rdata")
pbmc = subset(pbmc,cell %in% c("B cell","NK&T cell"))
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc = IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution =1.5)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "integrated.cca")

library(SingleR)
library(Seurat)
library(celldex)  # Contains pre-built reference datasets
library(SummarizedExperiment)
ref <- celldex::HumanPrimaryCellAtlasData()
singleR_results <- SingleR(test =GetAssayData(pbmc, slot = "data"), ref = ref, labels = ref$label.main)
pbmc$singleR <- singleR_results$labels
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","cell"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc$cell[pbmc$seurat_clusters %in% c(16)] = "Treg" #FOXP3, IL2RA
pbmc$cell[pbmc$seurat_clusters %in% c(13, 15)] = "NK cell" #KLRF1, NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(5)] = "CD8+T" #CD8A CD8B
pbmc$cell[pbmc$seurat_clusters %in% c(7,12)] = "Naive CD4+T" #CD4 CCR7 TCF7
pbmc$cell[pbmc$seurat_clusters %in% c(1)] = "Memory T" #S100A4
pbmc$cell[pbmc$seurat_clusters %in% c(6)] = "Th17" #RORC MAF 
pbmc$cell[pbmc$seurat_clusters %in% c(11)] = "SMC" #
pbmc$cell[pbmc$seurat_clusters %in% c(8,10,4,14,17)] = "Memory B cell" #CD79B CD79A CD27
pbmc$cell[pbmc$seurat_clusters %in% c(8,10,4,14,17)] = "Naive B cell" #CD19
save(pbmc,file = "Aorat_TB.Rdata")
load("Aorat_TB.Rdata")
pbmc = subset(pbmc,seurat_clusters %in% c(1,7,12,16,5,6,15,13))
save(pbmc,file="Aorta_T.Rdata")
FeaturePlot(pbmc,features = "IL7R")

library(Seurat)
setwd("E:\\AAA_rupture")
load("Aorta T.Rdata")
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc = IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE,k.weight = 50)
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution =1.5)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "integrated.cca")
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","cell"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc$cell[pbmc$seurat_clusters %in% c(7)] = "Treg" #FOXP3, IL2RA
pbmc$cell[pbmc$seurat_clusters %in% c(11,9)] = "NK cell" #KLRF1 NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(14)] = "Th2" #GATA3, MAF
pbmc$cell[pbmc$seurat_clusters %in% c(1,13)] = "Effect CD8+T" #CD8A CD8B GZMK
pbmc$cell[pbmc$seurat_clusters %in% c(8,10)] = "eMemory CD8+T" #GZMB KLRB1
pbmc$cell[pbmc$seurat_clusters %in% c(0,5,6,4)] = "Naive CD4+T" #CD4 CCR7 TCF7
pbmc$cell[pbmc$seurat_clusters %in% c(3,12,2)] = "Memory CD4+T" #GPR183

save(pbmc,file="Aorta T.Rdata")
library(ggplot2)
DotPlot(pbmc,group.by = "cell",features = c("CCR7","TCF7","SELL","GPR183","CD8A"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low = "#4A90E2",mid = "white",high= "#D95D39")

library(Seurat)
setwd("E:\\AAA_rupture")
load("Aorta SMC.Rdata")
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc = IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution =1.5)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "integrated.cca")
DimPlot(pbmc,label=T,group.by = c("cell"),reduction = "umap")
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc$cell[pbmc$seurat_clusters %in% c(18,9,8,11,13,1)] = "Fibroblast" #FOXP3, IL2RA
pbmc$cell[pbmc$seurat_clusters %in% c(10,2,15,16,20,21)] = "EC" #KLRF1 NKG7
pbmc$cell[pbmc$seurat_clusters %in% c(6,7,12)] = "Synthetic SMC" #FN1 COL1A1
pbmc$cell[pbmc$seurat_clusters %in% c(5,0,4,3,22,14)] = "Contractile SMC" #ACTA2 TAGLN
pbmc$cell[pbmc$seurat_clusters %in% c(16)] = "Lipid related SMC" #APOE CD36
pbmc$cell[pbmc$seurat_clusters %in% c(19)] = "Inflammatory SMC" #CXCL4 LYZ CD69 
save(pbmc,file="Aorta SMC.Rdata")

library(ggplot2)
DotPlot(pbmc,group.by = "cell",features = c("CCR7","TCF7","SELL","GPR183","CD8A"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.ticks.y = element_blank()
)+scale_color_gradient2(low = "#4A90E2",mid = "white",high= "#D95D39")
