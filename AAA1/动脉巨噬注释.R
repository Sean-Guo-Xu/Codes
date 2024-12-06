library(Seurat)
library(ggplot2)
load("E:/AAA_scRNA/8AAA.Rdata")
pbmc = subset(pbmc,cell == "Macrophage")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
DimPlot(pbmc,group.by = "seurat_clusters",label=T)
library(harmony)
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#??ʱ1min

pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.5) %>%
  identity()
DimPlot(pbmc,label = T,group.by = "orig.ident")
DimPlot(pbmc,label = T,group.by = c("orig.ident","seurat_clusters"))
DimPlot(pbmc,label = T)


TRM=c("LYVE1","MRC1","FOLR2","F13A1","CD163")
Inflammatorymacrophages=c("S100A8","IL1A","TNF","TLR2","NFKBIA")
Foamymacrophages=c("TREM2","SPP1","APOC1","APOE","CTSB","FABP5","PLIN2")
Proinf = c("NOS2","IL6","IL1B","CXCL9","CXCL10")
cDC1=c("XCR1","CLEC9A","BATF3","IRF8")
cDC2=c("CLEC10A","FCER1A","CD1C","CD1E")
Monocytes=c("HLA-DRA","CD14","CD64","FCGR3A","CST3","LYZ","FCN1","CD52")
marker = c(TRM,Inflammatorymacrophages,Foamymacrophages,cDC1,cDC2,Monocytes)
DotPlot(pbmc,features = marker)+theme(axis.text.x = element_text(angle = 315))

pbmc$cell = as.character(pbmc$seurat_clusters)
pbmc$cell[pbmc$seurat_clusters %in% c("17","16","1","0","19","21")] = "Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c("18")] = "cDC1"
pbmc$cell[pbmc$seurat_clusters %in% c("11")] = "cDC2"
pbmc$cell[pbmc$seurat_clusters %in% c("5")] = "Resident Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c("12","10","8","7","3")] = "Foamy Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c("9","14","2","13","20","15")] = "Inflammatory Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c("4","6")] = "Intermediate Macrophage"

DimPlot(pbmc,group.by = "cell",label = T)
DotPlot(pbmc,features = marker,group.by = "cell")+theme(axis.text.x = element_text(angle = 315))
save(pbmc,file="E:\\AAA_scRNA\\A_art_Macro.Rdata")
##############################################3
load("E:/AAA_scRNA/6Normal.Rdata")
pbmc = subset(pbmc,cell == "Macrophage")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
DimPlot(pbmc,group.by = "seurat_clusters",label=T)
library(harmony)
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#??ʱ1min

pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.5) %>%
  identity()
DimPlot(pbmc,label = T,group.by = "orig.ident")
DimPlot(pbmc,label = T,group.by = c("orig.ident","seurat_clusters"))
DimPlot(pbmc,label = T)
TRM=c("LYVE1","MRC1","FOLR2","F13A1","CD163")
Inflammatorymacrophages=c("S100A8","IL1A","TNF","TLR2","NFKBIA","AIF1","S100A9","CD74")
Foamymacrophages=c("TREM2","SPP1","APOC1","APOE","CTSB","FABP5","PLIN2")
cDC1=c("XCR1","CLEC9A","BATF3","IRF8")
cDC2=c("CLEC10A","FCER1A","CD1C","CD1E")
Monocytes=c("HLA-DRA","CD14","CD64","FCGR3A","CST3","LYZ","FCN1","CD52")
marker = c(TRM,Inflammatorymacrophages,Foamymacrophages,cDC1,cDC2,Monocytes)
DotPlot(pbmc,features = marker)+theme(axis.text.x = element_text(angle = 270))
pbmc$cell = as.character(pbmc$seurat_clusters)
pbmc$cell[pbmc$seurat_clusters %in% c("16","21")] = "Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c("12")] = "cDC2"
pbmc$cell[pbmc$seurat_clusters %in% c("3","5","8","17","20")] = "Resident Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c("15","14","13","6","10")] = "Foamy Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c("0","1","7")] = "Inflammatory Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c("2","4","18","19","11","9")] = "Intermediate Macrophage"

DimPlot(pbmc,group.by = "cell",label = T)
DotPlot(pbmc,features = marker,group.by = "cell")+theme(axis.text.x = element_text(angle = 315))
save(pbmc,file="E:\\AAA_scRNA\\N_art_Macro.Rdata")
