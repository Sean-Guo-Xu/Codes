setwd("E:\\AAA_scRNA")
load("all_AAA.Rdata")
pbmc = subset(pbmc,cell == "Myeloid")
library(Seurat)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
write.csv(t(as.matrix(pbmc@assays$RNA@counts[1:1000,])),file = "fibo_1000.csv")
library(harmony)
pbmc = RunHarmony(pbmc,"sample", plot_convergence = TRUE)#ºÄÊ±1min
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  identity()
DimPlot(pbmc,group.by = c("seurat_clusters","cell"),label = T)
load("aaa_hg_mye.Rdata")
markers = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
Infla = c("CXCL2", "CCL3", "CCL4", "IL1B","CLEC4E","IER3", "NFKBIA", "NR4A2")
Foamy = c("TREM2", "CD9", "GPNMB", "SPP1", "CTSL", "LIPA", "ACP5")
Res = c("LYVE1", "CD163", "SEPP1", "FOLR2", "F13A1", "MRC1")
Monocyte = c("VCAN",  "S100A9", "S100A8","PLAC8","FCN1" ,"APOBEC3A")
cDC1 = c("CLEC9A","THBD", "IRF8", "IDO1")
cDC2 = c("CLEC10A", "FCER1A", "CD1C","CD1E")
marker = c(Infla,Foamy,Res,Monocyte,cDC1,cDC2)
DotPlot(pbmc,features = marker)+RotatedAxis()
save(pbmc,file="aaa_hg_mye.Rdata")
pbmc$cell[pbmc$seurat_clusters %in% c(5)] = "cDC2"
pbmc$cell[pbmc$seurat_clusters %in% c(13)] = "cDC1"
pbmc$cell[pbmc$seurat_clusters %in% c(11,7,10,1)] = "Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(0,3,4,6,15,17)] = "Infla"
pbmc$cell[pbmc$seurat_clusters %in% c(2,18)] = "Resident"
pbmc$cell[pbmc$seurat_clusters %in% c(8,16,14)] = "Foamy"
pbmc$cell[pbmc$seurat_clusters %in% c(9)] = "T cell"
pbmc$cell[pbmc$seurat_clusters %in% c(12)] = "Plasma"
pbmc = subset(pbmc, cell != "Plasma")
pbmc = subset(pbmc, cell != "T cell")
DotPlot(pbmc,features = marker,group.by = "cell")+RotatedAxis()
