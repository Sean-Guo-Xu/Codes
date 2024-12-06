setwd("E:\\AAA_scRNA\\blood")
library(Seurat)
load("A_blood.Rdata")
library(harmony)

pbmc = subset(pbmc ,cell=="Macrophage")

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)

pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#ºÄÊ±1min

pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  identity()
DimPlot(pbmc,label = T,group.by = c("cell","seurat_clusters"))
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
)

pbmc$cell = as.integer(pbmc$seurat_clusters)
pbmc$cell[pbmc$cell %in% c(2,5,12)] = "Monocyte"
pbmc$cell[pbmc$cell %in% c(9,1,15,16)] = "DC"
pbmc$cell[!(pbmc$cell %in% c("DC","Monocyte"))]="Macrophage"
DimPlot(pbmc,group.by = "cell",label = T)
pbmc@active.ident = as.factor(pbmc$cell)
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
)
pbmc$cell[pbmc$seurat_clusters %in% c(7)] = "Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(9)] = "DC"
DimPlot(pbmc,group.by = c("cell","seurat_clusters"))
