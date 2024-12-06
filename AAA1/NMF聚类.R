library(NMF)
library(tidyverse)
setwd("E:\\AAA_scRNA")
library(Seurat)
load("all_AAA.Rdata")

pbmc = subset(pbmc,cell == "Myeloid")
meta  = pbmc@meta.data
pbmc <- CreateSeuratObject(pbmc@assays$RNA@counts, project = "pbmc")
pbmc <- NormalizeData(pbmc) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
vm <- pbmc@assays$RNA@scale.data
dims=15
pbmc =RunPCA(pbmc,npcs=dims)
pbmc <- AddMetaData(pbmc,metadata = meta)
## 高变基因表达矩阵的分解
res <- nmf(vm, dims, method = "snmf/r", seed = 'nndsvd') 
pbmc@reductions$nmf <- pbmc@reductions$pca
pbmc@reductions$nmf@cell.embeddings <- t(coef(res))    
pbmc@reductions$nmf@feature.loadings <- basis(res)

library(harmony)
pbmc <- RunHarmony(pbmc,"sample", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(seurat.obj, 'harmony')
harmony.data<-cbind(rownames(harmony_embeddings),harmony_embeddings)
pbmc<- RunUMAP(pbmc, reduction = 'harmony', dims = 1:dims) %>% 
  FindNeighbors(reduction = 'nmf', dims = 1:dims) %>% FindClusters()
   ## 使用nmf的分解结果降维聚类
p=DimPlot(pbmc,label = T)
res@fit
f <- extractFeatures(res, 30L)
f <- lapply(f, function(x) rownames(res)[x])
f <- do.call("rbind", f)
p = CellSelector(p)
pbmc$select = 1
######框选的是不要的
save(res,pbmc,file="nmf_aaa_hm.Rdata")
load("nmf_aaa_hm.Rdata")
pbmc$select[rownames(pbmc@meta.data) %in% p] = "remove"
pbmc$select[!(rownames(pbmc@meta.data) %in% p)]= "select"
pbmc@reductions$nmf = NULL
pbmc= subset(pbmc,select == "select")
markers = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct=0.2)
Infla = c("CXCL2", "CCL3", "CCL4", "IL1B","CLEC4E","IER3", "NFKBIA", "NR4A2")
Foamy = c("TREM2", "CD9", "GPNMB", "SPP1", "CTSL", "LIPA", "ACP5")
Res = c("LYVE1", "CD163", "SEPP1", "FOLR2", "F13A1", "MRC1")
Monocyte = c("VCAN",  "S100A9", "S100A8","PLAC8","FCN1" ,"APOBEC3A")
cDC1 = c("CLEC9A","THBD", "IRF8", "IDO1")
cDC2 = c("CLEC10A", "FCER1A", "CD1C","CD1E")
marker = c(Infla,Foamy,Res,Monocyte,cDC1,cDC2)
DotPlot(pbmc,features = marker)+RotatedAxis()
pbmc$cell[pbmc$seurat_clusters %in% c(5,17)] = "cDC"
pbmc$cell[pbmc$seurat_clusters %in% c(4,18,9)] = "Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(2,3,12,15,16)] = "Infla"
pbmc$cell[pbmc$seurat_clusters %in% c(1,10,11,19,20)] = "Foamy"
pbmc$cell[pbmc$seurat_clusters %in% c(0,6)] = "Res"
DimPlot(pbmc,label = T,group.by = "cell")
