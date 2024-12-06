setwd("E:\\AAA_scRNA")
library(Seurat)
load("8AAA.Rdata")
load("6Normal.Rdata")
pbmc <- con
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc)) 
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
DimPlot(pbmc,group.by = "orig.ident")
library(harmony)
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#耗时1min

#细胞聚类
###resolustion AAA 0.8 NORMAL 0.6#######
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.5) %>%
  identity()
DimPlot(pbmc,label = T)
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter
)
library(DoubletFinder)
pbmc_list <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
pbmc_stats <- summarizeSweep(pbmc_list, GT = FALSE)
bcmvn <- find.pK(pbmc_stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] 
pk_bcmvn = as.character(pK_bcmvn) 
pk_bcmvn = as.numeric(pk_bcmvn)
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)    
nExp_poi <- round(0.075*nrow(pbmc@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
memory.limit(size=1000000)
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = as.numeric(pk_bcmvn), nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.005_2998")
DimPlot(pbmc,label = T,group.by = "cell")
pbmc = subset(pbmc,DF.classifications_0.25_0.005_2998 == "Singlet" )
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref, labels = ref$label.main)
pbmc$singleR=cellpred$labels
pbmc = save(pbmc,file="6Normal.Rdata")
#######3AAA######
label=cbind(pbmc$cell_type,pbmc$seurat_clusters)
label[label[,2] %in% c(6,13),1] = "Endothelial cell"
label[label[,2] %in% c(7),1] = "B cell" 
label[label[,2] %in% c(10,5,22,1,20,18),1]="NK/T cell" 
label[label[,2] %in% c(2,4,12,8,19),1]="Macrophage" 
label[label[,2] %in% c(11,21,15),1]="Fibroblast"
label[label[,2] %in% c(3,14,9),1]= "Smooth muscle cell"
label[label[,2] %in% c(16,17),1]="Unsure" 
pbmc$cell=label[,1]
DimPlot(pbmc,label = T)                  #对细胞分组,优化标准模块化
DimPlot(pbmc,group.by = "cell")
save(pbmc,file="8AAA.Rdata")
#########Normal
label=cbind(pbmc$cell_type,pbmc$seurat_clusters)
label[label[,2] %in% c(26,9,17),1] = "Endothelial cell"
label[label[,2] %in% c(12),1] = "B cell" 
label[label[,2] %in% c(15),1]="NK/T cell" 
label[label[,2] %in% c(22,7,4,13,6,5,3,1,19,14),1]="Macrophage" 
label[label[,2] %in% c(11,10),1]="Fibroblast"
label[label[,2] %in% c(2,24,25,8,27,18,20,23,16),1]= "Smooth muscle cell"
label[label[,2] %in% c(21),1]="Unsure" 
pbmc$cell=label[,1]
DimPlot(pbmc,group.by = "cell")
DimPlot(pbmc,group.by = "orig.ident")
DimPlot(pbmc,group.by = "singleR")
save(pbmc,file="6Normal.Rdata")
###################

