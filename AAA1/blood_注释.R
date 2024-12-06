setwd("E:\\AAA_scRNA")
library(Seurat)
load("E:\\AAA_scRNA\\blood\\A_blood.Rdata")
ano = NULL
for(i in 1:12){
  see = read.csv(paste("E:\\AAA_scRNA\\data\\human_Blood_blood_part",i,".csv",sep = ""))
  ano = rbind(ano,see)
}
unique(pbmc$DF.classifications_0.25_0.18_3150)
pbmc= seurat_object
pbmc$cell = ano$cell_type
pbmc=subset(pbmc,type == "AAABlood")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc <- PercentageFeatureSet( pbmc,  "^MT-",col.name = "percent_mt")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent_mt","percent_hb"), ncol = 4,group.by = "type")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <3500 & percent_mt < 15 & percent_hb<5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:20)                #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
pbmc <- RunUMAP(pbmc,dims = 1:20)
DimPlot(pbmc,group.by = "seurat_clusters",label=T)
library(harmony)
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#耗时1min

#细胞聚类
###resolustion AAA 0.8 NORMAL 0.6#######
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  identity()
DimPlot(pbmc,label = T,group.by = "orig.ident")
DimPlot(pbmc,label = T,group.by = c("cell","seurat_clusters"))
DimPlot(pbmc,label = T)
library(DoubletFinder)
pbmc_list <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
pbmc_stats <- summarizeSweep(pbmc_list, GT = FALSE)
bcmvn <- find.pK(pbmc_stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)    
nExp_poi <- round(0.075*nrow(pbmc@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
memory.limit(size=1000000)
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.18_3150")
pbmc = subset(pbmc,DF.classifications_0.25_0.18_3150 == "Singlet" )
DimPlot(pbmc,label = T,group.by = "cell")

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE)
save(pbmc.markers,file="marker.Rdata")
load("marker.Rdata")

label=cbind(pbmc$cell,pbmc$seurat_clusters)
label[label[,2] %in% c(7),1] = "B cell"
label[label[,2] %in% c(19,11,9,6,5,1,23),1] = "T cell"
label[label[,2] %in% c(15),1] = "Dendritic cell"
label[label[,2] %in% c(14,12,3,10,21,2),1] = "Neutrophil"
label[label[,2] %in% c(4,18,24,8,13,16),1] = "Monocyte"
label[label[,2] %in% c(22),1] = "Unsure"
label[label[,2] %in% c(20),1] = "Plasma cell"
label[label[,2] %in% c(16),1] = "Megakaryocyte"
label[label[,2] %in% c(17),1] = "Mast cell"
pbmc$cell=label[,1]
DimPlot(pbmc,label = T,group.by = "cell") 
save(pbmc,file="Blood_AAA.Rdata")
pbmc = subset(pbmc,cell=="Monocyte")
save(pbmc,file="E:\\AAA_scRNA\\blood\\A_blood.Rdata")

library(SingleR)
library(celldex)

cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref_Hematopoietic, labels = ref_Hematopoietic$label.main)
pbmc$singleR=cellpred$labels
DimPlot(pbmc,group.by = c("singleR","cell"),label = T)
pbmc$cell[pbmc$seurat_clusters %in% c("0","15","1","16","20")] = "CD4 T cell"
pbmc$cell[pbmc$seurat_clusters %in% c("2","5","8","14","4")] = "CD8 T cell"
pbmc$cell[pbmc$seurat_clusters %in% c("3","7","9","18")] = "NK cell"
save(pbmc,file="E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref_Hematopoietic, labels = ref_Hematopoietic$label.main)
