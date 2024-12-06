setwd("E:\\AAA_scRNA")
library(Seurat)
load("all_AAA.Rdata")

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc$orig.ident = "all"
VlnPlot(pbmc,features = c("percent.mt","percent_hb","nFeature_RNA"),group.by = "orig.ident")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <6000 & percent.mt < 25 & percent_hb<5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))

library(harmony)
pbmc = RunHarmony(pbmc,"sample", plot_convergence = TRUE)#ºÄÊ±1min
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  identity()

library(DoubletFinder)
pbmc_list <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
pbmc_stats <- summarizeSweep(pbmc_list, GT = FALSE)
bcmvn <- find.pK(pbmc_stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] 
pK_bcmvn = as.character(pK_bcmvn) 
pK_bcmvn = as.numeric(pK_bcmvn)
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)    
nExp_poi <- round(0.075*nrow(pbmc@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
memory.limit(size=1000000)
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.17_2840")
pbmc = subset(pbmc,DF.classifications_0.25_0.17_2840 == "Singlet" )
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.17_2840")
save(pbmc,file = "7AAA.Rdata")

gc()

load("all_AAA.Rdata")
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref, labels = ref$label.main)
pbmc$singleR=cellpred$labels
DimPlot(pbmc,group.by =c("singleR","seurat_clusters"),label=T)
save(pbmc,file = "all_AAA.Rdata")

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,min.pct = 0.25,logfc.threshold = 1)
DimPlot(pbmc,group.by = c("singleR","seurat_clusters"),label = TRUE)
NK = c("KLRD1","NCAM1","NKG7") #7
T = c("CD3D","CD3E","CD3G") #0 17 26
Mast = c("KIT","CPA3","TPSB2") #16
Endothelial = c("PECAM1","PLVAP","PTPRB") #24 23 3
Fibroblast = c("DCN","COL1A1","COL1A2") #12 10 13
SMC = c("MYH11","ACTA2","MYL9") #4 25
Monocyte = c("FCN1", "APOBEC3A", "THBS1") #2 5 8 9 19 22
Macrophage = c("CD163", "CD68", "CD14")
B = c("CD79A","MZB1","MS4A1") #1 6 11 18 20 21 27
Proliferation = c("MKI67","TOP2A","CDC20") #15
#14 unsure
markers = c(T,B,NK,Mast,Endothelial,Fibroblast,SMC,Monocyte,Macrophage,Proliferation)
############ÊÖ¶¯×¢ÊÍ#########

DotPlot(pbmc,features = markers)
pbmc$cell=pbmc$orig.ident
pbmc$cell[pbmc$seurat_clusters %in% c(12)] = "NK cell"
pbmc$cell[pbmc$seurat_clusters %in% c(2)] = "T cell"
pbmc$cell[pbmc$seurat_clusters %in% c(15)] = "Mast cell"
pbmc$cell[pbmc$seurat_clusters %in% c(5,11)] = "Endothelial"
pbmc$cell[pbmc$seurat_clusters %in% c(8,9,13,18)] = "Fibroblast"
pbmc$cell[pbmc$seurat_clusters %in% c(4,10)] = "SMC"
pbmc$cell[pbmc$seurat_clusters %in% c(3,17)]="Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(1,7,19)]="Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c(0,6,20)] = "B cell"
pbmc$cell[pbmc$seurat_clusters %in% c(14)] = "Proliferation"
pbmc$cell[pbmc$seurat_clusters %in% c(16)] = "Unsure"

DimPlot(pbmc,group.by = c("cell"))
save(pbmc,file= "all_AAA.Rdata")

color=c("#b71f48","#F39B7F7F","#fbe18c","#e6f49a","#4fa5b2","#4f98c6","#5371b3","#7E61487F","#8591B47F","#DC00007F")
