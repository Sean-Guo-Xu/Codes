setwd("E:\\AAA_scRNA")
library(Seurat)
load("M_Normal.Rdata")

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
pbmc <- PercentageFeatureSet(pbmc, "^hb[^(P)]", col.name = "percent_hb")
VlnPlot(pbmc,features = c("percent.mt","percent_hb","nFeature_RNA"))
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <5500 & percent.mt < 25)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)


library(harmony)
pbmc = RunHarmony(pbmc,"sample", plot_convergence = TRUE)#ºÄÊ±1min
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.2) %>%
  identity()
  DimPlot(pbmc,label = T,group.by = c("seurat_clusters","orig.ident"))
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
memory.limit(size=10000000)
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.10, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.1_0.01_729")
pbmc = subset(pbmc,DF.classifications_0.1_0.01_729 == "Singlet" )
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.1_0.01_729")
save(pbmc,file = "M_Normal.Rdata")

gc()

load("M_Normal.Rdata")
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref, labels = ref$label.main)
pbmc$singleR=cellpred$labels
DimPlot(pbmc,group.by =c("singleR","seurat_clusters"),label=TRUE)
save(pbmc,file = "all_Normal.Rdata")

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,min.pct = 0.25,logfc.threshold = 1)
DimPlot(pbmc,group.by = c("singleR","seurat_clusters"),label = TRUE)

NK = c("KLRD1","NCAM1","NKG7") #9
T = c("CD3D","CD3E","CD3G") #9
Mast = c("KIT","CPA3","TPSB2") #25 14
Endothelial = c("PECAM1","PLVAP","PTPRB") #24 23 2 10
Fibroblast = c("DCN","COL1A1","COL1A2") #0 8 11 5 20 22
SMC = c("MYH11","ACTA2","MYL9") #6 7 18
Myeloid = c("CD14","CD68","LYZ") #1 3 4 13 15 19
B = c("CD79A","MZB1","MS4A1") #16
Proliferation = c("MKI67","TOP2A","CDC20") #21
Epithelial=c("KRT19", "KRT18","KRT8")
Unsure #12 17
#14 unsure
markers = c(T,B,NK,Mast,Endothelial,Fibroblast,SMC,Myeloid,Proliferation,Epithelial)
############ÊÖ¶¯×¢ÊÍ#########

DotPlot(pbmc,features = markers,group.by = "cell")+RotatedAxis()




pbmc$cell[pbmc$seurat_clusters %in% c(9)] = "T cell"
pbmc$cell[(pbmc$singleR %in% "NK cells" ) & (pbmc$seurat_clusters %in% 9)] = "NK cell"
pbmc$cell[pbmc$seurat_clusters %in% c(25,14)] = "Mast cell"
pbmc$cell[pbmc$seurat_clusters %in% c(24,23,2,10)] = "Endothelial"
pbmc$cell[pbmc$seurat_clusters %in% c(0,8,11,5,20,22)] = "Fibroblast"
pbmc$cell[pbmc$seurat_clusters %in% c(6,7,18)] = "SMC"
pbmc$cell[pbmc$seurat_clusters %in% c(1,3,4,13,15,19)]="Myeloid"
pbmc$cell[pbmc$seurat_clusters %in% c(16)] = "B cell"
pbmc$cell[pbmc$seurat_clusters %in% c(21)] = "Proliferation"
pbmc$cell[pbmc$seurat_clusters %in% c(12,17)] = "Epithelial"

DimPlot(pbmc,group.by = c("cell"))
save(pbmc,file= "all_AAA.Rdata")

color=c("#b71f48","#F39B7F7F","#fbe18c","#e6f49a","#4fa5b2","#4f98c6","#5371b3","#7E61487F","#8591B47F","#DC00007F")
