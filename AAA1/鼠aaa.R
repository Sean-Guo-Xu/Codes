setwd("D:")
library(Seurat)
load("mouse_aaa_all.Rdata")
pbmc=seurat_object
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
pbmc <- PercentageFeatureSet(pbmc, "^Hb[^(P)]", col.name = "percent_hb")
VlnPlot(pbmc,features = c("percent.mt","nFeature_RNA","percent_hb"))
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <6000 & percent.mt < 25 & percent_hb<5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))

library(harmony)
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#??ʱ1min
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  identity()
DimPlot(pbmc,group.by =c("seurat_clusters","orig.ident"),label=T)
library(DoubletFinder)
pbmc_list <- paramSweep(pbmc, PCs = 1:20, sct = FALSE)
pbmc_stats <- summarizeSweep(pbmc_list, GT = FALSE)
bcmvn <- find.pK(pbmc_stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] 
pK_bcmvn = as.character(pK_bcmvn) 
pK_bcmvn = as.numeric(pK_bcmvn)
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)    
nExp_poi <- round(0.075*nrow(pbmc@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
memory.limit(size=1000000)
pbmc <- doubletFinder(pbmc, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.22_1608")
pbmc = subset(pbmc,DF.classifications_0.25_0.22_1608 == "Singlet" )
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.22_1608")
save(pbmc,file = "mouse_aaa_all.Rdata")

gc()

load("M_AAA.Rdata")
library(SingleR)
library(celldex)
ref=celldex::ImmGenData()
cellpred <- SingleR(test =GetAssayData(pbmc, "RNA"), ref = ref, labels = ref$label.main)
pbmc$singleR=cellpred$labels
DimPlot(pbmc,group.by =c("singleR","seurat_clusters"),label=T)
save(pbmc,file = "M_AAA.Rdata")

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,min.pct = 0.25,logfc.threshold = 1)
save(pbmc.markers,file="A_mouse_marker.Rdata")
load("A_mouse_marker.Rdata")
DimPlot(pbmc,group.by = c("seurat_clusters"),label = T)

t <- c("Ccl5", "Tcf7", "Icos", "Il7r")
ec<- c("Cdh5","Pecam1","Tie1")
fb <- c("Pdgfra","Col1a1","Lum")
smc <- c("Acta2","Myh11","Cnn1")
myeloid <- c("Adgre1","Cd68","Cd14")
b <- c("Cd79a","Cd79b","Ly6d")
granulocyte <- c("S100a8","S100a9","Alox5", "Cd53")
Proliferation = c("Mki67","Top2a","Cdc20") 
#14 unsure
markers = c(t,ec,fb,smc ,myeloid,b,granulocyte,Proliferation)
############?ֶ?ע??#########

DotPlot(pbmc,features = markers)+RotatedAxis()
  pbmc$cell = pbmc$orig.ident
pbmc$cell[pbmc$seurat_clusters %in% c(11)] = "T cell"
pbmc$cell[pbmc$seurat_clusters %in% c(20,16,18)] = "Endothelial"
pbmc$cell[pbmc$seurat_clusters %in% c(0,2,3,5,8,9,15)] = "Fibroblast"
pbmc$cell[pbmc$seurat_clusters %in% c(4,19)] = "SMC"
pbmc$cell[pbmc$seurat_clusters %in% c(1,6,7,10,12)]="Myeloid"
pbmc$cell[pbmc$seurat_clusters %in% c(17)] = "B cell"
pbmc$cell[pbmc$seurat_clusters %in% c(14)] = "Granulocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(13)] = "Proliferation"
pbmc=subset(pbmc,seurat_clusters != 21)
pbmc=subset(pbmc,seurat_clusters != 22)

save(pbmc,file="mouse_aaa_all.Rdata")
DimPlot(pbmc,group.by = c("cell"))
save(pbmc,file= "M_AAA.Rdata")

color=c("#b71f48","#F39B7F7F","#fbe18c","#e6f49a","#4fa5b2","#4f98c6","#5371b3","#7E61487F","#8591B47F","#DC00007F")
