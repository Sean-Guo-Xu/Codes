setwd("E://AAA_rupture")
library(Seurat)
load("first_batch.Rdata")
pbmc1 = seurat_object
sample1 = pbmc1$orig.ident
load("second_batch.Rdata")
pbmc2 = seurat_object
sample2 = pbmc2$orig.ident
gene  =intersect(rownames(pbmc1),rownames(pbmc2))
pbmc1 = subset(pbmc1, features = gene)
pbmc2 = subset(pbmc2, features = gene)
pbmc1 = CreateSeuratObject(pbmc1@assays$RNA@counts)
pbmc1$orig.ident = sample1
pbmc2 = CreateSeuratObject(pbmc2@assays$RNA@counts)
pbmc2$orig.ident = sample2
pbmc1 = subset(pbmc1, orig.ident %in% c("control3","control1"),invert=T )
pbmc2 = subset(pbmc2, orig.ident %in% c("Normal3","AAA9" ) ,invert=T)
pbmc2$orig.ident[pbmc2$orig.ident %in% "AAA5"] = "AAA9"
pbmc = merge(pbmc1,pbmc2)

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc$all = "all"
VlnPlot(pbmc,features = c("percent.mt","percent_hb","nFeature_RNA","percent_hb"),group.by = "all")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& percent.mt < 25 & percent_hb<5)

pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc  = NormalizeData(pbmc)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce, samples="orig.ident",dbr = 0.01)
pbmc <- as.Seurat(sce)
DimPlot(pbmc,  group.by = "scDblFinder.class")
table(pbmc$scDblFinder.class)
#singlet doublet 
#58548    2674  
pbmc = subset(pbmc,scDblFinder.class == "singlet")

pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:30, reduction = "pca")
pbmc <- FindClusters(pbmc, resolution = 1, cluster.name = "unintegrated_clusters")

pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(pbmc,label=T,group.by = c("seurat_clusters","orig.ident"))




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
DimPlot(pbmc,label = T,group.by =c( "seurat_clusters","orig.ident"),reduction = c("umap"))
save(pbmc,file="Aorta.Rdata")
load("Aorta.Rdata")
marker  =  FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc$cell = "cell"
pbmc$cell[pbmc$seurat_clusters %in% c(2,14)] = "NK&T cell"
pbmc$cell[pbmc$seurat_clusters %in% c(21,16,7)] = "EC"
pbmc$cell[pbmc$seurat_clusters %in% c(18)] = "Mast cell"
pbmc$cell[pbmc$seurat_clusters %in% c(1,12,19)] = "B cell"
pbmc$cell[pbmc$seurat_clusters %in% c(0,3,9,17)] = "SMC"
pbmc$cell[pbmc$seurat_clusters %in% c(4,13)] = "Fibroblast"
pbmc$cell[pbmc$seurat_clusters %in% c(20)] = "Proliferation"
pbmc$cell[pbmc$seurat_clusters %in% c(11,5,6,10,8,15)] = "Myeloid"
pbmc$cell[pbmc$seurat_clusters %in% c(22)] = "Plasmacytoid DC"
DimPlot(pbmc,group.by = "cell",label = T)

pbmc$sample = "AAA"
pbmc$sample[pbmc$orig.ident %in% c("AAA9","AAA2")] ="rAAA"
pbmc$sample[pbmc$orig.ident %in% c("control2","Normal1","Normal2")] ="Normal"
DimPlot(pbmc,split.by = "sample",group.by = "cell",label = T)



load("Aorta_Myeloid.Rdata")
mye = pbmc
load("Aorat_TB.Rdata") 
TB = pbmc
load("Aorta T.Rdata")
TC = pbmc
load("Aorta SMC.Rdata")
smc = pbmc
load("Aorta.Rdata")

pbmc$subcell = pbmc$cell
pbmc$subcell[colnames(pbmc) %in% colnames(mye)] = mye$cell
pbmc$subcell[colnames(pbmc) %in% colnames(TB)] = TB$cell
pbmc$subcell[colnames(pbmc) %in% colnames(TC)] = TC$cell
pbmc$subcell[colnames(pbmc) %in% colnames(smc)] = smc$cell
DimPlot(pbmc,label = TRUE,group.by =c( "seurat_clusters","subcell"),reduction = c("umap"))
pbmc = subset(pbmc,subcell != "Myeloid")
DimPlot(pbmc,label = T,group.by =c( "seurat_clusters","subcell"),reduction = c("umap"))
###去除部分双细胞簇###

load("Aorta_score.Rdata")
score = t(score)
score = score[,colnames(score) %in% colnames(pbmc)]
pbmc@assays$Hallmark = score
save(pbmc,file="Aorta.Rdata")
