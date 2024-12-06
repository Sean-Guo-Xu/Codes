library(Seurat)
setwd("E:\\AAA_rupture")
pbmc1 = Read10X("Blood1")
pbmc1 = CreateSeuratObject(pbmc1)
pbmc1$orig.ident = "Blood1"
pbmc2 = Read10X("Blood2")
pbmc2 = CreateSeuratObject(pbmc2)
pbmc2$orig.ident = "Blood2"
pbmc3 = Read10X("Blood3")
pbmc3 = CreateSeuratObject(pbmc3)
pbmc3$orig.ident = "Blood3"
pbmc4 = Read10X("Blood4")
pbmc4 = CreateSeuratObject(pbmc4)
pbmc4$orig.ident = "Blood4"
pbmc6 = Read10X("Blood6")
pbmc6 = CreateSeuratObject(pbmc6)
pbmc6$orig.ident = "Blood6"
gene2 = rownames(pbmc6)
pbmc7 = Read10X("Blood7")
pbmc7 = CreateSeuratObject(pbmc7)
pbmc7$orig.ident = "Blood7"
pbmc = merge(pbmc1,pbmc2)
pbmc = merge(pbmc,pbmc3)
pbmc = merge(pbmc,pbmc4)
pbmc = merge(pbmc,pbmc6)
pbmc = merge(pbmc,pbmc7)

gene1 = rownames(pbmc1)
gene2 = rownames(pbmc2)
gene3 = rownames(pbmc3)
geneset1 = union(gene1,gene2)
geneset1 = union(geneset1,gene3)

gene4 = rownames(pbmc4)
gene6 = rownames(pbmc6)
gene7 = rownames(pbmc7)
geneset2 = union(gene4,gene6)
geneset2 = union(geneset1,gene7)



pbmc1=read.csv("H2_expression_counts.csv",row.names = 1)
pbmc1 = CreateSeuratObject(pbmc1)
pbmc1$orig.ident = "H2"
pbmc2=read.csv("H3_expression_counts.csv",row.names = 1)
pbmc2 = CreateSeuratObject(pbmc2)
pbmc2$orig.ident = "H3"
pbmc3=read.csv("H4_expression_counts.csv",row.names = 1)
pbmc3 = CreateSeuratObject(pbmc3)
pbmc3$orig.ident = "H4"
pbmc4=read.csv("H5_expression_counts.csv",row.names = 1)
pbmc4 = CreateSeuratObject(pbmc4)
pbmc4$orig.ident = "H5"
pbmc5=read.csv("H6_expression_counts.csv",row.names = 1)
pbmc5 = CreateSeuratObject(pbmc5)
pbmc5$orig.ident = "H6"
pbmc6=read.csv("H7_expression_counts.csv",row.names = 1)
pbmc6 = CreateSeuratObject(pbmc6)
pbmc6$orig.ident = "H7"
pbmc = merge(pbmc,pbmc1)
pbmc = merge(pbmc,pbmc2)
pbmc = merge(pbmc,pbmc3)
pbmc = merge(pbmc,pbmc4)
pbmc = merge(pbmc,pbmc5)
pbmc = merge(pbmc,pbmc6)

gene1 = rownames(pbmc1)
gene2 = rownames(pbmc2)
gene3 = rownames(pbmc3)
gene4 = rownames(pbmc4)
gene5 = rownames(pbmc5)
gene6 = rownames(pbmc6)
geneset3 = union(gene1,gene2)
geneset3 = union(geneset3,gene3)
geneset3 = union(geneset3,gene4)
geneset3 = union(geneset3,gene5)
geneset3 = union(geneset3,gene6)
gene = intersect(geneset1,geneset2)
gene = intersect(gene,geneset3)

pbmc = subset(pbmc,features = gene)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc$all = "all"
VlnPlot(pbmc,features = c("percent.mt","percent_hb","nFeature_RNA","percent_hb"),group.by = "all")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& percent.mt < 25 & percent_hb<5)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:30, reduction = "pca")
pbmc <- FindClusters(pbmc, resolution = 1, cluster.name = "unintegrated_clusters")
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(pbmc,label=T,group.by = c("seurat_clusters","orig.ident"))
pbmc = IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                verbose = FALSE)

pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])

pbmc <- FindNeighbors(pbmc, reduction = "INTEGRATED.CCA", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "INTEGRATED.CCA")
DimPlot(pbmc,label = T,group.by =c( "seurat_clusters","orig.ident"), reduction = "umap")

library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce, samples="orig.ident",dbr = 0.01)
pbmc <- as.Seurat(sce)
DimPlot(pbmc,  group.by = "scDblFinder.class")

table(pbmc$scDblFinder.class)
#singlet doublet 
#96336    4441
pbmc = subset(pbmc,scDblFinder.class == "singlet")
pbmc@active.ident = pbmc$seurat_clusters
marker = FindAllMarkers(pbmc,min.pct = 0.2,logfc.threshold = 0.5)
pbmc$cell = "cell"
pbmc$cell[pbmc$seurat_clusters %in% c(17,23,10,5,32,19,27,21,30,18,22,30,29,2,3,15,20)] = "Myeloid"
pbmc$cell[pbmc$seurat_clusters %in% c(0,7,26,4,11,14,9,6,1,8,16)] = "NK&T cell"
pbmc$cell[pbmc$seurat_clusters %in% c(13,28,12,31)] = "B cell"
#pbmc$cell[pbmc$seurat_clusters %in% c(18,22,30)] = "Megakaryocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(25)] = "Plasma"
#pbmc$cell[pbmc$seurat_clusters %in% c(2,3,15,20,21)] = "Granulocytes"
#pbmc$cell[pbmc$seurat_clusters %in% c(29)] = "Mast cell"
pbmc$cell[pbmc$seurat_clusters %in% c(24)] = "Cycling cell"
DimPlot(pbmc,group.by = "subcell",label = T,reduction = "umap")
unique(pbmc$orig.ident)
pbmc$sample = "cell"
pbmc$sample[pbmc$orig.ident %in% c("Blood1","Blood3","Blood6","Blood7")] = "AAA"
pbmc$sample[pbmc$orig.ident %in% c("Blood2","Blood4")] = "rAAA"
pbmc$sample[pbmc$orig.ident %in% c("H2","H3", "H4","H5" ,"H6","H7" )] = "Normal"
DimPlot(pbmc,label = T,reduction = "umap",group.by = "cell",split.by = "sample")
save(pbmc, file= "Blood.Rdata")
unique(pbmc$cell)

load("Blood_Myeloid.Rdata")
mye = pbmc

load( "Blood.Rdata")
unique(pbmc$cell)
mye = subset(mye,cell == "NK&T cell")
pbmc$cell[colnames(pbmc) %in% colnames(mye)] = mye$cell
save(pbmc,file="Blood.Rdata")


load("Blood_Myeloid.Rdata")
mye = pbmc
load("Blood_NK&T.Rdata")
nkt = pbmc
load( "Blood.Rdata")
pbmc$subcell = pbmc$cell 
pbmc$subcell[colnames(pbmc) %in% colnames(mye)] = mye$cell
pbmc$subcell[colnames(pbmc) %in% colnames(nkt)] = nkt$cell
pbmc$subcell[pbmc$seurat_clusters %in% c(13,28)] = "Naive B"
pbmc$subcell[pbmc$seurat_clusters %in% c(12)] = "Memory B"
pbmc$subcell[pbmc$seurat_clusters %in% c(31)] = "CD14 Monocyte"
DimPlot(pbmc,group.by = "subcell",label = T,reduction = "umap")
save(pbmc,file="Blood.Rdata")
pbmc@active.ident = factor(pbmc$subcell)
marker = FindAllMarkers(pbmc,min.pct = 0.2,logfc.threshold = 0.5)
pbmc = subset(pbmc,cell %in% c("EC","SMR","Fibroblast"))
save(pbmc,file = "Aorta SMC.Rdata")
