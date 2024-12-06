library(Seurat)
setwd("E:\\AAA_scRNA\\AAA")
dir = list.dirs()
pbmc=Read10X(dir[2])
gene = rownames(pbmc)
pbmc=CreateSeuratObject(pbmc)
pbmc$sample = gsub("./", "", dir[2])

for (i in 3:23) {
  part = Read10X(dir[i])
  gene= intersect(gene,rownames(part))
  part=CreateSeuratObject(part)
  part$sample = gsub("./", "", dir[i])
  pbmc = merge(pbmc,part)
}

pbmc = JoinLayers(pbmc)
pbmc = subset(pbmc,features = gene)
save(pbmc,file="E:\\AAA_scRNA\\21_AAA.Rdata")

setwd("E:\\AAA_scRNA\\Control")
dir = list.dirs()
pbmc=Read10X(dir[2])
gene = rownames(pbmc)
pbmc=CreateSeuratObject(pbmc)
pbmc$sample = gsub("./", "", dir[2])

for (i in 3:12) {
  part = Read10X(dir[i])
  gene= intersect(gene,rownames(part))
  part=CreateSeuratObject(part)
  part$sample = gsub("./", "", dir[i])
  pbmc = merge(pbmc,part)
}
pbmc = JoinLayers(pbmc)
pbmc = subset(pbmc,features = gene)
save(pbmc,file="E:\\AAA_scRNA\\11_Control.Rdata")

setwd("E:\\AAA_scRNA")
load("22_AAA.Rdata")
pbmc$type = "AAA"
AAA = pbmc

load("12_Control.Rdata")
pbmc$type = "Control"
gene = intersect(pbmc@assays$RNA@features[[1]],AAA@assays$RNA@features[[1]])
pbmc=merge(AAA,pbmc)
pbmc = JoinLayers(pbmc)
pbmc = subset(pbmc,features = gene)
VlnPlot(pbmc,features = c("LDAH"),group.by = "cell")
save(pbmc,file="33_all.Rdata")
library(ggpubr)
VlnPlot(pbmc,features = "LDAH",raster = F,split.by ="type",cols = c("#E07B54","#6BB7CA"))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")

load("32_all.Rdata")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,raster=FALSE)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5500)
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$sample)
pbmc = NormalizeData(pbmc)
pbmc = FindVariableFeatures(pbmc)
pbmc = ScaleData(pbmc)
pbmc  = RunPCA(pbmc)
pbmc = FindNeighbors(pbmc, dims=1:30,reduction="pca")
pbmc <- FindClusters(pbmc, resolution = 2, cluster.name = "unintegrated_clusters")
pbmc <- IntegrateLayers(object = pbmc, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc = RunUMAP(pbmc, dims = 1:30, reduction = "integrated.cca")
save(pbmc,file="integrated.Rdata")
markers = FindAllMarkers(pbmc,min.pct = 0.2, logfc.threshold = 1)
save(markers,file = "all_markers.Rdata")
DimPlot(pbmc,group.by = "seurat_clusters",label = TRUE)

NK = c("KLRD1","NCAM1","NKG7") #26
T = c("CD3D","CD3E","CD3G") #12,13,8
DC= c("CD1C","FCER1A","HLA-DRA")#23
Mast = c("KIT","CPA3","TPSB2") #24
Endothelial = c("PECAM1","PLVAP","PTPRB") #5,30,15,19
Fibroblast = c("DCN","COL1A1","COL1A2") #3,6,31,25
SMC = c("MYH11","ACTA2","MYL9",,"LMOD1") #2,7,21,17
Monocyte = c("FCN1", "APOBEC3A", "THBS1") #16,18
Macrophage = c("CD163", "CD68", "CD14") #1,4,11,27,10
B = c("CD79A","MZB1","MS4A1") #0,9
Proliferation = c("MKI67","TOP2A","CDC20") #22
Plasma = c("CD38","XBP1") #20
Erythrocytes=c("HBA1","HBA2","HBD") #28

pbmc$cell=pbmc$sample
pbmc$cell[pbmc$seurat_clusters %in% c(26)] = "NK cell"
pbmc$cell[pbmc$seurat_clusters %in% c(13,12,8)] = "T cell"
pbmc$cell[pbmc$seurat_clusters %in% c(24)] = "Mast cell"
pbmc$cell[pbmc$seurat_clusters %in% c(5,30,15,19)] = "Endothelial"
pbmc$cell[pbmc$seurat_clusters %in% c(3,6,31,25)] = "Fibroblast"
pbmc$cell[pbmc$seurat_clusters %in% c(14,2,7,21,17)] = "SMC"
pbmc$cell[pbmc$seurat_clusters %in% c(16,18)]="Monocyte"
pbmc$cell[pbmc$seurat_clusters %in% c(1,4,27,11,10)]="Macrophage"
pbmc$cell[pbmc$seurat_clusters %in% c(0,9,29)] = "B cell"
pbmc$cell[pbmc$seurat_clusters %in% c(22)] = "Proliferation"
pbmc$cell[pbmc$seurat_clusters %in% c(28)] = "Erythrocytes"
pbmc$cell[pbmc$seurat_clusters %in% c(20)] = "Plasma"
pbmc$cell[pbmc$seurat_clusters %in% c(23)] = "DCs"

DimPlot(pbmc,group.by=c("seurat_clusters","cell"),label = TRUE)

pbmc$set = pbmc$sample
pbmc$set = pbmc$set[pbmc$sample %in% c("Control1","Control2","Control3","Control4","Control5","Control6",
                                       "GSE226492_con1","GSE226492_con2","GSE226492_con3",
                                       "GSM5077731","GSM5077732")] = "Control"
  pbmc$set = pbmc$set[ !(pbmc$sample %in% c("Control1","Control2","Control3","Control4","Control5","Control6",
                                           "GSE226492_con1","GSE226492_con2","GSE226492_con3",
                                           "GSM5077731","GSM5077732")) ] = "AAA" 