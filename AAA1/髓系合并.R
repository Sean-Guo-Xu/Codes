setwd("F:\\AAA_scRNA")
library(Seurat)
load("mouse_aaa_myeloid.Rdata")
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
DimPlot(pbmc,group.by = "cell",label=T)
musGenes  = pbmc@assays$RNA@counts@Dimnames[[1]]

genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
               values = musGenes, 
               mart = mouse, 
               attributesL = c("hgnc_symbol"), 
               martL = human, uniqueRows=T) 
11781/13079
genes = genes[!duplicated(genes$MGI.symbol),]
genes = genes[!duplicated(genes$HGNC.symbol),]
pbmc = subset(pbmc,features = genes$MGI.symbol)
count = GetAssayData(pbmc, slot = "counts")
count@Dimnames[[1]] = genes$HGNC.symbol
mousecount = count
n=0
pbmc$sample = pbmc$orig.ident
for (i in unique(pbmc$sample)) {
  n=n+1
  pbmc$orig.ident[pbmc$sample %in% i] = paste("mouse_","AAA",n,sep = "")
}
pbmc$species = "mouse"
mousemeta= pbmc@meta.data


load("human_aaa_myeloid.Rdata")
pbmc$sample = pbmc$orig.ident
n=0
for (i in unique(pbmc$sample)) {
  n=n+1
  pbmc$orig.ident[pbmc$sample %in% i] = paste("human_","AAA",n,sep = "")
}
pbmc$species = "human"
pbmc = subset(pbmc,features = genes$HGNC.symbol)
humancount = GetAssayData(pbmc,slot = "counts")
mousepbmc=CreateSeuratObject(mousecount,meta.data = mousemeta)
mousepbmc = subset(mousepbmc,features = pbmc@assays$RNA@counts@Dimnames[[1]])
humanpbmc = CreateSeuratObject(humancount,meta.data = pbmc@meta.data)
11611/13079
11611/29874
pbmc= merge(humanpbmc,mousepbmc)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
#############################
human.mouse.list <- SplitObject(pbmc, split.by = "orig.ident")
for (i in 1:length(human.mouse.list)) {
  human.mouse.list[[i]] <- NormalizeData(human.mouse.list[[i]], verbose = FALSE)
  human.mouse.list[[i]] <- FindVariableFeatures(human.mouse.list[[i]], selection.method = "vst", 
                                                nfeatures = 1000, verbose = FALSE)
}
human.mouse.anchors <- FindIntegrationAnchors(object.list =human.mouse.list , anchor.features = 1000,dims = 1:30, k.filter=50)
human.mouse.integrated <- IntegrateData(anchorset = human.mouse.anchors, dims = 1:30,k.weight = 30)
DefaultAssay(human.mouse.integrated) <- "integrated"
integrated_MPC=human.mouse.integrated
integrated_MPC <- ScaleData(integrated_MPC, verbose = FALSE)
integrated_MPC <- RunPCA(integrated_MPC, npcs = 30, verbose = FALSE)
integrated_MPC <- RunUMAP(integrated_MPC, reduction = "pca", dims = 1:20)
integrated_MPC  <- FindNeighbors(integrated_MPC , dims = 1:20)
integrated_MPC  <- FindClusters(integrated_MPC , resolution = 0.8)
DimPlot(integrated_MPC,group.by = c("species","cell"),label = T)
DimPlot(subset(integrated_MPC,species == "human"),label = T,group.by = "orig.ident")
DimPlot(subset(integrated_MPC,sample=="mouse"),group.by = "sample",label = T)
save(integrated_MPC,file="mouse_human.Rdata")
