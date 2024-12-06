color11 =c("#FFFF00", "#9B59B6", "#16A085", "#F1C40F", "#56B4E9", "#E67E22", "#C0392B","#2980B9", "#27AE60","#8E44AD", "#2C3E50")
cellnames = c("Immune cells","Neuron", "Oligodendrocyte" , "Astrocyte" ,"GP4","GP4", "GP3or4", "SHH","WNT","Endothelial","Pericytes")
names(color11) = cellnames

setwd("E:\\brain")
library(Seurat)
mat = read.csv("GSE155446_human_raw_counts.csv",row.names=1)
pbmc = CreateSeuratObject(mat)
metadata = read.csv("GSE155446_human_cell_metadata.csv")
pbmc@meta.data=cbind(pbmc@meta.data,metadata)
pbmc$cell = pbmc$coarse_cell_type 
pbmc$sample = pbmc$subgroup
pbmc1=pbmc
load("E:\\brain\\Group4\\mydata.Rdata")
pbmc = merge(pbmc1,pbmc)
unique(pbmc$sample)
pbmc  =JoinLayers(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc$all = "all"
VlnPlot(pbmc,features = c("percent.mt","percent_hb","nFeature_RNA","percent_hb"),group.by = "all")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <7000 & percent.mt < 30& percent_hb<5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 2, cluster.name = "unintegrated_clusters")
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
cancermakrer = FindAllMarkers(pbmc,logfc.threshold = 1,min.pct = 0.2)
DimPlot(pbmc, group.by = c("seurat_clusters","orig.ident"),label = T)
pbmc$cell[pbmc$cell %in% "malignant"] = pbmc$subgroup[pbmc$cell %in% "malignant"]
save(pbmc,file="cancer.Rdata")

library(SingleR)
library(celldex)
load("E:\\brain\\cancer.Rdata")
pbmc = subset(pbmc,coarse_cell_type %in% c( "lymphocytes" ,"macrophage_monocytes" ))
save(pbmc,file="cancer.Rdata")
neuron = subset(pbmc,seurat_clusters %in% c(49))

#65 Endothelial
#63 Pericytes
save(neuron,file="neuron cell.Rdata")
im = subset(pbmc,seurat_clusters %in% c(13,58,23,54,47,39,26,48,60))
save(im,file = "im cell.Rdata")
ref = subset 
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =pbmc@assays$RNA@data, ref = ref, labels = ref$label.main)
pbmc$SingleR = cellpred$labels

pbmc$cell[pbmc$seurat_clusters %in% c(9,23,34)] ="Interneuron" #KIT LHX5 PAX2
pbmc$cell[pbmc$seurat_clusters %in% c(26,35)] ="Molecular Layer Interneuron" #NXPH1 CXCL12
pbmc$cell[pbmc$seurat_clusters %in% c(10.29,32,21,27,28)] ="Neuron Stem cell" #SOX2 NES
pbmc$cell[pbmc$seurat_clusters %in% c(1,18,3,6,15)] ="Granule cell" ## GAP43 GRIK2 SEMA6A
pbmc$cell[pbmc$seurat_clusters %in% c(16)] ="Middle Brain Neuron" ##
pbmc$cell[pbmc$seurat_clusters %in% c(2,4,5,8,11,31)] ="GCP" # ATOH1
pbmc$cell[pbmc$seurat_clusters %in% c(13,17,24,30)] ="Purkinje cell" # RORA PCP4
pbmc$cell[pbmc$seurat_clusters %in% c(0,25,14,20)] ="Unipolar Brush cell" #OTX2 RSPO3 EOEMS 
pbmc$cell[pbmc$seurat_clusters %in% c(5,12,22)] ="GCP_progenitor" #RFC3
pbmc$cell[pbmc$seurat_clusters %in% c(19,33)] ="UBC_progenitor" #RFC3
pbmc$cell[pbmc$seurat_clusters %in% c(10,29)] ="Neuron Stem cell" #SOX11 CTNNB1 HNRNPH1
pbmc$cell[pbmc$seurat_clusters %in% c(7)] ="choroid plexus progenitor cells" #SOX11 CTNNB1 HNRNPH1

library(ggplot2)
DimPlot(pbmc,group.by = "cell",label = T)
gene = c("KIT", "LHX5" ,"PAX2","NXPH1","SLC1A3","SOX2","NES","GAP43","GRIK2","MGP","ATOH1",
         "RORA","PCP4","OTX2","RSPO3","EOMES","RFC3","SOX11","CTNNB1","SPP1")
DotPlot(pbmc,group.by = "seurat_clusters",features = gene) +
  theme(axis.text.x = element_text(angle= 45 , vjust= 1 , hjust= 1 )) 

########去除双细胞######
library(scDblFinder)
library(SingleCellExperiment)
library(Seurat)
library(BiocParallel)
sce <- as.SingleCellExperiment(pbmc)
setwd("E:\\brain")
load("E:\\brain\\cancer.Rdata")
sce <- scDblFinder(sce, samples="orig.ident")
pbmc <- as.Seurat(sce)
DimPlot(pbmc,  group.by = "scDblFinder.class")
table(pbmc$scDblFinder.class)
#singlet doublet 
#87585    7254 
pbmc = subset(pbmc,scDblFinder.class == "singlet")
save(pbmc,file = "cancer.Rdata")

see=available_outcomes()
gwas= see[grep("HLA",see$trait),]
ebi-a-GCST90002106
exposure_dat <-extract_instruments(outcomes='ebi-a-GCST90002106')
