################### h5ad&seurat transform######
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
library(Seurat) 
load("F:\\AAA_scRNA\\7AAA.Rdata")
pbmc$tissue="Aorta"
aaa = pbmc


load("F:\\AAA_scRNA\\blood\\A_blood.Rdata")
pbmc$tissue = "PBMC"
pbmc = merge(aaa,pbmc)

SaveH5Seurat(subset(pbmc,tissue =="Aorta"),filename="Aorta.h5seurat", overwrite = TRUE)
Convert("Aorta.h5seurat", dest = "h5ad", overwrite = TRUE)

SaveH5Seurat(subset(pbmc,tissue =="PBMC"),filename="Blood.h5seurat", overwrite = TRUE)
Convert("Blood.h5seurat", dest = "h5ad", overwrite = TRUE)
load("F:\\AAA_scRNA\\3Normal.Rdata")
pbmc[["RNA4"]] <- as(object = pbmc[["RNA"]], Class = "Assay")
DefaultAssay(pbmc) = "RNA4"
pbmc[["RNA"]] = NULL
SaveH5Seurat(pbmc,filename="Aorta2.h5seurat", overwrite = TRUE)
Convert("Aorta2.h5seurat", dest = "h5ad", overwrite = TRUE)
setwd("D:\\AAA\\scDRS")



Convert("7AAA.h5ad", dest = "7AAA.h5seurat", overwrite = T)
pbmc <- LoadH5Seurat("7AAA.h5seurat",meta.data = T)
