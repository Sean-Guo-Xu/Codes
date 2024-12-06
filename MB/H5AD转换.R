
library(anndata)
library(reticulate)
setwd("E:\\brain\\Normal2")
use_condaenv(condaenv = "r-h5ad", required = TRUE)
ad <- import("anndata")
final_ad <- ad$read_h5ad("aldinger20.processed.h5ad")
colnames(final_ad$X) <- final_ad$var$features

pbmc<- Seurat::CreateSeuratObject(counts = final_ad$X, assay = "RNA",
                                         meta.data = final_ad$obs,
                                         min.cells = 3, min.features = 200)
pbmc@meta.data = final_ad$obs
save(pbmc,file="normal2.Rdata")

