setwd("E:\\brain")
library(SeuratDisk)
library(Seurat)
Convert("aldinger20.processed.h5ad", "h5seurat",
        overwrite = TRUE,assay = "RNA")
 LoadH5Seurat("aldinger20.processed.h5seurat", assays="RNA")
 