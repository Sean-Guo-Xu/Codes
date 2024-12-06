setwd("D:\\AAA\\scenic")
load("E:\\AAA_scRNA\\aaa_hg_mye.Rdata")
library(Seurat)
exprMat = as.matrix(pbmc@assays$RNA@data)

library(SCENIC)
data(list="motifAnnotations_hgnc", package="RcisTarget")
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir="./cisTarget_databases", nCores=10) 
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

###3
setwd("D:\\mouse")
load("Mmouse16.Rda")
library(Seurat)
exprMat = as.matrix(pbmc@assays$RNA@data)
cellInfo <-  pbmc@meta.data

head(cellInfo)
table(cellInfo$CellType)

library(SCENIC)
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org="mgi", 
                                  dbDir="./", nCores=6)
scenicOptions@inputDatasetInfo$cellInfo =cellInfo
genesKept <- geneFiltering(exprMat, scenicOptions,npart=20)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

##############
setwd("D:\\human")
load("aaa_hg_mye.Rdata")
library(Seurat)
exprMat = as.matrix(pbmc@assays$RNA@data)

library(SCENIC)
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir="./cisTarget_databases", nCores=8) 
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
rm(exprMat)
rm(pbmc)
runGenie3(exprMat_filtered_log, scenicOptions,nPart = 50)


memory.limit()



                   