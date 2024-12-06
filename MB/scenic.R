library(SCENIC)
library(Seurat)
library(GENIE3)
library(AUCell)
library(data.table)
library(RcisTarget)
library(loomR) 
setwd("E:\\scdata")
exprMat = as.matrix(pbmc@assays$RNA@data)
data(list="motifAnnotations_hgnc_v9",package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
hg19_dbs <- list('500bp'= 'hg19-500bp-upstream-7species.mc9nr.feather', 
                 '10kb' = 'hg19-tss-centered-10kb-7species.mc9nr.feather')

scenicOptions<- initializeScenic(org='hgnc',dbs =hg19_dbs,
                                 datasetTitle='Analysis title', 
                                 nCores=4,dbDir = 'E:/scdata/cisTarget_databases')
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
library(doParallel)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1) 
save(exprMat_filtered_log ,file="log_exp.Rdata")
load("log_exp.Rdata")
rm(pbmc,exprMat_filtered ,exprMat)
source("gene3.R")
gc()
closeAllConnections()
gene3(exprMat_filtered_log, finished = 36,scenicOptions,nParts = 100)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
source("step2.R")
scenicOptions <- step2(scenicOptions, coexMethod=c("top10perTarget"))
exprMat_log <- log2(exprMat+1)
library(BiocParallel)
library(doMC)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
