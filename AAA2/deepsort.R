library(Seurat)
setwd("D:\\AAA\\scRNA")
load("7AAA.Rdata")
pbmc = pbmc[["RNA"]]@data
pbmc =as.matrix(pbmc)
pbmc = pbmc[,1:2000]
write.csv(pbmc,file = '7AAA.csv')
