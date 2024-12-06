setwd("E:\\AAA_scRNA")
library(Seurat)
load("first_batch.Rda")
unique(seurat_object$orig.ident)
meta=seurat_object@meta.data
meta$index=rownames(meta)
cell=read.table("E:\\AAA_scRNA\\scRNA\\first_batch.txt",header=T,sep="\t")
meta=merge(meta,cell,by="index")
rownames(meta)=meta$index
seurat_object@meta.data=meta
save(seurat_object,file="first_batch.Rdata")
################

load("second_batch.Rda")
unique(seurat_object$orig.ident)
meta=seurat_object@meta.data
meta$index=rownames(meta)
cell=read.table("E:\\AAA_scRNA\\scRNA\\second_batch.txt",header=T,sep="\t")
meta=merge(meta,cell,by="index")
rownames(meta)=meta$index
seurat_object@meta.data=meta
save(seurat_object,file="second_batch.Rdata")
DimPlot(seurat_object,group.by = "cell_type")

##############################
load("first_batch.Rdata")
aaa1=subset(seurat_object,orig.ident %in% unique(seurat_object$orig.ident)[1:3])
con1=subset(seurat_object,orig.ident %in% unique(seurat_object$orig.ident)[4:6])
load("second_batch.Rdata")
aaa2=subset(seurat_object,orig.ident %in% unique(seurat_object$orig.ident)[1:5])
con2=subset(seurat_object,orig.ident %in% unique(seurat_object$orig.ident)[6:8])
con=merge(con1,con2)
aaa=merge(aaa1,aaa2)
con$sample = rep("Normal",length(con$index))
aaa$sample = rep("AAA",length(aaa$index))
save(con,file="6Normal.Rdata")
save(aaa,file="8AAA.Rdata")
