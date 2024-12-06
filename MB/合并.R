library(Seurat)
library(org.Hs.eg.db)
setwd("E:\\brain\\Group4")
dir = list.dirs()
pbmc=Read10X(dir[3],gene.column = 2)
pbmc=CreateSeuratObject(pbmc)
pbmc$orig.ident = "GP4-1"
dir = dir[c(5,7,9,11,13,15,17,19,21)]
for (i in 1:9) {
  part = Read10X(dir[i],gene.column = 2)
  part = CreateSeuratObject(part)
  part$orig.ident = paste("GP4-",i,sep = "")                            
  pbmc = merge(pbmc,part)
}
pbmc$sample = "GP4"
unique(pbmc$orig.ident)
save(pbmc,file="gp4.Rdata")

G3=readRDS("g3.RDS")
load("gp4.Rdata")
G3 = CreateSeuratObject(G3)
G3$sample = "GP3"
pbmc = merge(pbmc,G3)
pbmc = JoinLayers(pbmc)
unique(pbmc$orig.ident)
pbmc$orig.ident[pbmc$orig.ident %in% "016012"] = "GP3-1"
pbmc$orig.ident[pbmc$orig.ident %in% "018034"] = "GP3-2"
pbmc$orig.ident[pbmc$orig.ident %in% "058"] = "GP3-3"
pbmc$orig.ident[pbmc$orig.ident %in% "018086"] = "GP3-4"
pbmc$orig.ident[pbmc$orig.ident %in% "7"] = "GP3-5"
pbmc$orig.ident[pbmc$orig.ident %in% "019011"] = "GP3-6"
save(pbmc,file="mydata.Rdata")
