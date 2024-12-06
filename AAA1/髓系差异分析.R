setwd("E:\\AAA_scRNA")
library(Seurat)
load("6Normal.Rdata")
con = pbmc
load("8AAA.Rdata")
pbmc = merge(pbmc,con)
pbmc = subset(pbmc,cell %in% c("Monocyte","Macrophage"))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
artery = FindMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2,ident.1 = "AAA" ,ident.2 = "Normal" ,group.by = "sample")
library(Seurat)
load("6Normal.Rdata")
con = pbmc
setwd("E:\\AAA_scRNA\\blood")
load("N_blood.Rdata")
pbmc$sample = "Normal"
con = pbmc
load("A_blood.Rdata")
pbmc$sample = "AAA"
pbmc = merge(pbmc,con)
pbmc = subset(pbmc,cell %in% c("Monocyte","Macrophage"))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
blood = FindMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2,ident.1 = "AAA" ,ident.2 = "Normal" ,group.by = "sample")
artery$gene = rownames(artery)
blood$gene = rownames(blood)
arteryup = artery[artery$avg_log2FC>0,]
bloodup =blood[blood$avg_log2FC>0, ]
arterydown = artery[artery$avg_log2FC<0,]
blooddown =blood[blood$avg_log2FC<0, ]

up  = merge(arteryup, bloodup,by="gene")
down = merge(arterydown, blooddown ,by="gene")
write.csv(up,"up.csv",col.names = T,row.names = F)
write.csv(down,"down.csv",row.names = F)
