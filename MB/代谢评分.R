.libPaths("C:\\Users\\56988\\Documents\\R\\win-library\\4.1")
setwd("E:\\brain")
library(phangorn)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(AUCell)
library(GSEABase)
library(GSVA)
library(Seurat)
load("4type_withmeta.Rdata")
memory.limit(size = 10000000)
score = pbmc@assays$METABOLISM$score
score = t(score)
load("meta.Rdata")
pheatmap::pheatmap(score,show_rownames = F,show_colnames = F)
part = subset(pbmc,cell %in% mycell)
partscore = score[]

i = 1
data =data.frame(AUcell = score[,i],Week = pbmc$batch)
data = as.data.frame(data)
my_comparisons = list(c("PCW8","PCW9"),c("PCW9","PCW12"),c("PCW12","PCW13"),c("PCW13","PCW14"),c("PCW14","PCW15"),c("PCW15","PCW16"),c("PCW16","PCW17"))
library(ggpubr)
ggplot(data=data,aes(x=Week,
                         y=AUcell))+geom_boxplot()  +
  stat_compare_means( comparisons = my_comparisons,method = "wilcox.test") +ggtitle(colnames(score)[i])
table(data$AUcell,data$Week)

mycell = c( "Neuron Stem cell", "GCP_progenitor" , "UBC_progenitor")
mycell = unique(pbmc$cell)
part = subset(pbmc,cell %in% mycell)
data =data.frame(AUcell = score[pbmc$cell %in% mycell,i],Cell= part$cell)
ggplot(data=data,aes(x=Cell,
                     y=AUcell))+geom_boxplot() 
pbmc$meata = data$AUcell
FeaturePlot(pbmc,features = "meata")

part@active.ident = as.factor(part$cell)
diff = FindAllMarkers(part,min.pct = 0.25,logfc.threshold = 1)
diff=FindMarkers(part,ident.1 = "GCP_progenitor" ,ident.2 = "UBC_progenitor",min.pct = 0.25,logfc.threshold = 1,group.by = "cell")
