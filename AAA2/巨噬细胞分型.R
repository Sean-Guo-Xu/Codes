library(Seurat)
library(clusterProfiler) 
library(AUCell)
library("EnrichmentBrowser")
library("KEGGREST")
setwd("D:\\AAA\\scRNA")
load("Macro seurat.Rdata")
geneset = list(TREM2=c("TREM2","CD9","LGALS3","CTSB","SPP1"))
memory.limit(size=100000)
cells_rankings <- AUCell_buildRankings(pbmc@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- as.numeric(getAUC(cells_AUC))
pbmc$TREM2  <- aucs
library(ggraph)
library(paletteer)
library(ggsci)

ggplot(data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=TREM2)
) + geom_point( size=1.5
) + scale_color_viridis(option="C")  + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
VlnPlot(pbmc,features = c('TREM2'),  pt.size = 0, adjust = 2,group.by = "type",combine = T)+scale_fill_npg()+theme(legend.position = "none")+stat_summary(fun.data = "mean_sdl",
                                                                                                                                                        fun.args = list(mult = 1), geom = "pointrange", color = "white")
library(monocle)       
load("Macro cds.Rdata")
cds$TREM2 = pbmc$TREM2

plot_cell_trajectory(cds,color_by = "TREM2",show_branch_points = F)

data = cbind(pbmc$sphingolipid,pbmc$TREM2,cds$type)
data = as.data.frame(data)
colnames(data) = c("sphingolipid","Pseudotime","Cluster")
data$Cluster = as.factor(data$Cluster)
pbmc$seurat_clusters = as.factor(pbmc$seurat_clusters)
data$sphingolipid=as.numeric(data$sphingolipid)
data$Pseudotime =as.numeric(data$Pseudotime)
library(ggpubr)
ggplot(data, aes(x=Pseudotime, y=sphingolipid,color= Cluster))+geom_point()+scale_color_npg()+theme_bw()+ 
  stat_smooth(aes(x=Pseudotime, y=sphingolipid,color=c("Line")),data,method=lm,formula=y~x,color="black")+stat_cor(aes(x =Pseudotime, y =sphingolipid,color =c("Line")),color = "black")

library(Seurat)
setwd("D:\\AAA\\scRNA")
load("Macro seurat.Rdata")
pbmc = subset(pbmc,type %in% c("M2","M0"))
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 13
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.45)                  #??????????,??????????????
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)   
DimPlot(pbmc,label = T)
DimPlot(pbmc,label = T,group.by = "state")
VlnPlot(pbmc,features = c('TREM2'),  pt.size = 0, adjust = 2,group.by = "seurat_clusters",combine = T)+scale_fill_npg()+theme(legend.position = "none")+stat_summary(fun.data = "mean_sdl",
                                                                                                                                                      fun.args = list(mult = 1), geom = "pointrange", color = "white")
pbmc$type[pbmc$seurat_clusters %in% 2]=""   
ggplot(data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=sphingolipid)
) + geom_point( size=1.5
) + scale_color_viridis(option="C")  + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

######################¹ì¼£##############
library(monocle)  
library(Seurat)
library(ggpubr)
setwd("D:\\AAA\\scRNA")
load("Macro seurat.Rdata")

data = cbind(pbmc$sphingolipid,pbmc$TREM2,pbmc$seurat_clusters)
data = as.data.frame(data)
colnames(data) = c("sphingolipid","Pseudotime","Cluster")
data$Cluster = as.factor(data$Cluster)
pbmc$seurat_clusters = as.factor(pbmc$seurat_clusters)
plot_cell_trajectory(cds, color_by = "TREM2",show_branch_points =F) 
ggplot(data, aes(x=Pseudotime, y=sphingolipid,color= Cluster))+geom_point()+scale_color_npg()+theme_bw()+ 
  stat_smooth(aes(x=Pseudotime, y=sphingolipid,color=c("Line")),data,method=lm,formula=y~x,color="black")+stat_cor(aes(x =Pseudotime, y =sphingolipid,color =c("Line")),color = "black")



