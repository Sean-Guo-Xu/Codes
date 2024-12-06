
library(Seurat) 
library(SCENIC)
library(SummarizedExperiment)
library(RColorBrewer)
library(doParallel)
library(SCopeLoomR)
setwd("F:\\AAA_scRNA\\scenic\\mouse")
load("F:\\AAA_scRNA\\mouse_aaa_myeloid.Rdata")
cellInfo <-  pbmc@meta.data[,c("cell","nFeature_RNA","nCount_RNA")]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)

regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

 
identical(colnames(regulonAUC), colnames(pbmc))
FeaturePlot(pbmc,features = rownames(regulonAUC)[1])
regulonActivity_byCellType <- sapply(split(rownames(pbmc@meta.data), pbmc$cell), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
max(regulonActivity_byCellType_Scaled)
min(regulonActivity_byCellType_Scaled)
pdf("mouse_scenic_heatmap.pdf",width = 3.6,height = 11)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,  color=colorRampPalette(c("#90d1a4","white", "#FB9A99"))(100), breaks=seq(-3, 3, length.out = 100),
                                     treeheight_row=10, treeheight_col=10, border_color=NA,cluster_rows = T)
dev.off()

setwd("F:\\AAA_scRNA\\scenic\\human")
load("F:\\AAA_scRNA\\human_aaa_myeloid.Rdata")
cellInfo <-  pbmc@meta.data[,c("cell","nFeature_RNA","nCount_RNA")]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom)

regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)


identical(colnames(regulonAUC), colnames(pbmc))
FeaturePlot(pbmc,features = rownames(regulonAUC)[1])
regulonActivity_byCellType <- sapply(split(rownames(pbmc@meta.data), pbmc$cell), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
max(regulonActivity_byCellType_Scaled)
min(regulonActivity_byCellType_Scaled)
pdf("human_scenic_heatmap.pdf",width = 3.6,height = 8)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,  color=colorRampPalette(c("#90d1a4","white", "#FB9A99"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,cluster_rows = T)
dev.off()

setwd("F:\\AAA_scRNA\\scenic\\mouse")
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
mouse_auc = get_regulons_AUC(loom)
setwd("F:\\AAA_scRNA\\scenic\\human")
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
human_auc = get_regulons_AUC(loom)
mouse_auc = assay(mouse_auc)
human_auc = assay(human_auc)
humanname = cbind(rownames(human_auc),gsub(" \\([0-9]*g\\)","",rownames(human_auc)))
mousename = cbind(rownames(mouse_auc),gsub(" \\([0-9]*g\\)","",rownames(mouse_auc)))
mousename[,2] = toupper(mousename[,2])
humanname[,2] = toupper(humanname[,2])
overlap = merge(mousename,humanname,by="V2")
overlap = overlap[,2:3]
colnames(overlap) = c("mouse","human")
human_auc = human_auc[rownames(human_auc) %in% overlap$human,]
mouse_auc = mouse_auc[rownames(mouse_auc) %in% overlap$mouse,]
setwd("F:\\AAA_scRNA\\scenic")
write.table(overlap,"h_m_overlapTF.txt",sep = "\t",quote = F,row.names = F,col.names = F)
#############featureplot#
library(ggplot2)
library(gcookbook)
library(tidyverse)

setwd("F:\\AAA_scRNA\\scenic\\human")
load("F:\\AAA_scRNA\\human_aaa_myeloid.Rdata")
identical(colnames(human_auc), colnames(pbmc))
pbmc@meta.data = cbind(pbmc@meta.data,t(human_auc))
for (i in 1:length(rownames(human_auc))) {

umap = pbmc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = pbmc@meta.data$cell) 
p=FeaturePlot(pbmc,features = rownames(human_auc)[i])
p=p+ geom_segment(data =umap,aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                                 xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
                  colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(data = umap,aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                               xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ,angle=90) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p = p +scale_colour_gradient(low = "#90d1a4", high ="#FB9A99" )+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标???
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景???
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title ,
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=10), #设置legend标签的大???
        legend.key.size=unit(1,'cm'))
ggsave(paste("human",rownames(human_auc)[i],".tiff",sep = "_"),p,width =5.5 ,height = 5)
}
setwd("E:\\AAA_scRNA\\scenic\\mouse")
load("E:\\AAA_scRNA\\mouse_aaa_myeloid.Rdata")
identical(colnames(mouse_auc), colnames(pbmc))
pbmc@meta.data = cbind(pbmc@meta.data,t(mouse_auc))
for (i in 1:length(rownames(mouse_auc))) {
  
  umap = pbmc@reductions$umap@cell.embeddings %>%  
    as.data.frame() %>% 
    cbind(cell_type = pbmc@meta.data$cell) 
  p=FeaturePlot(pbmc,features = rownames(mouse_auc)[i])
  p=p+ geom_segment(data =umap,aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
                    colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
    geom_segment(data = umap,aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                                 xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
             color="black",size = 3, fontface="bold" ,angle=90) + 
    annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
             color="black",size = 3, fontface="bold" ,angle=90) 
  p = p +scale_colour_gradient(low = "#90d1a4", high ="#FB9A99" )+
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标???
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景???
          plot.background=element_rect(fill="white"),
          legend.title = element_blank(), #去掉legend.title ,
          legend.key=element_rect(fill='white'), #
          legend.text = element_text(size=10), #设置legend标签的大???
          legend.key.size=unit(1,'cm'))
  ggsave(paste("mouse",rownames(human_auc)[i],".tiff",sep = "_"),p,width =5.5 ,height = 5)
}
tf = read.table("F:\\AAA_scRNA\\scenic\\Cebpb.txt",sep = "\t")

col = c("#f57245","#e6f49a","#A6CEE3","#fbe18c","#90d1a4","#33A02C","#FB9A99","#feac60","#B2DF8A")
names(col) = c("Mono","Foamy/TREM2hi-Mø","TRM/FLOR2hi-Mø","Inf-Mø","cDC1","cDC2","IFNIC-Mø","ProInf-Mø","Foamy/SPP1hi-Mø")
library(Seurat)
setwd("F:\\AAA_scRNA\\scenic")
load("F:\\AAA_scRNA\\human_aaa_myeloid.Rdata")
pbmc@meta.data = cbind(pbmc@meta.data,t(human_auc))
library(ggplot2)
for(i in overlap$human)
{
  RidgePlot(pbmc, features = i ,group.by="cell") +scale_fill_manual(values = col)+ theme(legend.position='none')
  ggsave(paste(i,".tiff"),width =7,height = 4)
  }
load("F:\\AAA_scRNA\\mouse_aaa_myeloid.Rdata")

col = c("#f57245","#e6f49a","#A6CEE3","#fbe18c","#90d1a4","#33A02C","#FB9A99","#feac60","#B2DF8A")
names(col) = c("Mono","Foamy-Mø/Trem2hi","TRM","Inf-Mø","cDC1","cDC2","IFNIC-Mø","Pro-Inf","Foamy-Mø/Spp1hi")

pbmc@meta.data = cbind(pbmc@meta.data,t(mouse_auc))
for(i in overlap$mouse)
{
  RidgePlot(pbmc, features = i ,group.by="cell") +scale_fill_manual(values = col)+ theme(legend.position='none')
  ggsave(paste(i,".tiff"),width =7,height = 4)
  
}
