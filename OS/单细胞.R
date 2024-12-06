library(Seurat)
setwd("E:\\bladder")

library(SeuratDisk)
Convert("GSE169379_MIBC_snSeq.h5ad", dest="h5seurat",
        assay = "RNA",
        overwrite=F)
pbmc = LoadH5Seurat("GSE169379_MIBC_snSeq.h5seurat")
library(Seurat)
pbmc = UpdateSeuratObject(pbmc)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #�����ڽӾ���
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #��ϸ������,�Ż���׼ģ�黯
pbmc <- RunUMAP(pbmc,dims = 1:20)

DimPlot(pbmc,group.by = c("seurat_clusters","celltype"),label=T)
library(harmony)
pbmc = RunHarmony(pbmc,"batch", plot_convergence = TRUE)#��ʱ1min
#ϸ������
###resolustion AAA 0.8 NORMAL 0.6#######
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  identity()
DimPlot(pbmc,label = T,group.by = "orig.ident")
DimPlot(pbmc,label = T,group.by = c("batch","celltype","seurat_clusters"))

library(ggsci)
library(ggplot2)
p1=DimPlot(pbmc,group.by = c("celltype"),label = T)+scale_color_frontiers()+labs(title = "")
p2 =FeaturePlot(pbmc,features = c("EFEMP1","SVIL"),col=c("#6F286AFF","#EFD500FF"))
p3 = VlnPlot(pbmc,group.by = "celltype",features = c("EFEMP1"))+scale_fill_frontiers()
p1
ggsave("sc_tumor_all.tiff",width = 6.5,height = 5)
p2
ggsave("sc_tumor_feature.tiff",width = 10,height = 5)
p3
bladder = pbmc
ggsave("sc_tumor_violin.tiff",width = 5,height = 3)
Convert("GSE169379_non_tumor_snSeq.h5ad", dest="h5seurat",
        assay = "RNA",
        overwrite=F)
pbmc = LoadH5Seurat("GSE169379_non_tumor_snSeq.h5seurat")
library(Seurat)
pbmc = UpdateSeuratObject(pbmc)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #�����ڽӾ���
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #��ϸ������,�Ż���׼ģ�黯
pbmc <- RunUMAP(pbmc,dims = 1:20)

DimPlot(pbmc,group.by = c("seurat_clusters","celltype"),label=T)
library(harmony)
pbmc = RunHarmony(pbmc,"batch", plot_convergence = TRUE)#��ʱ1min
#ϸ������
###resolustion AAA 0.8 NORMAL 0.6#######
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  identity()
DimPlot(pbmc,label = T,group.by = c("batch","celltype","seurat_clusters"))
pbmc$celltype = as.character(pbmc$SubType_Normal)
pbmc = subset(pbmc, SubType_Normal != "Normal_Unk")
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Fibroblast_1","Normal_Fibroblast_2")] = "Fibroblast"
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Endothelial_1","Normal_Endothelial_2")] = "Endothelial"
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Basal_Intermediate","Normal_CDH12","Normal_KRT7_KRT13_Basal","Normal_Umbrella_Intermediate")] ="Epithelial"
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Myeloid")] = "Myeloid"
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Bcell","Normal_Mast_cell","Normal_Tcell")] = "Lymphocyte"
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Neuronal")] = "Neuronal"
pbmc$celltype[pbmc$SubType_Normal %in% c("Normal_Smooth_muscle")]="SMC"

library(ggsci)
library(ggplot2)
p1=DimPlot(pbmc,group.by = c("celltype"),label = T)+scale_color_frontiers()+labs(title = "")
p2=FeaturePlot(pbmc,features = c("EFEMP1","SVIL"),col=c("#6F286AFF","#EFD500FF"))
p3 = VlnPlot(pbmc,group.by = "celltype",features = c("EFEMP1"))+scale_fill_frontiers()
p1
ggsave("sc_normal_all.tiff",width = 6.5,height = 5)
p2
ggsave("sc_normal_feature.tiff",width = 10,height = 5)
p3
ggsave("sc_normal_violin.tiff",width = 6,height = 3)
bladder$sample = rep("Tumor",nrow(bladder@meta.data))
pbmc$sample = rep("Normal",nrow(pbmc@meta.data))
load("bladder_scRNA.Rdata")
all=merge(bladder,pbmc)
library(AUCell)
gene = c("EFEMP1","SVIL")
for (i in gene) {
  

n = cbind(FetchData(pbmc,vars = c(paste("rna",i,sep = "_"))),pbmc$celltype,pbmc$sample)
t = cbind(FetchData(bladder,vars =c(paste("rna",i,sep = "_"))),bladder$celltype,bladder$sample)
colnames(n) = c("GeneExp","CellType","Sample")
colnames(t) = c("GeneExp","CellType","Sample")
data = rbind(n,t)
library(ggplot2)
library(ggpubr)
data = data[!(data$CellType %in% c("Neuronal","SMC")),]
ggplot(data, aes(x=CellType, y=GeneExp,fill=Sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                             method = "t.test")+
  geom_violin(trim=FALSE) + geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_manual(values = c("#0094CDFF","#D51317FF"))+ 
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显�?(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属�?
                                 size=12),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属�?
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  xlab("")
ggsave(paste(i,"N_T_Violin.tiff"),width = 5,height = 3)}
