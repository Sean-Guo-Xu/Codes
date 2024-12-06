library(Seurat)
setwd("D:\\AAA\\qtl_scRNA")
library(ggplot2)


###############feature plot#########################
allgene = c("EIF2B2","MLH3","NEK9")
myfeatureplot<-function(pbmc,allgene,lowcolor,highcolor,tissue){
  for (gene in allgene) {
    conver = function(x){
      x[x %in% 0] = NA
      return(x)
    }
  data = data.frame(FetchData(pbmc, vars=gene), pbmc@reductions$umap@cell.embeddings)
  data=apply(data, 2, conver)
  data = as.data.frame(data)
  colnames(data)[1]="gene"
  ggplot() +  geom_point(data=data[is.na(data[,1]),],aes(UMAP_1, UMAP_2, color=gene), size=1.5) +geom_point(data=data[!is.na(data[,1]),],aes(UMAP_1, UMAP_2, color=gene), size=1.5)+scale_colour_gradient(low = lowcolor, high = highcolor)  + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+labs(color = "Expression")+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank(),axis.text.x = element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) + xlab("")+ylab("")+ggtitle(gene)
  ggsave(paste(gene,tissue,"featureplot.tiff"),width = 5.5,height = 5)
  }}

load("E:\\AAA_scRNA\\8AAA.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"), "#6BB7CA",  "#E07B54","Artery A")
load("E:\\AAA_scRNA\\6Normal.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"), "#6BB7CA",  "#E07B54","Artery N")

load("E:\\AAA_scRNA\\Blood\\A_blood.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"), "#6BB7CA",  "#E07B54","Blood A")
load("E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"), "#6BB7CA",  "#E07B54","Blood N")


########################umapplot#####################3

myumapplot<-function(pbmc,featurename,plotname){
  library(Seurat)
  library(ggplot2)
  library(ggraph)
  library(ggrepel)
  library(gcookbook)
  library(paletteer)
  library(ggsci)
  pbmc$cell = pbmc@meta.data[,featurename]
data=data.frame(pbmc@reductions$umap@cell.embeddings,pbmc$cell)
colnames(data) = c("UMAP_1","UMAP_2","cell")
labeldata = NULL
for(i in unique(pbmc$cell)){
  subpbmc = subset(pbmc ,cell == i)
  x = mean(subpbmc@reductions$umap@cell.embeddings[,1])
  y = mean(subpbmc@reductions$umap@cell.embeddings[,2])
  labeldata  = rbind(labeldata,cbind(x,y,i))
}
labeldata = as.data.frame(labeldata)
colnames(labeldata) = c("UMAP_1","UMAP_2","CellType")
labeldata$UMAP_1  = as.numeric(labeldata$UMAP_1)
labeldata$UMAP_2 = as.numeric(labeldata$UMAP_2)
ggplot() + geom_point(data=data,aes(UMAP_1, UMAP_2, color=cell), size=1.5)+ scale_color_d3() + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+labs(color = "Cell Type")+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")+  geom_label_repel( aes(UMAP_1, UMAP_2,label = CellType,color=CellType), data=labeldata,  nudge_y = 2,  alpha = 2 )
ggsave(paste(plotname,"umap.tiff"),width = 5,height = 5)
}

load("E:\\AAA_scRNA\\6Normal.Rdata")
myumapplot(pbmc,"cell","Artery_N")
load("E:\\AAA_scRNA\\8AAA.Rdata")
myumapplot(pbmc,"cell","Artery_A")
load("E:\\AAA_scRNA\\Blood\\A_blood.Rdata")
myumapplot(pbmc,"cell","Blood_A")
load("E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
myumapplot(pbmc,"cell","Blood_N")
#####################vlnplot###########
load("E:\\AAA_scRNA\\8AAA.Rdata")

allgene = c("EIF2B2","MLH3","NEK9")
pbmc$cell[pbmc$cell %in% c("T cell","NK cell")]="NK/T cell"
aaa=pbmc
library(ggpubr)
load("E:\\AAA_scRNA\\6Normal.Rdata")
pbmc=merge(aaa,pbmc)
pbmc = subset(pbmc,cell != "Unsure")
for (i in allgene) {
VlnPlot(pbmc,features = i,group.by = "cell",split.by ="sample",cols = c("#E07B54","#6BB7CA"))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")
 ggsave(paste(i,"Artery.tiff"),width = 5,height = 5)                                                                                                                      
}

load("E:\\AAA_scRNA\\Blood\\A_blood.Rdata")
aaa=pbmc
aaa$sample = "AAA"
load("E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
pbmc$sample="Normal"
pbmc=merge(aaa,pbmc)
pbmc = subset(pbmc,cell != "Unsure")
for (i in allgene) {
  VlnPlot(pbmc,features = i,group.by = "cell",split.by ="sample",cols = c("#E07B54","#6BB7CA"))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")
  ggsave(paste(i,"blood.tiff"),width = 5,height = 5)                                                                                                                      
}
 #########################

NK = c("KLRD1","NCAM1","NKG7") #7
T = c("CD3D","CD3E","CD3G") #0 17 26
Mast = c("KIT","CPA3","TPSB2") #16
Endothelial = c("PECAM1","PLVAP","PTPRB") #24 23 3
Fibroblast = c("DCN","COL1A1","COL1A2") #12 10 13
SMC = c("MYH11","ACTA2","MYL9") #4 25
Monocyte = c("FCN1","APOBEC3A","THBS1") #2 5 8 9 19 22
Macrophage = c("CD163","CD68","CD14")
B = c("CD79A","MZB1","MS4A1") #1 6 11 18 20 21 27
Proliferation = c("MKI67","TOP2A","CDC20") #15

rep("")
load("E:\\AAA_scRNA\\8AAA.Rdata")
markers = c(T,B,NK,Mast,Endothelial,Fibroblast,SMC,Monocyte,Macrophage,Proliferation)
cluster = c(rep("NK/T",3),rep("B cell",3),rep("NK/T",3),rep("Mast cell",3)
            ,rep("EC",3),rep("Fibroblast",3),rep("SMC",3),rep("Monocyte",3),rep("Macrophage",3),rep("Proliferation",3))
library(ggplot2)
DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#6BB7CA","#E07B54"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("Artery_A_dotplot.tiff",width = 15,height = 5)

load("E:\\AAA_scRNA\\6Normal.Rdata")
markers = c(T,B,NK,Mast,Endothelial,Fibroblast,SMC,Monocyte,Macrophage)
cluster = c(rep("NK/T",3),rep("B cell",3),rep("NK/T",3),rep("Mast cell",3)
            ,rep("EC",3),rep("Fibroblast",3),rep("SMC",3),rep("Monocyte",3),rep("Macrophage",3))
library(ggplot2)
DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#6BB7CA","#E07B54"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("Artery_N_dotplot.tiff",width = 14,height = 5)

#################blood dimplot#############

NK = c("KLRD1","NCAM1","NKG7") #7http://127.0.0.1:18957/graphics/plot_zoom_png?width=1536&height=705
T = c("CD3D","CD3E","CD3G") #0 17 26
Mast = c("KIT","CPA3","TPSB2") #16
Monocyte = c("FCN1","CD14","VCAN") #2 5 8 9 19 22
DC = c("HLA-DQA1","HLA-DPB1","CLEC10A")
B = c("CD79A","CD79B","MS4A1") #1 6 11 18 20 21 27
Plasma= c("IGHG1","MZB1", "SDC1")
Megakaryocyte =c("PLA2G12A","PF4","PPBP")
Granulocytes = c("MSLN","ITLN1", "PRG4")

load("E:\\AAA_scRNA\\Blood\\A_blood.Rdata")
markers = c(T,B,NK,Mast,Plasma,Monocyte,Megakaryocyte,DC,Granulocytes)
cluster = c(rep("T cell",3),rep("B cell",3),rep("NK",3),rep("Mast cell",3)
            ,rep("Plasma",3),rep("Monocyte",3),rep("Megakaryocyte",3),rep("DC",3),rep("Granulocytes",3))

DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#6BB7CA","#E07B54"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("Blood_A_dotplot.tiff",width = 15,height = 5)

load("E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
markers = c(T,B,NK,Mast,Plasma,Monocyte,Megakaryocyte,DC,Granulocytes)
cluster = c(rep("T cell",3),rep("B cell",3),rep("NK",3),rep("Mast cell",3)
            ,rep("Plasma",3),rep("Monocyte",3),rep("Megakaryocyte",3),rep("DC",3),rep("Granulocytes",3))
DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#6BB7CA","#E07B54"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("Blood_N_dotplot.tiff",width = 14,height = 5)

pbmc$cell[pbmc$cell %in% "CD4 T cell"] = "T cell"
pbmc$cell[pbmc$cell %in% "CD8 T cell"] = "T cell"
pbmc$cell[pbmc$cell %in% "Neutrophil"] = "Granulocytes"
pbmc$cell[pbmc$cell %in% "HSC"] = "Mast cell"
pbmc = subset(pbmc, cell != "Unsure")
save(pbmc,file="E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
