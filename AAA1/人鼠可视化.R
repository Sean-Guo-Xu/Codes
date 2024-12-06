library(Seurat)
setwd("E:\\AAA_scRNA")
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
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"), "#90d1a4", "#FB9A99" ,"Artery A")
load("E:\\AAA_scRNA\\6Normal.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"),  "#90d1a4", "#FB9A99" ,"Artery N")

load("E:\\AAA_scRNA\\Blood\\A_blood.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"),  "#90d1a4", "#FB9A99" ,"Blood A")
load("E:\\AAA_scRNA\\Blood\\N_blood.Rdata")
myfeatureplot(pbmc,c("EIF2B2","MLH3","NEK9"),  "#90d1a4", "#FB9A99" ,"Blood N")


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
  VlnPlot(pbmc,features = i,group.by = "cell",split.by ="sample",cols = c( "#90d1a4", "#FB9A99" ))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")
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
  VlnPlot(pbmc,features = i,group.by = "cell",split.by ="sample",cols = c( "#90d1a4", "#FB9A99" ))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")
  ggsave(paste(i,"blood.tiff"),width = 5,height = 5)                                                                                                                      
}
#########################

NK = c("KLRD1","NCAM1","NKG7") #7
T = c("CD3D","CD3E","CD3G") #0 17 26
Mast = c("KIT","CPA3","TPSB2") #16
Endothelial = c("PECAM1","PLVAP","PTPRB") #24 23 3
Fibroblast = c("DCN","COL1A1","COL1A2") #12 10 13
SMC = c("MYH11","ACTA2","MYL9") #4 25
Monocyte = c("FCN1","APOBEC3A","CD163","CD68","CD14") #2 5 8 9 19 22
B = c("CD79A","MZB1","MS4A1") #1 6 11 18 20 21 27
Proliferation = c("MKI67","TOP2A","CDC20") #15

load("E:\\AAA_scRNA\\8AAA.Rdata")
pbmc$cell[pbmc$cell %in% c("Monocyte","Macrophage")] = "Myeloid"
pbmc$cell[pbmc$cell %in% c("NK cell","T cell")] = "NK/T"
markers = c(T,B,Mast,Endothelial,Fibroblast,SMC,Monocyte,Proliferation)
cluster = c(rep("NK/T",3),rep("B cell",3),rep("Mast cell",3)
            ,rep("EC",3),rep("Fibroblast",3),rep("SMC",3),rep("Myeloid",5),rep("Proliferation",3))
library(ggplot2)
DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#90d1a4", "#FB9A99"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("Human_aaa_all_dotplot.tiff",width = 12,height = 4)

load("E:\\AAA_scRNA\\6Normal.Rdata")
pbmc$cell[pbmc$cell %in% c("Monocyte","Macrophage")] = "Myeloid"
pbmc$cell[pbmc$cell %in% c("NK cell","T cell")] = "NK/T"
markers = c(T,B,Mast,Endothelial,Fibroblast,SMC,Monocyte)
cluster = c(rep("NK/T",3),rep("B cell",3),rep("Mast cell",3)
            ,rep("EC",3),rep("Fibroblast",3),rep("SMC",3),rep("Myeloid",5))
library(ggplot2)
DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#90d1a4", "#FB9A99"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("Human_normal_all_dotplot.tiff",width = 14,height = 5)
#################mouse dimplot#############

t <- c("Ccl5", "Tcf7", "Icos", "Il7r")
ec<- c("Cdh5","Pecam1","Tie1")
fb <- c("Pdgfra","Col1a1","Lum")
smc <- c("Acta2","Myh11","Cnn1")
myeloid <- c("Adgre1","Cd68","Cd14","Ccr2","Lgals3")
b <- c("Cd79a","Cd79b","Ly6d")
granulocyte <- c("S100a8","S100a9","Alox5", "Cd53")
Proliferation = c("Mki67","Top2a","Cdc20") 
markers = c(t,ec,fb,smc ,myeloid,b,granulocyte,Proliferation)
load("E:\\AAA_scRNA\\mouse_aaa_all.Rdata")

cluster = c(rep("NK/T cell",4),rep("EC",3),rep("Fibroblast",3)
            ,rep("SMC",3),rep("Myeloid",5),rep("B cell",3),rep("Granulocyte",4),rep("Proliferation",3))
DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#90d1a4", "#FB9A99"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("mouse_aaa_all_dotplot.tiff",width = 15,height = 5)

load("E:\\AAA_scRNA\\mouse_normal_all.Rdata")
markers = c(t,ec,fb,smc ,myeloid,b)
cluster = c(rep("NK/T cell",4),rep("EC",3),rep("Fibroblast",3)
            ,rep("SMC",3),rep("Myeloid",5),rep("B cell",3))

DotPlot(pbmc,features = split(markers,cluster),group.by = "cell",cols = c("#90d1a4", "#FB9A99"))+RotatedAxis()+  theme(
  panel.border = element_rect(color = "black"),
  panel.spacing = unit(1, "mm"),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
ggsave("mouse_normal_all_dotplot.tiff",width = 14,height = 5)

pbmc$cell[pbmc$cell %in% "CD4 T cell"] = "T cell"
pbmc$cell[pbmc$cell %in% "CD8 T cell"] = "T cell"
pbmc$cell[pbmc$cell %in% "Neutrophil"] = "Granulocytes"
pbmc$cell[pbmc$cell %in% "HSC"] = "Mast cell"
pbmc = subset(pbmc, cell != "Unsure")
save(pbmc,file="E:\\AAA_scRNA\\Blood\\N_blood.Rdata")

#############heatmap#############
library(Seurat)
setwd("E:\\AAA_scRNA")
library(ggplot2)
library(ComplexHeatmap)
load("mouse_human.Rdata")
pbmc = integrated_MPC
pbmc$cell[pbmc$cell %in% c("Foamy/TREM2hi-Mø","Foamy-Mø/Trem2hi")] = "Foamy-TREM2hi"
pbmc$cell[pbmc$cell %in% c("Foamy/SPP1hi-Mø","Foamy-Mø/Spp1hi")] = "Foamy-SPP1hi"
pbmc$cell[pbmc$cell %in% c("TRM","TRM/FLOR2hi-Mø")] = "TRM"
pbmc$cell[pbmc$cell %in% c("Pro-Inf","ProInf-Mø")] = "ProInf"
pbmc@active.ident  = as.factor(pbmc$cell)
markers=FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
markers = markers[order(markers$avg_log2FC,decreasing = T),]
top10 = NULL
for (i in unique(pbmc$cell)) {
  part= markers[markers$cluster %in% i,]
  part = part[part$avg_log2FC>0,]
  if(nrow(part) >= 10){
    top10 = rbind(top10,part[1:10,])}
  else{
    top10 = rbind(top10,part)
  }
}
unique(pbmc$lab[order(pbmc$lab)])
pbmc$lab = factor(pbmc$lab,levels =unique(pbmc$lab[order(pbmc$lab)]))
pbmc$lab = paste(pbmc$cell,pbmc$species,sep = "_")

library(gcookbook)
col = rep(c("#f57245","#e6f49a","#A6CEE3","#fbe18c","#90d1a4","#33A02C","#FB9A99","#feac60","#B2DF8A"),2)
cellname =  c("Mono","Foamy-TREM2hi","TRM","Inf-Mø","cDC1","cDC2","IFNIC-Mø","ProInf",
               "Foamy-SPP1hi")
names(col) = c(paste(cellname,"human",sep = "_"),paste(cellname,"mouse",sep = "_"))

DoHeatmap(pbmc,features = top10$gene,group.by = "lab",group.colors = col,size=4)+scale_fill_gradientn(colors = c("#90d1a4","white","#FB9A99"))
ggsave("top10_heatmap.tiff",width = 8,height = 10)
DoHeatmap(subset(pbmc,species =="mouse"),features = top10$gene,group.by = "cell",group.colors = col,size=4)+scale_fill_gradientn(colors = c("#90d1a4","white","#FB9A99"))
ggsave("top10_mouse.tiff",width = 6,height = 10)
Mono	S100A8	#f57245
Mono	VCAN	#f57245
Foamy/TREM2hi-Mø	C1QB	#e6f49a
Foamy/TREM2hi-Mø	APOC1	#e6f49a
Foamy/TREM2hi-Mø	TREM2	#e6f49a
TRM	F13A1	#A6CEE3
TRM	LYVE1	#A6CEE3
Inf-Mø	CD74	#fbe18c
Inf-Mø	CCL4	#fbe18c
Inf-Mø	CCL3	#fbe18c
cDC1	CLEC9A	#90d1a4
cDC2	CLEC10A	#33A02C
IFNIC-Mø	ISG15	#FB9A99
IFNIC-Mø	IFIT2	#FB9A99
ProInf-Mø	DUSP2	#feac60
ProInf-Mø	NFKBIA	#feac60
ProInf-Mø	THBS1	#feac60
Foamy/SPP1hi-Mø	SPP1	#B2DF8A
Foamy/SPP1hi-Mø	CTSB	#B2DF8A
Foamy/SPP1hi-Mø	SDC2	#B2DF8A
