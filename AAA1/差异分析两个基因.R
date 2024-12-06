setwd("F:\\AAA_scRNA")
library(ggplot2)
library(Seurat)
library("Nebulosa")
mygenes <- c("CSK","PLEKHJ1","AMH","LDAH","NEK9","PSRC1","SCAPER","HLA-DRB1","FAM66C")  # 示例基因
load("7AAA.Rdata")
myfeatureplot<-function(pbmc,allgene,lowcolor,midcolor,highcolor,tissue){
  for (gene in allgene) {
    conver = function(x){
      x[x %in% 0] = NA
      return(x)
    }
    data = data.frame(FetchData(pbmc, vars=gene,slot="data"), pbmc@reductions$umap@cell.embeddings)
    data=apply(data, 2, conver)
    data = as.data.frame(data)
    colnames(data)=c("gene","UMAP_1","UMAP_2")
    ggplot() +  geom_point(data=data[is.na(data[,1]),],aes(UMAP_1, UMAP_2, color=gene), size=1.5) +geom_point(data=data[!is.na(data[,1]),],aes(UMAP_1, UMAP_2, color=gene), size=1.5)+scale_colour_gradient2(low = lowcolor,mid = midcolor,high = highcolor,midpoint = mean(data[!is.na(data[,1]),1]))  + theme_light(base_size = 15)+labs(title = "")+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+labs(color = "Expression")+
      theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank(),axis.text.x = element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) + xlab("")+ylab("")+ggtitle(gene)
    ggsave(paste(gene,tissue,"featureplot.tiff"),width = 5.5,height = 5)
  }}
library(ggsci)
load("7AAA.Rdata")
myfeatureplot(pbmc,mygenes,"#62B197","#fbe18c","#E18E6D","Aorta")
load("3Normal.Rdata")
myfeatureplot(pbmc,mygenes,"#62B197","#fbe18c","#E18E6D","Aorta_N")
load("blood/A_blood.Rdata")
myfeatureplot(pbmc,mygenes,"#62B197","#fbe18c","#E18E6D","PBMC_A")
load("blood/N_blood.Rdata")
DimPlot(pbmc,group.by = "cell")+scale_color_d3()+  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+ggtitle("")
ggsave("nblood.tiff",width = 6,height = 5)
myfeatureplot(pbmc,mygenes,"#62B197","#fbe18c","#E18E6D","PBMC_N")
mygenes <- c("CSK","PLEKHJ1","AMH","LDAH","NEK9","PSRC1","SCAPER","HLA-DRB1","FAM66C")  # 示例基因
plot_density(aaa, c("CSK","PLEKHJ1","AMH","LDAH"))
ggsave("",)
plot_density(aaa,mygenes)
ggsave("9_genes.tiff",width = 12,height = 10)
load("human_aaa_myeloid.Rdata")
plot_density(pbmc, c("CD44","SPP1"))
ggsave("human_myeloid_all.tiff",width = 12,height = 5)

load("mouse_aaa_all.Rdata")
plot_density(pbmc, c("Cd44","Spp1"))
ggsave("mouse_aaa_all.tiff",width = 12,height = 5)
load("mouse_aaa_myeloid.Rdata")
plot_density(pbmc, c("Cd44","Spp1"))
ggsave("mouse_aaa_myeloid.tiff",width = 12,height = 5)
pbmc = subset(pbmc,cell != "Unsure")

load("7AAA.Rdata")
plot_density(pbmc,mygenes)
ggsave("9_genes_aaa_Artery.tiff",width = 12,height = 10)
unique(pbmc$cell)
aaa=pbmc
aaa$sample = "AAA"

library(ggpubr)
load("3Normal.Rdata")
unique(pbmc$orig.ident)
DimPlot(pbmc,label = T,reduction = "UMAP")
mygenes[3] ="CD14"
plot_density(pbmc,mygenes,reduction = "UMAP")
ggsave("9_genes_normal_Artery.tiff",width = 12,height = 10)
pbmc$cell[pbmc$cell %in% c("Macrophage","Monocyte")] = "Myeloid"
pbmc@assays$METABOLISM=NULL
celltype = unique(intersect(pbmc$cell,aaa$cell))
pbmc=merge(aaa,pbmc)
pbmc = NormalizeData(pbmc)
pbmc=ScaleData(pbmc)
pbmc$all="all"
deg_artery = FindMarkers(pbmc,ident.1 ="AAA" ,ident.2 ="Normal",group.by = "sample",  logfc.threshold = 0, min.pct = 0,features = mygenes)

VolcanoPlot=function(dif, log2FC=0.5, padj=0.05, 
                     label.symbols=NULL, label.max=30,
                     cols=c("#62B197", "#E18E6D"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[dif$log2FoldChange > log2FC & dif$padj <padj,]$threshold="up";
  if(nrow(dif[dif$log2FoldChange < (-log2FC) & dif$padj < padj,]) !=0 ){
  dif[dif$log2FoldChange < (-log2FC) & dif$padj < padj,]$threshold="down";}
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=4) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("DEG:", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print((dd_text))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5, 
                       colour="black",alpha=1)
}
rownames(deg_blood) = paste(rownames(deg_blood),"(B)",sep = "")
rownames(deg_artery) = paste(rownames(deg_artery),"(A)",sep = "")
deg_all = rbind(deg_artery,deg_blood)
dif=data.frame(
  symbol=rownames(deg_all),
  log2FoldChange=deg_all$avg_log2FC,
  padj=deg_all$p_val_adj
)
dif$padj[dif$padj == 0] = 1e-250
VolcanoPlot(dif, log2FC = 0.1,padj=0.05, title="AAAs vs Controls", label.max = 50)
ggsave("volcano.tiff",width = 6,height = 5.5)
pbmc$all = "all"
library(ggpubr)
for (i in mygenes) {
  VlnPlot(pbmc,slot = "data",features = i,group.by = "sample",raster = F,cols = c("#FB9A99","#90d1a4"),pt.size = 0)+xlab("")+theme(legend.position = 'none',)
  ggsave(paste(i,"Artery_all.tiff"),width = 2,height = 4) 
  #VlnPlot(pbmc,slot = "data",features = i,group.by = "cell",raster = F,split.by ="sample",cols = c("#FB9A99","#90d1a4"))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")
 # ggsave(paste(i,"Artery_cell.tiff"),width = 6,height = 4)                                                                                                                      
}
marker = list()
marker[["total"]] = FindMarkers(pbmc, group.by = "sample",ident.1 = "AAA",ident.2 = "Normal" ,features = mygenes,min.pct=0,logfc.threshold=0)
for (i in unique(celltype)) {
  part = subset(pbmc,cell == i)
  marker[[i]] =  FindMarkers(part, group.by = "sample",ident.1 = "AAA",ident.2 = "Normal" ,features = mygenes,min.pct=0,logfc.threshold=0)
}
save(marker,file="Arery_diff.Rdata")
library(ggsci)
pbmc = NormalizeData(pbmc)
pbmc = ScaleData(pbmc)
DoHeatmap(pbmc,features =mygenes,group.by = "sample" ,split.by="cell", group.colors =pal_d3()(length(unique(pbmc$cell))),size = 4,angle = 45)+scale_fill_gradient(low="#6BB7CA",high="#E07B54")
ggsave("Aorta_heatmap.tiff",width = 5,height = 5)


library(batchelor)
library(SingleCellExperiment)
load("blood/A_blood.Rdata")
plot_density(pbmc,mygenes)
ggsave("9_genes_aaa_blood.tiff",width = 12,height = 10)
aaa = pbmc  
aaa$sample = "AAA"
load("blood/N_blood.Rdata")
plot_density(pbmc,mygenes)
ggsave("9_genes_normal_blood.tiff",width = 12,height = 10)
pbmc$sample = "Normal"
pbmc = merge(aaa,pbmc)
pbmc=NormalizeData(pbmc)
pbmc = ScaleData(pbmc)
library(Seurat)
library(ggpubr)
deg_blood = FindMarkers(pbmc,ident.1 ="AAA" ,ident.2 ="Normal",group.by = "sample",  logfc.threshold = 0, min.pct = 0,features = mygenes)
RidgePlot(pbmc,features = "PLEK",group.by = "cell")
for (i in mygenes) {
  VlnPlot(pbmc,slot = "data",features = i,pt.size = 0,group.by = "sample",raster = F,cols = c("#FB9A99","#90d1a4"))+xlab("")+theme(legend.position = 'none')
  ggsave(paste(i,"Blood_all.tiff"),width = 2,height = 4) 
  #VlnPlot(pbmc,slot = "data",features = i,group.by = "cell",raster = F,split.by ="sample",cols = c("#FB9A99","#90d1a4"))+stat_compare_means( aes(label = ..p.signif..),  method = "wilcox.test")+xlab("")
  #ggsave(paste(i,"Blood_cell.png"),width = 6,height = 4)                                                                                                                      
}

marker = list()
marker[["total"]] = FindMarkers(pbmc, group.by = "sample",ident.1 = "AAA",ident.2 = "Normal" ,features = mygenes,min.pct=0,logfc.threshold=0)
for (i in unique(celltype)) {
  part = subset(pbmc,cell == i)
  marker[[i]] =  FindMarkers(part, group.by = "sample",ident.1 = "AAA",ident.2 = "Normal" ,features = mygenes,min.pct=0,logfc.threshold=0)
}
save(marker,file="Blood_diff.Rdata")
library(ggsci)
pbmc = NormalizeData(pbmc)
pbmc = ScaleData(pbmc)
DoHeatmap(pbmc,features =mygenes,group.by = "sample" , group.colors =pal_d3()(length(unique(pbmc$cell))),size = 4,angle = 45)+scale_fill_gradient(low="#6BB7CA",high="#E07B54")
ggsave("Blood_heatmap.tiff",width = 5,height = 5)

library(scRNAtoolVis)
load("blood/A_blood.Rdata")
aaa =pbmc
aaa$sample = "AAA"
load("blood/N_blood.Rdata")
pbmc$sample ="Normal"
pbmc = merge(aaa,pbmc)
pbmc@active.ident = factor(pbmc$cell)
allmarker = NULL
for (i in unique(pbmc$cell)) {
marker <- FindMarkers(subset(pbmc,cell == i), group.by ="sample",only.pos = FALSE,ident.1 = "AAA",ident.2 = "Normal",
                               min.pct = 0.2,
                               logfc.threshold = 0)
marker$cluster = i 
marker$gene = rownames(marker)
allmarker = rbind(allmarker,marker)
}
allmarker$cluster = as.factor(allmarker$cluster)
allmarker = allmarker[allmarker$pct.1 != 0,]
allmarker = allmarker[allmarker$pct.2 != 0,]
write.csv(allmarker,"PBMC_diff.csv",row.names = F)
source("cellvolplot.R")
cellvolplot(diffData = allmarker,log2FC.cutoff = 0.5,aesCol = c("#62B197", "#E18E6D"),  myMarkers = mygenes, tile.col=pal_d3()(length(unique(allmarker$cluster))),labelsize = 3)
ggsave("cell_vol.tiff",width = 10,height = 6)

load("7AAA.Rdata")
aaa =pbmc
aaa$sample = "AAA"
load("3Normal.Rdata")
celltype = intersect(pbmc$cell,aaa$cell)
pbmc$sample ="Normal"
pbmc = merge(aaa,pbmc)
pbmc@active.ident = factor(pbmc$cell)
pbmc = JoinLayers(pbmc)
pbmc = subset(pbmc, cell %in% celltype)
allmarker = NULL
for (i in unique(pbmc$cell)) {
  marker <- FindMarkers(subset(pbmc,cell == i), group.by ="sample",only.pos = FALSE,ident.1 = "AAA",ident.2 = "Normal",
                        min.pct = 0.2,
                        logfc.threshold = 0)
  marker$cluster = i 
  marker$gene = rownames(marker)
  allmarker = rbind(allmarker,marker)
}
allmarker$cluster = as.factor(allmarker$cluster)
allmarker = allmarker[allmarker$pct.1 != 0,]
allmarker = allmarker[allmarker$pct.2 != 0,]
write.csv(allmarker,"aorta_diff.csv",row.names = F)
source("cellvolplot.R")
library(ggsci)
cellvolplot(diffData = allmarker,log2FC.cutoff = 0.5,aesCol = c("#62B197", "#E18E6D"),  myMarkers = mygenes, tile.col=pal_d3()(length(unique(allmarker$cluster))),labelsize = 3)
ggsave("cell_vol.tiff",width = 10,height = 6)
