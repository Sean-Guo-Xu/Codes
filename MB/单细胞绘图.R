library(ggsci)
library(ggplot2)
library(Seurat)
setwd("E:\\brain")
load("neuron.Rdata")
DimPlot(pbmc,group.by = c("cell","celltype"),label=T)
ggsave("umap.tiff",width = 14,height = 5)
DimPlot(pbmc,group.by = "cell",label = T,split.by = "batch",ncol = 4)
ggsave("batch_umap.tiff",width = 20,height = 8)
pbmc$orig.ident = pbmc$batch
unique(pbmc$orig.ident)
sample_table <- as.data.frame(table(pbmc$orig.ident,pbmc$cell))
names(sample_table) <- c("Samples","celltype","CellNumber")
ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_d3()(length(unique(pbmc$cell)))) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
ggsave("proportion.tiff",height = 6,width = 5)
gene = c("KIT", "LHX5" ,"PAX2","SOX2","NES","GAP43","GRIK2","SEMA6A","ATOH1",
         "RORA","PCP4","OTX2","RSPO3","EOMES","RFC3","SOX11","CTNNB1","SPP1")
DotPlot(pbmc,group.by = "cell",features = gene) +
  theme(axis.text.x = element_text(angle= 45 , vjust= 1 , hjust= 1 )) 
ggsave("dotplot_marker.tiff", width = 8, height = 5)
