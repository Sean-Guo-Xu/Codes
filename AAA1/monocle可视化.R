setwd("F:\\AAA_scRNA\\monocle")
library(Seurat)
library(monocle3)
library(tidyverse)
library(SCENIC)
library(patchwork)
library(gcookbook)
load("human_mouse_cds.Rdata")
load("F:\\AAA_scRNA\\mouse_human.Rdata")
save(pbmc,file="E:\\AAA\\mouse_human.Rdata")
umap = pbmc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = pbmc@meta.data$cell) 
p=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
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
p +scale_colour_gradient(low = "#90d1a4", high ="#FB9A99" )+
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

ggsave("human_mouse_tra.tiff",width = 6,height = 5)
DimPlot(pbmc,group.by = "cell",label=T)
load("human_mouse_cds.Rdata")
Track_genes <- graph_test(cds2, neighbor_graph="principal_graph", cores=4)
Track_genes <-Track_genes[Track_genes$p_value<0.001,]
Track_genes <-Track_genes[!is.na(Track_genes$p_value),]
Track_genes <- Track_genes[Track_genes$gene_short_name %in% marker$V2,]
gene = read.table("gene.txt",sep="\t",comment.char="", encoding = "UTF-8")
load("time_marker_genes.Rdata")
color = c(unique(gene$V3))
names(color) = unique(gene$V1)
plot_genes_in_pseudotime(cds[gene$V2,], color_cells_by="cell", 
                              min_expr=0.5, ncol = 4)+scale_color_manual(values = color)
 ggsave("time_cell.tiff",height = 10,width = 16)
 plot_genes_in_pseudotime(cds[gene$V2,], color_cells_by="species", 
                          min_expr=0.5, ncol = 4)+scale_color_manual(values = c( "#90d1a4","#FB9A99"))
 ggsave("time_species.tiff",height = 10,width = 16)
 cds2[as.character(Track_genes$gene_short_name),]
pbmc = integrated_MPC
cds2@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- pbmc@reductions$pca@feature.loadings
cds2@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- pbmc@reductions$pca@stdev
marker = read.table("gene.txt",sep="\t",encoding = "UTF-8")
save(Track_genes,file="time_marker_genes.Rdata")


cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds2, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
ggsave("Genes_Module.pdf", plot = p, width = 8, height = 8)

color=c("#b71f48","#e45d52","#f57245","#feac60","#fbe18c","#e6f49a","#90d1a4","#4fa5b2","#4f98c6","#5371b3")
color=c("grey","#b71f48","#f57245","#feac60","#FB9A99","#fbe18c","#e6f49a","#B2DF8A","#90d1a4","#33A02C","#A6CEE3","#4f98c6","#5371b3",
        "#1F78B4","#FDBF6F","#CAB2D6", "#6A3D9A","#35978F")