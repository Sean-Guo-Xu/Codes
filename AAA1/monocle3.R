setwd("D:\\monocle")
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
load("mouse_human.Rdata")
#pbmc  =integrated_MPC
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
library(SeuratWrappers)
cds <- as.cell_data_set(pbmc)
cds@clusters$UMAP$clusters <- Idents(pbmc)[rownames(colData(cds))]
cds@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(cds)))), levels = 1)
names(cds@clusters$UMAP$partitions) <- rownames(colData(cds))
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(cds)
rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds2 <- learn_graph(cds, use_partition = F)
cds2 <- order_cells(cds2)
plot_cells(cds2, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
DimPlot(pbmc,group.by = "cell",label=T)
save(cds2,file="human_mouse_cds.Rdata")
