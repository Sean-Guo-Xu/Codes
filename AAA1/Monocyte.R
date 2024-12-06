setwd("E:\\AAA_scRNA\\blood")
library(Seurat)
library(monocle3)
load("E:\\AAA_scRNA\\Blood\\A_blood.Rdata")
blood = pbmc
load("Macro.Rdata")
art = pbmc 
art$tissue = "Artery"
blood$tissue = "Blood"
pbmc=merge(art,blood)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc)) 
ElbowPlot(pbmc, ndims=20, reduction="pca")
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
DimPlot(pbmc,group.by = "cell")
pbmc = RunUMAP(pbmc,dims = 1:20)
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
cds <-learn_graph(cds)
plot_cells(cds, color_cells_by = "cell", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)

library(harmony)
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#耗时1min

#细胞聚类
###resolustion AAA 0.8 NORMAL 0.6#######
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  identity()
DimPlot(pbmc,label = T,group.by = "cell")

data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
cds <-learn_graph(cds)
plot_cells(cds, color_cells_by = "cell", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)
