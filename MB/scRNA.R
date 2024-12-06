library(Seurat)
library(monocle3)
setwd("E:\\brain")
pbmc = readRDS("GSM5952337_Fetal_cerebellum_final_cluster_1.rds")
pbmc = subset(pbmc,seurat_clusters %in% c(35,29,21,17,28,38))
save(pbmc,file="immune_normal.Rdata")
   library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =pbmc@assays$RNA@data, ref = ref, labels = ref$label.main)
pbmc$SingleR = cellpred$labels
remotes::install_version("matrixStats", version="1.1.0")
BiocManager::install("SingleR", update = TRUE, ask = FALSE)
pbmc$singleR=cellpred$labels
DimPlot(pbmc,label=T,group.by = c("SingleR","celltype","seurat_clusters"))


ggsave("cluster.png",width = 7.5,height = 5) 
unique(pbmc$batch)
pbmc@active.ident =  pbmc$celltype
marker = FindAllMarkers(pbmc,min.pct = 0.2,logfc.threshold = 0.5,)
save(marker,file= "nromal_marker.Rdata")
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
BiocManager::install("SeuratWrappers")
library(SeuratWrappers)
cds <- as.cell_data_set(pbmc)
cds@clusters$UMAP$clusters <- Idents(pbmc)[rownames(colData(cds))]
cds@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(cds)))), levels = 1)
names(cds@clusters$UMAP$partitions) <- rownames(colData(cds))
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(cds)
rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds)
library(ggplot2)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
ggsave("pseudotime.png",width = 5.5,height = 5)
clusterlist = c(40,30,26,41,39,45,44,21,17,28,38,7,29) 
pbmc = subset(pbmc, seurat_clusters %in% clusterlist, invert=T)
DimPlot(pbmc,label=T,group.by = c("batch","seurat_clusters"))
save(pbmc,file = "E:\\neuron.Rdata")
load( "E:\\neuron.Rdata")
                                                                    load("neuron.Rdata")
load("marker.Rdata")
pbmc<-FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 2000)
pbmc<-RunPCA(pbmc,features = VariableFeatures(object=pbmc))
 pbmc<-FindNeighbors(pbmc)
pbmc<-FindClusters(pbmc,resolution = 1.5)
pbmc<-RunUMAP(pbmc,dims = 1:30)
DimPlot(pbmc,group.by=c("seurat_clusters"),label=T)
FeaturePlot(pbmc,features = c("SORCS3","PTPRK"))
pbmc@active.ident 
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
library(ggplot2)
cds <- cluster_cells(cds)
cds <- learn_graph(cds,close_loop = F)
 cds <- order_cells(cds)
 save(cds,file="E:\\monocle3_cds.Rdata")
 plot_cells(cds, color_cells_by = "seurat_clusters", label_cell_groups=F, label_leaves=FALSE, labels_per_group = T,label_branch_points=FALSE, graph_label_size=1.5)
ggsave("monocle3.png",width = 5,height = 5)
 write.table(marker,file="D:\\marker.txt",quote = F,sep = "\t",row.names = F)
load("E:\\neuron.Rdata")
pbmc$cell = as.numeric(pbmc$seurat_clusters)
pbmc$cell[pbmc$seurat_clusters %in% c(11,29,24)] ="Interneuron" #KIT LHX5 PAX2
#pbmc$cell[pbmc$seurat_clusters %in% c(26,35)] ="Molecular Layer Interneuron" #NXPH1 CXCL12
pbmc$cell[pbmc$seurat_clusters %in% c(28,25,22,13,18,26,27,31,32)] ="Neuron Stem cell" #SOX2 NES
pbmc$cell[pbmc$seurat_clusters %in% c(1,10,8,12,9)] ="Granule cell" ## GAP43 GRIK2 SEMA6A
pbmc$cell[pbmc$seurat_clusters %in% c(16)] ="Middle Brain Neuron" ##
pbmc$cell[pbmc$seurat_clusters %in% c(4,5,2,30)] ="GCP" # ATOH1
pbmc$cell[pbmc$seurat_clusters %in% c(30,32,21,15,17)] ="Purkinje cell" # RORA PCP4
pbmc$cell[pbmc$seurat_clusters %in% c(0,23,14,19,20)] ="Unipolar Brush cell" #OTX2 RSPO3 EOEMS 
pbmc$cell[pbmc$seurat_clusters %in% c(3,7)] ="GCP_progenitor" #RFC3
pbmc$cell[pbmc$seurat_clusters %in% c(20)] ="UBC_progenitor" #RFC3
#pbmc$cell[pbmc$seurat_clusters %in% c(10,29)] ="Neuron Stem cell" #SOX11 CTNNB1 HNRNPH1
pbmc$cell[pbmc$seurat_clusters %in% c(6)] ="Choroid Plexus Progenitor" #SOX11 CTNNB1 HNRNPH1
save(pbmc,file="E:\\neuron.Rdata")
pbmc = subset(pbmc,cell != "Middle Brain Neuron")
DimPlot(pbmc,group.by = "cell",label = T)
gene = c("KIT", "LHX5" ,"PAX2","NXPH1","SLC1A3","SOX2","NES","GAP43","GRIK2","MGP","ATOH1",
         "RORA","PCP4","OTX2","RSPO3","EOMES","RFC3","SOX11","CTNNB1","SPP1")
library(ggplot2)
DotPlot(pbmc,group.by = "cell",features = gene) +
  theme(axis.text.x = element_text(angle= 45 , vjust= 1 , hjust= 1 )) 
pbmc@active.ident = as.factor(pbmc$batch)
timemarker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc@active.ident 
save(marker,file = "marker.Rdata")
ggsave("my_anno.tiff",width = 6.5,height = 5)
####VEXTOR#####
VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings

source('vector.R')
PCA=vector.rankPCA(PCA)
# Define pixel
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)

# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

# Calculate Quantile Polarization (QP) score
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

# Get pixel's QP score
OUT=vector.gridValue(OUT,SHOW=TRUE)

# Find starting point
#OUT=vector.selectCenter(OUT)
OUT = vector.autoCenter(OUT)
# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
library(ggplot2)
ggsave("monocle.tiff",width = 5.5,height = 5)

########top10######
setwd("E:\\brain")
load("E:\\brain\\time_marker.Rdata")
load("neuron.Rdata")
timemarker = timemarker[order(timemarker$avg_log2FC,decreasing = T),]
top10 = NULL
for (i in unique(pbmc$batch)) {
  part = timemarker[timemarker$cluster %in% i,]
  part = part[1:10,]
  top10 = rbind(top10,part)
}
write.csv(top10,row.names = F,file = "time_top10.csv")
pbmc$batch= factor(pbmc$batch,levels = c("PCW8","PCW9","PCW12","PCW13","PCW14","PCW15","PCW16","PCW17"))
pbmc$cell[pbmc$cell %in% "choroid plexus progenitor cells"] = "choroid plexus cells"
save(pbmc,file="neuron.Rdata")
DimPlot(pbmc,group.by = "cell",label = T,split.by = "batch",ncol = 4)

pbmc = subset(pbmc,cell != "Middle Brain Neuron" )
pbmc = subset(pbmc,cell != "choroid plexus cells" )
save(pbmc,file="4type.Rdata")

library(ggsci)
library(ggplot2)
pbmc$batch= factor(pbmc$batch,levels = c("PCW8","PCW9","PCW12","PCW13","PCW14","PCW15","PCW16","PCW17"))
DimPlot(pbmc,split.by = "batch",label = F,ncol = 4,group.by = "cell",cols = color)+ggtitle("")
ggsave("8 weeks.png",width = 20,height = 8)
color = c("#FF0000","#E69F00",  "#56B4E9",  "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#8E44AD",
          "#1ABC9C", "#F39C12", "#2ECC71", "brown4")
DimPlot(pbmc,label = F,group.by = "cell",cols = color)+ggtitle("")+theme(legend.position = "None")
ggsave("All.png",width = 5,height = 5)
