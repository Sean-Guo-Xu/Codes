library(Seurat)
library(harmony)
library(snow)
library(monocle)
load("E:\\AAA_scRNA\\blood\\A_blood.Rdata")
blood = subset(pbmc,cell==c("Monocyte"))
blood$cell = "pbmc_Monocyte"
load("E:\\AAA_scRNA\\A_art_Macro.Rdata")
artery = pbmc
pbmc=merge(blood,artery)
bcounts <- GetAssayData(blood, assay = "RNA")
acounts<- GetAssayData(artery, assay = "RNA")
name = rownames(bcounts)[rownames(bcounts) %in% rownames(acounts)]
pbmc = subset(pbmc,features = name)

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)

pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
pbmc = RunHarmony(pbmc,"orig.ident", plot_convergence = TRUE)#耗时1min

pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  identity()


DimPlot(pbmc,label = T,group.by = c("cell","seurat_clusters"))
pbmc$celltype = pbmc$cell
pbmc$celltype[pbmc$celltype %in% "pbmc_Monocyte"] = "Monocyte"


expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
p_data$celltype <- pbmc$celltype
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
save(cds,pbmc,file="MM cds.Rdata")
setwd("E:\\AAA_scRNA")
load("MM cds.Rdata")
pbmc$celltype=pbmc$cell
pbmc$celltype[pbmc$celltype %in% "pbmc_Monocyte"] = "Monocyte"
pbmc@active.ident = as.factor(pbmc$celltype)
memory.limit(size = 10000000)
marker=FindAllMarkers(pbmc,logfc.threshold = 1,min.pct = 0.25)
cds <- setOrderingFilter(cds, unique(marker$gene))
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
save(cds,pbmc,file="MM cds.Rdata")
plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE) 
cds <- orderCells(cds)
cds$celltype[!(cds$celltype %in% "AAA")] = "Blood" 

cds$celltype = pbmc$celltype
plot_ordering_genes(cds)
plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE) 

save(cds,pbmc,file="E:\\AAA_scRNA\\MM trajectory.Rdata")
BEAM_res <- BEAM(cds, branch_point = 1, cores = 4) 
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf("point heatmap.pdf",width = 8,height = 15)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 6,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T)
dev.off()
Time_diff <- differentialGeneTest(cds, cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff = Time_diff[Time_diff$qval<0.001,]
Time_diff = Time_diff[order(Time_diff$qval),]


plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points = F,markers="AKT",markers_linear = T) 
cds$gene =c( FetchData(pbmc,"rna_CD163"))
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
load("E:\\AAA_scRNA\\MM.Rdata")
FeaturePlot(pbmc,features = c("CD163","FCGR3A"))
pbmc = subset(pbmc,seurat_clusters != 7)
pbmc$pseudotime = sc$slingPseudotime_2
pbmc = subset(pbmc,pseudotime %in% pbmc$pseudotime[!is.na(pbmc$pseudotime)] )
sc$seurat_clusters = droplevels(sc$seurat_clusters)
save(pbmc,cds,file="E:\\AAA_scRNA\\MM.Rdata")

load("E:\\AAA_scRNA\\MM trajectory.Rdata")

cds$gene=c( FetchData(pbmc,vars="rna_UGCG"))[[1]]
plot_cell_trajectory(cds,color_by = "gene")
subcds = cds["UGCG",]
plot_genes_in_pseudotime(subcds,color_by = "Pseudotime")
plot_cell_trajectory(cds,color_by = "sample")
plot_cell_trajectory(cds,color_by="State",show_branch_points = F,markers="PPARG",markers_linear = T,show_tree = T)
pbmc$state = as.numeric(cds$State)
library(ggpubr)
library(ggplot2)
VlnPlot(pbmc,group.by = "state",features = "HLA-DQB2")+geom_signif(comparisons = com,test = wilcox.test)
com = list(c(1,2),c(2,3),c(1,3))
###############
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(RColorBrewer)
pbmc$cell = as.factor(pbmc$cell)
pbmc$seurat_clusters = droplevels(pbmc$seurat_clusters)
pbmc <- as.SingleCellExperiment(pbmc)
pbmc <- slingshot(pbmc, clusterLabels = "seurat_clusters", reducedDim = "UMAP",end.clus=c(1,12,7))

sc$cell = as.factor(sc$cell)
sc$seurat_clusters = droplevels(sc$seurat_clusters)
sc <- as.SingleCellExperiment(sc)
sc <- slingshot(sc, clusterLabels = "seurat_clusters", reducedDim = "UMAP")

line <- getLineages(sc, 
                    clusterLabels = "seurat_clusters", 
                    #start.clus = 'Pi16',#可指定起始细胞簇
                   # end.clus=c(1,12,7),#可指定终点细胞簇
                    reducedDim = "UMAP"
                    )

cor =as.character(sc$cell)
cor[cor %in% "Monocyte"] = "#E64B35FF"
cor[cor %in% "Macrophage"] = "#4DBBD5FF"

pdf("E:\\Tra.pdf",width = 5,height = 5)
plot(reducedDims(sc)$UMAP,col = cor,pch=16)
FeaturePlot(pbmc,features = "pseudotime")
lines(SlingshotDataSet(line), lwd=2,col = 'black',type = 'lineages')
dev.off()
see=SlingshotDataSet(line)

