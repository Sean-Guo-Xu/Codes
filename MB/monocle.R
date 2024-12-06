
library(monocle)
library(Seurat)

setwd("E:\\brain")
load("UGCP_GCPP_GCP_GCP_NSC_UGC_cds.Rdata")
load("4type.Rdata")
pbmc = subset(pbmc, cell %in% c( "GCP_progenitor","UBC_progenitor","Neuron Stem cell"))
#cds = importCDS(p,import_all = T)
pbmc@active.ident = as.factor(pbmc$cell)
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
rm(pbmc,expr_matrix)
gc()
cds <- estimateSizeFactors(cds)
gc()
cds <- estimateDispersions(cds)
gc()
load("real_time_gene.Rdata")
allresult = allresult[allresult$para.pv < 0.0001,]
cds <- setOrderingFilter(cds, allresult$gene)  
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5) 
plot_cell_trajectory(cds,color_by=c("cell"), size=1,show_backbone=TRUE) 
ggsave("UCP_GCP_NSC_GCPP_cell.png",width = 6,height = 5)
plot_cell_trajectory(cds,color_by=c("Pseudotime"), size=1,show_backbone=TRUE) 
ggsave("UCP_GCP_PC_NSC_time.png",width = 6,height = 5)
save(cds,file = "UCP_GCP_PC_NSC_cds.Rdata")

########################

library(monocle)
setwd("E:\\brain")
load("4type.Rdata")
pbmc = subset(pbmc,cell != "Middle Brain Neuron" )
pbmc = subset(pbmc,cell != "choroid plexus cells" )
pbmc = subset(pbmc, cell %in% c("Neuron Stem cell" ,"GCP_progenitor","UBC_progenitor"))
#cds = importCDS(p,import_all = T)
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
rm(pbmc,expr_matrix)
gc()
cds <- estimateSizeFactors(cds)
gc()
cds <- estimateDispersions(cds)
load("real_time_gene.Rdata")
diff = differentialGeneTest(cds = cds,fullModelFormulaStr = "~cell",cores = 4)
diff = diff[diff$pval<0.05,]
diff = diff[diff$qval<0.05,]
gene = rownames(diff)[diff$use_for_ordering %in% "TRUE"]
cds <- setOrderingFilter(cds, gene)  
gc()
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5) 
plot_cell_trajectory(cds,color_by=c("cell"), size=1,show_backbone=TRUE) 
ggsave("GCP_GCPP_UCB_UCBP_NSC_cell.png",width = 5,height = 5)
plot_cell_trajectory(cds,color_by=c("Pseudotime"), size=1,show_backbone=TRUE) 
ggsave("GCP_GCPP_UCB_UCBP_NSC_time.png",width = 5,height = 5)
save(cds,file = "UCP_GCP_NSC_cds.Rdata")

data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds3 <- new_cell_data_set(data,
                          cell_metadata = cell_metadata,
                          gene_metadata = gene_annotation)
int.embed=cds@reducedDimS
int.embed = t(int.embed)
colnames(int.embed) = c("UMAP_1","UMAP_2")
cds3@int_colData$reducedDims$UMAP <- int.embed
cds3 <- cluster_cells(cds3)
cds3 <- learn_graph(cds3)
cds3 = order_cells(cds3)
plot_cells(cds3, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
pseudotime <- pseudotime(cds3, reduction_method = 'UMAP')
pbmc$pseudotime = pseudotime
cds$pseudotime = pseudotime

save(pbmc,file="GCP_GCPP_UCB_UCBP_NSC.Rdata")
save(cds,file="GCP_GCPP_UCB_UCBP_NSC_cds.Rdata")

load("GCP_GCPP_UCB_UCBP_NSC_cds.Rdata")
source("ordercell.R")
cds=ordercell(cds)
plot_genes_branched_pseudotime(cds["CD14",],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)
branchTest(cds)
cds = cds[,cds$cell %in% c("Neuron Stem cell" ,"GCP_progenitor","UBC_progenitor")]
all.vars(terms(as.formula("~sm.ns(Pseudotime, df = 3)*Branch")))
cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = NULL, 
                                     branch_point = 1, branch_labels = NULL )
?buildBranchCellDataSet

#########抽样
library(monocle)
library(Seurat)

setwd("E:\\brain")
load("4type.Rdata")
pbmc = subset(pbmc, cell %in% c( "GCP","GCP_progenitor","UBC_progenitor","Neuron Stem cell","Unipolar Brush cell"))
pbmc$cell
name = names(table(pbmc$cell))
num = as.integer(table(pbmc$cell)/2)
pbmc$index = colnames(pbmc)
for(i in 1:length(name)){

cellname = sample(colnames(pbmc)[pbmc$cell %in% name[i]],num[i])
pbmc = subset(pbmc,index %in% cellname,invert = T)
}
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
rm(pbmc,expr_matrix)
gc()
cds <- estimateSizeFactors(cds)
gc()
cds <- estimateDispersions(cds)
gc()
load("real_time_gene.Rdata")
allresult = allresult[allresult$para.pv < 0.0001,]
cds <- setOrderingFilter(cds, allresult$gene)  
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
save(cds,file="cds2.Rdata")
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3) 
plot_cell_trajectory(cds,color_by=c("cell"), size=1,show_backbone=TRUE) 
ggsave("5type_cell.png",width = 6,height = 5)
plot_cell_trajectory(cds,color_by=c("Pseudotime"), size=1,show_backbone=TRUE) 
ggsave("5type_time.png",width = 6,height = 5)

setwd("E:\\brain")
load("cds2.Rdata")
BEAM_res <- BEAM(cds[allresult$gene,], branch_point = 1, cores = 4, progenitor_method = "duplicate") 
load("5type_BEAM.Rdata")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
library(tidyverse)
BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
pdf("heatplot.pdf",height = 10,width = 10)
p=plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 1, 
                                 num_clusters = 4, show_rownames = T, return_heatmap = T)
dev.off()
