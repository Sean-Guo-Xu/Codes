setwd("E:\\brain")
library(Seurat)
library(PseudotimeDE)
library(tibble)
library(dplyr)
library(scales)
library(RColorBrewer)
library(SingleCellExperiment)
library(irlba)
library(slingshot)
library(parallel)
library(BiocParallel)
library(scales)
library(paletteer) 
library(viridis)
source("persudoDE.R")
source("run.R")
load("neuron.Rdata")
load("time_marker.Rdata")

timemarker=timemarker[abs(timemarker$avg_log2FC)>1, ]

pbmc$time  = gsub("PCW","",pbmc$batch)
pbmc$time = as.integer(pbmc$time)
counts <- as.matrix(pbmc@assays$RNA@counts)
counts = counts[rownames(counts) %in% unique(timemarker$gene),]
counts = counts[401:500,]
cell_info <- data.frame(cell = colnames(counts),pseudotime = pbmc$time)

sub.tbl <- lapply(1:50, function(i) {
  sample_idx <- sample(nrow(cell_info), size = 0.8 * nrow(cell_info))
  see=cell_info[sample_idx, ]
})
rm(pbmc)
# Run PseudotimeDE
res <- run(gene.vec = rownames(counts),
                       ori.tbl = cell_info,
                       sub.tbl = sub.tbl,
                       mat = counts,
                       model = "nb",
                       mc.cores=1)
save(res,file="E:\\brain\\res5.Rdata")
  
setwd("E:\\brain")
allresult = NULL
for(i in 1:7)
{
  load(paste("res",i,".Rdata",sep = ""))
  res =as.data.frame(cbind(res$gene,res$fix.pv,res$emp.pv,res$para.pv,res$ad.pv,res$rank,res$test.statistics))
  colnames(res)=c("gene","fix.pv","emp.pv","para.pv","ad.pv","rank","test.statistics")
  allresult = rbind(allresult,res)
}
allresult[,2:7]=apply(allresult[,2:7], 2, as.numeric) 

library(dyno)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
setwd("E:\\brain")
load("neuron.Rdata")
pbmc = subset(pbmc,cell != "Middle Brain Neuron" )
pbmc = subset(pbmc,cell != "choroid plexus cells" )
sce <- as.SingleCellExperiment(pbmc)
sce <- slingshot(sce, clusterLabels = "cell", reducedDim = 'UMAP')
sce$cell = as.factor(sce$cell)
for(i in 1:4){
pdf(paste("slingshot",i,".pdf"),width = 5,height = 5)
plot(reducedDims(sce)$UMAP, col =brewer.pal(10, "Set3")[sce$cell] , pch = 16, cex = 0.5,asp=1) #col =brewer.pal("qual", "Set2")[sds$cell_type] 
lines(SlingshotDataSet(sce), linInd=i,lwd = 2, type = 'lineages', col = 'red')
dev.off()
}
dataset <- wrap_expression(
  expression =t(pbmc@assays$RNA@data),
  counts =t(pbmc@assays$RNA@counts)
)
time = as.numeric(gsub("PCW","",pbmc$batch))
names(time) = dataset$cell_ids 
group=cbind(dataset$cell_ids,pbmc$cell)
group=as.data.frame(group)
colnames(group) = c("cell_id","group_id")

