library(Seurat)
library(CellChat)
load("E:\\AAA_scRNA\\blood\\A_blood.Rdata")
blood = subset(pbmc,cell==c("Monocyte"))
blood$cell = "Monocyte"
load("E:\\AAA_scRNA\\8AAA.Rdata")
artery = subset(pbmc , cell != "Unsure")
pbmc = merge(blood,artery)
bcounts <- GetAssayData(blood, assay = "RNA")
acounts<- GetAssayData(artery, assay = "RNA")
name = rownames(bcounts)[rownames(bcounts) %in% rownames(acounts)]
pbmc = subset(pbmc,features = name)
pbmc = NormalizeData(pbmc)
data  <- pbmc@assays$RNA@data
A_cellchat <-  createCellChat(pbmc@assays$RNA@data, meta = pbmc@meta.data, group.by = "cell")

load("E:\\AAA_scRNA\\blood\\N_blood.Rdata")
blood = subset(pbmc,cell==c("Monocyte"))
blood$cell = "Monocyte"
load("E:\\AAA_scRNA\\6Normal.Rdata")
artery = subset(pbmc , cell != "Unsure")
pbmc = merge(blood,artery)
bcounts <- GetAssayData(blood, assay = "RNA")
acounts<- GetAssayData(artery, assay = "RNA")
name = rownames(bcounts)[rownames(bcounts) %in% rownames(acounts)]
pbmc = subset(pbmc,features = name)
pbmc = NormalizeData(pbmc)
data  <- pbmc@assays$RNA@data
N_cellchat <-  createCellChat(pbmc@assays$RNA@data, meta = pbmc@meta.data, group.by = "cell")

cellchat <- N_cellchat 
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
N_cellchat <- cellchat

cellchat <- A_cellchat 
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
A_cellchat <- cellchat
chatlist <- list(N = N_cellchat,A = A_cellchat)
cellchat <- mergeCellChat(chatlist, add.names = names(chatlist), cell.prefix = TRUE)

save(A_cellchat,N_cellchat,cellchat,chatlist,file="E:\\AAA_scRNA\\cellchat.Rdata")
load("E:\\AAA_scRNA\\cellchat.Rdata")
setwd("E:\\AAA_scRNA\\cellchat")

pdf("all_circle.pdf",width = 6,height =8)
netVisual_diffInteraction(cellchat, weight.scale = T,color.edge =c("#D62728FF","#1F77B4FF"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.edge =c("#D62728FF","#1F77B4FF"))
dev.off()

pdf("all_heatmap.pdf",width = 12,height =6)
gg1 <- netVisual_heatmap(cellchat,color.heatmap=c("#1F77B4FF","#D62728FF"))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.heatmap=c("#1F77B4FF","#D62728FF"))
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

library(ggsci)
pal_d3()(2)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("#1F77B4FF","#D62728FF"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#1F77B4FF","#D62728FF"))
gg1 + gg2
ggsave("allsignal_bar.tiff",width = 8,height = 5)
library(ComplexHeatmap)

num.link <- sapply(chatlist, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(chatlist)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(chatlist[[i]], title = names(chatlist)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
ggsave("all_cell torch.tiff",width = 10,height = 5)

for(i in unique(cellchat@meta$cell)){
 netAnalysis_signalingChanges_scatter(cellchat, idents.use = i,color.use=c("#BCBD22FF","#1F77B4FF","#D62728FF"))
ggsave(paste(i,"torch.tiff"),height = 6,width = 8)
}
library(Seurat)
library("umap-learn")
library(uwot)
library(ggsci)
cellchat <- computeNetSimilarityPairwise(cellchat, type = c("functional"))
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,pathway.remove.show = F)+scale_color_d3()+scale_fill_d3()
ggsave("all_signal_torch.tiff",width = 6,height = 5)
pos.dataset = "A"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "A",ligand.logFC = 0.2, receptor.logFC = 0.2)
net.down <- subsetCommunication(cellchat, net = net, datasets = "N",ligand.logFC = -0.2, receptor.logFC = -0.2)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(chatlist)[2]))
gg1 + gg2
ggsave("point_signal.tiff",width = 10,height = 6)

pdf("up_down_gene.pdf",width =6,height = 6)
netVisual_chord_gene(chatlist[[2]], sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))
netVisual_chord_gene(chatlist[[1]], sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(chatlist)[2]))
computeEnrichmentScore(net.up, species = 'human')
computeEnrichmentScore(net.down, species = 'human')
dev.off()


cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE)
netVisual_embeddingPairwise(cellchat, type =  "structural", label.size = 3.5,pathway.remove.show = F)+scale_color_d3()+scale_fill_d3()



netVisual_bubble(cellchat, sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"),  comparison = c(1, 2), angle.x = 45)
