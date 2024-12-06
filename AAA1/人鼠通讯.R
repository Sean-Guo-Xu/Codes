library(Seurat)
library(CellChat)
load("F:\\AAA_scRNA\\human_aaa_myeloid.Rdata")
pbmc$sample
aaa = pbmc
load("F:\\AAA_scRNA\\human_normal_myeloid.Rdata")
con = pbmc
aaa$cell[aaa$cell %in% c( "Foamy/TREM2hi-Mø","Foamy/SPP1hi-Mø" )] ="Foamy-Mø"
con = subset(con, cell != "Proliferation")
con$cell[con$cell %in% c("Foamy/SPP1hi-Mø")] = "Foamy-Mø"
aaa$cell[aaa$cell %in% c( "TRM/FLOR2hi-Mø" )] ="TRM-Mø"
con$cell[con$cell %in% c("TRM/FLOR2hi-Mø")] = "TRM-Mø"
A_cellchat <-  createCellChat(aaa@assays$RNA@data, meta = aaa@meta.data, group.by = "cell")
N_cellchat <-  createCellChat(con@assays$RNA@data, meta = con@meta.data, group.by = "cell")

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

save(A_cellchat,N_cellchat,cellchat,chatlist,file="E:\\AAA_scRNA\\human_myeloid_cellchat.Rdata")
load("E:\\AAA_scRNA\\human_myeloid_cellchat.Rdata")
setwd("E:\\AAA_scRNA\\cellchat\\human")
col = rep(c("#f57245","#e6f49a","#A6CEE3","#fbe18c","#90d1a4","#33A02C","#FB9A99","#feac60"),1)
names(col) =  c("Mono","Foamy-Mø","TRM-Mø","Inf-Mø","cDC1","cDC2","IFNIC-Mø","ProInf-Mø")            
col = col[order(names(col))]
pdf("all_circle.pdf",width = 6,height =8)
netVisual_diffInteraction(cellchat, weight.scale = T,color.edge =c("#FB9A99","#90d1a4"),color.use = col,shape = "crectangle")
netVisual_diffInteraction(cellchat, weight.scale = T,color.edge =c("#FB9A99","#90d1a4"),measure = "weight",color.use = col,shape = "crectangle")
dev.off()

pdf("all_heatmap.pdf",width = 8,height =4)
gg1 <- netVisual_heatmap(cellchat,color.heatmap=c("#90d1a4","#FB9A99"),color.use = col)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.heatmap=c("#90d1a4","#FB9A99"),color.use = col)
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

library(ggsci)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("#90d1a4","#FB9A99"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#90d1a4","#FB9A99"))
gg1 + gg2
ggsave("allsignal_bar.tiff",width = 8,height = 5)
library(ComplexHeatmap)

num.link <- sapply(chatlist, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(chatlist)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(chatlist[[i]], title = names(chatlist)[i], weight.MinMax = weight.MinMax)+scale_x_continuous(limits = c(0, 0.3))+scale_y_continuous(limits = c(0, 0.32))+scale_color_manual(values = col)
}
patchwork::wrap_plots(plots = gg)
ggsave("all_cell torch.tiff",width = 10,height = 5)

for(i in unique(cellchat@meta$cell)){
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = i,color.use=c("#fbe18c","#90d1a4","#FB9A99"))
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
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = unique(con$cell), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = unique(con$cell), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(chatlist)[2]))
gg1 + gg2
ggsave("point_signal.tiff",width = 10,height = 6)

pdf("up_down_gene.pdf",width =6,height = 6)
netVisual_chord_gene(chatlist[[2]], sources.use = 4, targets.use = unique(con$cell), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))
netVisual_chord_gene(chatlist[[1]], sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(chatlist)[2]))
computeEnrichmentScore(net.up, species = 'human')
computeEnrichmentScore(net.down, species = 'human')
dev.off()


cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE)
netVisual_embeddingPairwise(cellchat, type =  "structural", label.size = 3.5,pathway.remove.show = F)+scale_color_d3()+scale_fill_d3()



netVisual_bubble(cellchat, sources.use = 4, targets.use = c("Macrophage","Monocyte","Endothelial cell"),  comparison = c(1, 2), angle.x = 45)
