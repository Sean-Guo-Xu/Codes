library(Seurat)
library(CellChat)

load("E:\\AAA_scRNA\\mouse_aaa_all.Rdata")
aaa = pbmc
aaa$index = rownames(aaa@meta.data)
load("E:\\AAA_scRNA\\mouse_aaa_myeloid.Rdata")
pbmc$index = rownames(pbmc@meta.data)
aaa$cell[aaa$index %in% pbmc$index] = pbmc$cell

load("E:\\AAA_scRNA\\mouse_normal_all.Rdata")
con = pbmc
con$index = rownames(con@meta.data)

load("E:\\AAA_scRNA\\mouse_normal_myeloid.Rdata")
pbmc$index = rownames(pbmc@meta.data)
con$cell[con$index %in% pbmc$index] = pbmc$cell


con$cell[con$cell %in% c("Foamy/Spp1hi-Mø" , "Foamy/Trem2hi-Mø")] = "Foamy-Mø"
con$cell[con$cell %in% c("TRM")] = "TRM-Mø"
aaa$cell[aaa$cell %in% c( "TRM"   )] ="TRM-Mø"
aaa$cell[aaa$cell %in% c(  "Foamy-Mø/Trem2hi", "Foamy-Mø/Spp1hi" )] ="Foamy-Mø"
aaa$cell[aaa$cell %in% c(  "Pro-Inf" )] ="Inf-Mø" 
aaa$cell[aaa$cell %in% c( "Smooth muscle cell")] = "SMC"
aaa$cell[aaa$cell %in% c( "Endothelial cell")] = "Endothelial"
aaa$cell[aaa$cell %in% c( "T cell","NK cell")] = "NK&T cell"
con$cell[con$cell %in% c( "NK/T cell")] = "NK&T cell"
aaa  = subset(aaa,cell != "Mast cell")
con = subset(con,cell != "Myeloid")
aaa  = subset(aaa,cell !=  "Granulocyte" )
aaa  = subset(aaa,cell != "Proliferation")

unique(con$cell)
unique(aaa$cell)
which(CellChatDB.mouse[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.mouse[["interaction"]] <- CellChatDB.mouse[["interaction"]][-1887,]
which(CellChatDB.mouse[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.mouse[["interaction"]] <- CellChatDB.mouse[["interaction"]][-1900,]

A_cellchat <-  createCellChat(aaa@assays$RNA@data, meta = aaa@meta.data, group.by = "cell")
N_cellchat <-  createCellChat(con@assays$RNA@data, meta = con@meta.data, group.by = "cell")

cellchat <- N_cellchat 
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
N_cellchat <- cellchat

cellchat <- A_cellchat 
cellchat@DB <- CellChatDB.mouse
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

save(A_cellchat,N_cellchat,cellchat,chatlist,file="E:\\AAA_scRNA\\mouse_all_cellchat.Rdata")
load("E:\\AAA_scRNA\\mouse_all_cellchat.Rdata")
setwd("E:\\AAA_scRNA\\cellchat\\out_mouse")
col = rep(c("#f57245","#e6f49a","#A6CEE3","#fbe18c","#90d1a4","#33A02C","#FB9A99","#feac60"),1)
names(col) =  c("Mono","Foamy-Mø","TRM-Mø","Inf-Mø","cDC1","cDC2","IFNIC-Mø","ProInf-Mø")            
col = col[order(names(col))]
pdf("all_circle.pdf",width = 6,height =8)
netVisual_diffInteraction(cellchat, weight.scale = T,color.edge =c("#FB9A99","#90d1a4"),shape = "crectangle")
netVisual_diffInteraction(cellchat, weight.scale = T,color.edge =c("#FB9A99","#90d1a4"),measure = "weight",shape = "crectangle")
dev.off()

pdf("all_heatmap.pdf",width = 8,height =4)
gg1 <- netVisual_heatmap(cellchat,color.heatmap=c("#90d1a4","#FB9A99"))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.heatmap=c("#90d1a4","#FB9A99"))
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

library(ggsci)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("#90d1a4","#FB9A99"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#90d1a4","#FB9A99"))
gg1 + gg2
ggsave("allsignal_bar.tiff",width = 10,height = 10)
library(ComplexHeatmap)

num.link <- sapply(chatlist, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(chatlist)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(chatlist[[i]], title = names(chatlist)[i], weight.MinMax = weight.MinMax)+scale_x_continuous(limits = c(0, 0.1))+scale_y_continuous(limits = c(0, 0.1))
}
patchwork::wrap_plots(plots = gg)
ggsave("all_cell torch.tiff",width = 10,height = 5)

for(i in unique(cellchat@meta$cell)){
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = i,color.use=c("#fbe18c","#90d1a4","#FB9A99"))
  ggsave(paste(i,"torch.tiff"),height = 6,width = 8)
}
pathway.union <- union(chatlist[[1]]@netP$pathways, chatlist[[2]]@netP$pathways)
source("E:\\AAA_scRNA\\net_heatmap.R")
pdf("oun_in_heatmap.pdf",width = 8,height =20)
ht1 = net_heatmap(chatlist[[1]], pattern = "all", signaling = pathway.union, title = names(chatlist)[1], width = 3.5, height = 24)
ht2 = net_heatmap(chatlist[[2]], pattern = "all", signaling = pathway.union, title = names(chatlist)[2], width = 3.5, height = 24, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

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

library(gcookbook)
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "A",ligand.logFC = 0.2, receptor.logFC = 0.2)
net.down <- subsetCommunication(cellchat, net = net, datasets = "N",ligand.logFC = -0.2, receptor.logFC = -0.2)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,sources.use = c("Foamy-Mø","Inf-Mø"),color.text =c("#90d1a4","#FB9A99") , pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))+scale_color_gradientn(colours = c("#90d1a4","white","#fbe18c","#FB9A99"))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,targets.use = c("Foamy-Mø","Inf-Mø"),color.text =c("#90d1a4","#FB9A99") , pairLR.use = pairLR.use.up,, comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(chatlist)[2]))+scale_color_gradientn(colours = c("#90d1a4","white","#fbe18c","#FB9A99"))
gg1 + gg2
ggsave("point_signal.tiff",width = 20,height = 11)

