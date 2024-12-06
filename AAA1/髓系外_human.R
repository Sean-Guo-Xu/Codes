library(Seurat)
library(CellChat)

load("F:\\AAA_scRNA\\7AAA.Rdata")
aaa = pbmc
load("F:\\AAA_scRNA\\human_aaa_myeloid.Rdata")
pbmc$index = rownames(pbmc@meta.data)
aaa$index = rownames(aaa@meta.data)
aaa$cell[aaa$index %in% pbmc$index] = pbmc$cell
aaa = subset(aaa, cell != "Macrophage")
aaa = subset(aaa, cell != "Monocyte")
aaa = subset(aaa, cell != "Proliferation")
aaa$cell[aaa$cell %in% c("NK cell","T cell")] ="NK/T cell"
load("F:\\AAA_scRNA\\6Normal.Rdata")
con = pbmc

load("F:\\AAA_scRNA\\human_normal_myeloid.Rdata")
con$cell[con$index %in% pbmc$index] = pbmc$cell

aaa$cell[aaa$cell %in% c( "Foamy/TREM2hi-Mø","Foamy/SPP1hi-Mø" )] ="Foamy-Mø"
con = subset(con, cell != "Proliferation")
con$cell[con$cell %in% c("Foamy/SPP1hi-Mø")] = "Foamy-Mø"
aaa$cell[aaa$cell %in% c( "TRM/FLOR2hi-Mø" )] ="TRM"
con$cell[con$cell %in% c("TRM/FLOR2hi-Mø")] = "TRM"
con = subset(con,cell != "Proliferation")
aaa = subset(aaa, cell != "Plasma cell")
aaa = subset(aaa, cell != "Myeloid")
aaa = subset(aaa, cell != "Plasmacytoid DC")
con$cell[con$cell %in% c("cDC1","cDC2")] = "cDC"
aaa$cell[aaa$cell %in% c("cDC1","cDC2")] = "cDC"
aaa$cell[aaa$cell %in% c("SMC")] = "Smooth muscle cell"
aaa$cell[aaa$cell %in% c("IFNIC")] = "IFNIC-Mø"
aaa$cell[aaa$cell %in% c("TRM/FLOR2hi-Mø")] = "TRM"

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

save(A_cellchat,N_cellchat,cellchat,chatlist,file="F:\\AAA_scRNA\\human_all_cellchat.Rdata")
load("F:\\AAA_scRNA\\human_all_cellchat.Rdata")
setwd("F:\\AAA_scRNA\\cellchat\\out_human")
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
ggsave("allsignal_bar.tiff",width = 10,height = 6)
library(ComplexHeatmap)

num.link <- sapply(chatlist, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()

  gg[[1]] <- netAnalysis_signalingRole_scatter(chatlist[[1]], title = names(chatlist)[1], weight.MinMax = weight.MinMax)+scale_x_continuous(limits = c(0, 0.1))+scale_y_continuous(limits = c(0, 0.1))
  gg[[2]] <- netAnalysis_signalingRole_scatter(chatlist[[2]], title = names(chatlist)[2], weight.MinMax = weight.MinMax)+scale_x_continuous(limits = c(0, 0.1))+scale_y_continuous(limits = c(0, 0.1))
  
patchwork::wrap_plots(plots = gg)
ggsave("all_cell torch.tiff",width = 10,height = 5)

for(i in unique(cellchat@meta$cell)){
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = i,color.use=c("#fbe18c","#90d1a4","#FB9A99"))
  ggsave(paste(i,"torch.tiff"),height = 6,width = 8)
}
source("F:\\AAA_scRNA\\net_heatmap.R")
pdf("oun_in_heatmap.pdf",width = 8,height =20)

pathway.union <- union(chatlist[[1]]@netP$pathways, chatlist[[2]]@netP$pathways)
ht1 = net_heatmap(chatlist[[1]], pattern = "all", signaling = pathway.union, title = names(chatlist)[1], width = 3.5, height = 24)
ht2 = net_heatmap(chatlist[[2]], pattern = "all", signaling = pathway.union, title = names(chatlist)[2], width = 3.5, height = 24, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
library(gcookbook)

pos.dataset = "A"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "A",ligand.logFC = 0.2, receptor.logFC = 0.2)
net.down <- subsetCommunication(cellchat, net = net, datasets = "N",ligand.logFC = -0.2, receptor.logFC = -0.2)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,sources.use = c("Foamy-Mø","Inf-Mø"),color.text =c("#90d1a4","#FB9A99") , pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))+scale_color_gradientn(colours = c("#90d1a4","white","#fbe18c","#FB9A99"))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,targets.use = c("Foamy-Mø","Inf-Mø"),color.text =c("#90d1a4","#FB9A99") , pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(chatlist)[2]))+scale_color_gradientn(colours = c("#90d1a4","white","#fbe18c","#FB9A99"))
gg1 + gg2
ggsave("point_signal.tiff",width = 20,height = 7)

