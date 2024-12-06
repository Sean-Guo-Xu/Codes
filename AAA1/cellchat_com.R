
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
setwd("F:\\AAA_scRNA")
load("7AAA.Rdata")
aaa = pbmc
load("3Normal.Rdata")
cellname = intersect(pbmc$cell,aaa$cell)
pbmc = subset(pbmc, cell %in% cellname)
aaa = subset(aaa, cell %in% cellname)



check = table(pbmc$cancer,pbmc$cell)

aaa =NormalizeData(aaa)
cellchat <- createCellChat(object = aaa, group.by = "cell")
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat1 = cellchat
gc()
pbmc =NormalizeData(pbmc)
cellchat <- createCellChat(object = pbmc, group.by = "cell")
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat2 = cellchat
color2=c( "#6BB7CA",  "#E07B54")
library(ggsci)
color = pal_d3()(9)
names(color) = unique(cellchat@meta$cell)
color = color[order(names(color))]
cellchatlist <- list(Normal = cellchat2, AAA = cellchat1)
cellchat <- mergeCellChat(cellchatlist, add.names = names(cellchatlist))
g1<-compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = color2)
g2<-compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use = color2)
g1+g2
ggsave("sum_barplot.png",width = 4,height = 4)

pdf("compare_all_circle.pdf",width = 10,height = 10)
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color,label.edge = T,margin=c(1,1,1,1))
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color,measure = "weight",label.edge = T,margin=c(1,1,1,1))
dev.off()
pdf("compare_all_heatmap.pdf",width =6.5,height = 6)

netVisual_heatmap(cellchat, color.heatmap = color2,color.use = color)
netVisual_heatmap(cellchat, color.heatmap =  color2,measure = "weight",color.use = color)
dev.off()

rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,color.use = color2)
ggsave("compare_barplot.png",width = 4,height = 6)
netVisual_bubble(cellchat,  signaling = "MHC-II", comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in AAA", angle.x = 45, remove.isolate = T,color.heatmap = "viridis",color.text =color2 )
ggsave("compare_In_pathway.png",width = 14,height = 5)
netVisual_bubble(cellchat, signaling = "MHC-II",comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in AAA", angle.x = 45, remove.isolate = T,color.heatmap = "viridis",color.text =color2 )
ggsave("compare_De_pathway.png",width = 5,height = 5)

pos.dataset = "AAA"
db = cellchat@DB$interaction
db  = db[db$pathway_name %in% "MHC-II",]
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "AAA",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[net.up$interaction_name %in% db$interaction_name, "interaction_name", drop = F]
netVisual_bubble(cellchat, angle.x = 45,pairLR.use = pairLR.use.up,  comparison = c(1, 2),   remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(cellchatlist)[2]),color.heatmap = "viridis",color.text =color2)
ggsave("compare_Up_pathway.png",width = 10,height = 5)
pairLR.use.down = net.down[net.down$interaction_name %in% db$interaction_name, "interaction_name", drop = F]
netVisual_bubble(cellchat, angle.x = 45,pairLR.use = pairLR.use.down,  comparison = c(1, 2),   remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(cellchatlist)[2]),color.heatmap = "viridis",color.text =color2)
ggsave("compare_Down_pathway.png",width = 10,height = 4)
weight.max <- getMaxWeight(cellchatlist, slot.name = c("netP"), attribute = "MHC-II") 
pdf("com_path_circle.pdf",width = 10,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchatlist)) {
  netVisual_aggregate(cellchatlist[[i]],color.use = color ,signaling = "MHC-II", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste("MHC-II", names(cellchatlist)[i]),scale = T)
}
dev.off()
pdf("com_path_hearmap.pdf",width = 12,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cellchatlist)) {
  ht[[i]] <- netVisual_heatmap(cellchatlist[[i]], signaling = "MHC-II",title.name = paste("MHC-II", "signaling ",names(cellchatlist)[i]),color.use = color,color.heatmap = color2)
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

plotGeneExpression(cellchat, signaling = "MHC-II", split.by = "datasets", colors.ggplot = T, type = "violin",color.use =  color2)
ggsave("compare_path_vlnplot.png",width=5,height = 6)
save(cellchat,cellchatlist,file="Aorta_cellchat.Rdata")
load('Aorta_cellchat.Rdata')
pdf("Aorta_circle_com.pdf",width = 6,height = 7)
for (i in 1:2) {
  netVisual_aggregate(cellchatlist[[i]], signaling = "MHC-II", layout = "chord", signaling.name = paste("MHC-II", names(cellchatlist)[i]),color.use = color,remove.isolate = T,cell.order = unique(cellchatlist[[i]]@meta$cell))
}
dev.off()
###############blood#########

library(CellChat)
library(patchwork)
library(Seurat)
library(ggsci)
options(stringsAsFactors = FALSE)
setwd("F:\\AAA_scRNA")
load("blood/A_blood.Rdata")
aaa = pbmc
load("blood/N_blood.Rdata")

cellname = intersect(pbmc$cell,aaa$cell)
pbmc = subset(pbmc, cell %in% cellname)
aaa = subset(aaa, cell %in% cellname)


aaa =NormalizeData(aaa)
cellchat <- createCellChat(object = aaa, group.by = "cell")
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat1 = cellchat
gc()
pbmc =NormalizeData(pbmc)
cellchat <- createCellChat(object = pbmc, group.by = "cell")
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat2 = cellchat
color2=c( "#6BB7CA",  "#E07B54")

library(ggsci)
color = pal_d3()(9)
names(color) = unique(cellchat@meta$cell)
color = color[order(names(color))]
cellchatlist <- list(Normal = cellchat2, AAA = cellchat1)
cellchat <- mergeCellChat(cellchatlist, add.names = names(cellchatlist))
g1<-compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = color2)
g2<-compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use = color2)
g1+g2
ggsave("sum_barplot.png",width = 4,height = 4)

pdf("compare_all_circle.pdf",width = 10,height = 10)
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color,label.edge = T,margin=c(1,1,1,1))
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color,measure = "weight",label.edge = T,margin=c(1,1,1,1))
dev.off()
pdf("compare_all_heatmap.pdf",width =6.5,height = 6)

netVisual_heatmap(cellchat, color.heatmap = color2,color.use = color)
netVisual_heatmap(cellchat, color.heatmap =  color2,measure = "weight",color.use = color)
dev.off()

rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,color.use = color2)
ggsave("compare_barplot.png",width = 4,height = 6)
netVisual_bubble(cellchat,  signaling = "MHC-II", comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in AAA", angle.x = 45, remove.isolate = T,color.heatmap = "viridis",color.text =color2 )
ggsave("compare_In_pathway.png",width = 13,height = 5)
netVisual_bubble(cellchat,targets.use = c("Dendritic cell",  "Monocyte" ), signaling = "MHC-II",,comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in AAA", angle.x = 45, remove.isolate = T,color.heatmap = "viridis",color.text =color2 )
ggsave("compare_De_pathway.png",width = 7,height = 5)

pos.dataset = "AAA"
db = cellchat@DB$interaction
db  = db[db$pathway_name %in% "MHC-II",]
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "AAA",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[net.up$interaction_name %in% db$interaction_name, "interaction_name", drop = F]
netVisual_bubble(cellchat, angle.x = 45,pairLR.use = pairLR.use.up,  comparison = c(1, 2),   remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(cellchatlist)[2]),color.heatmap = "viridis",color.text =color2)
ggsave("compare_Up_pathway.png",width = 10,height = 5)
pairLR.use.down = net.down[net.down$interaction_name %in% db$interaction_name, "interaction_name", drop = F]
netVisual_bubble(cellchat, angle.x = 45,pairLR.use = pairLR.use.down,  comparison = c(1, 2),   remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(cellchatlist)[2]),color.heatmap = "viridis",color.text =color2)
ggsave("compare_Down_pathway.png",width = 10,height = 4)
weight.max <- getMaxWeight(cellchatlist, slot.name = c("netP"), attribute = "MHC-II") 
pdf("com_path_circle.pdf",width = 10,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchatlist)) {
  netVisual_aggregate(cellchatlist[[i]],color.use = color ,signaling = "MHC-II", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste("MHC-II", names(cellchatlist)[i]),scale = T)
}
dev.off()
pdf("com_path_hearmap.pdf",width = 12,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cellchatlist)) {
  ht[[i]] <- netVisual_heatmap(cellchatlist[[i]], signaling = "MHC-II",title.name = paste("MHC-II", "signaling ",names(cellchatlist)[i]),color.use = color,color.heatmap = color2)
}
ComplexHeatmap::draw(ht[[2]] + ht[[1]], ht_gap = unit(0.5, "cm"))
dev.off()

plotGeneExpression(cellchat, signaling = "MHC-II", split.by = "datasets", colors.ggplot = T, type = "violin",color.use =  color2)
ggsave("compare_path_vlnplot.png",width=5,height = 6)
save(cellchat,cellchatlist,file="PBMC_cellchat.Rdata")
load("PBMC_cellchat.Rdata")

pdf("circle_com.pdf",width = 6,height = 7)
for (i in 1:2) {
netVisual_aggregate(cellchatlist[[i]], signaling = "MHC-II", layout = "chord", signaling.name = paste("MHC-II", names(cellchatlist)[i]),color.use = color,remove.isolate = T,cell.order = unique(cellchatlist[[i]]@meta$cell))
  }
dev.off()

