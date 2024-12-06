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
aaa$cell[aaa$cell %in% c( "cDC1","cDC2" )] ="cDC"
con$cell[con$cell %in% c("cDC1","cDC2" )] = "cDC"
con$cell[con$cell %in% c("IFNIC-Mø"  )] = "IFNIC"
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
setwd("F:\\AAA_scRNA\\cellchat\\human")
col = rep(c("#f57245","#e6f49a","#A6CEE3","#fbe18c","#90d1a4","#FB9A99","#feac60"),1)
names(col) =  c("Mono","Foamy-Mø","TRM-Mø","Inf-Mø","cDC","IFNIC-Mø","ProInf-Mø")            
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
source("F:\\AAA_scRNA\\net_heatmap.R")
pathway.union <- union(chatlist[[1]]@netP$pathways, chatlist[[2]]@netP$pathways)
pdf("oun_in_heatmap.pdf",width = 8,height =20)
ht1 = net_heatmap(chatlist[[1]], pattern = "all", signaling = pathway.union, title = names(chatlist)[1], width = 3.5, height = 24,color.use = col)
ht2 = net_heatmap(chatlist[[2]], pattern = "all", signaling = pathway.union, title = names(chatlist)[2], width = 3.5, height = 24, color.heatmap = "OrRd",color.use = col)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
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

