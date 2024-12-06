
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
setwd("E:\\brain")
load("cancer.Rdata")
unique(pbmc$subgroup)
groupname = "SHH"

pbmc = subset(pbmc,sample == groupname)
pbmc = subset(pbmc, cell %in% c("Endothelial","Pericytes","Neuron","Neuronal progenitor cell","Cycling cells"),invert=T)
color24 = c("#E69F00", "#56B4E9", "#34495E", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#8E44AD", "#1ABC9C", "#F39C12", "#2ECC71", "#E74C3C", "#3498DB", "#9B59B6", "#16A085", "#F1C40F", "#56B4E9", "#E67E22", "#C0392B","#2980B9", "#27AE60","#8E44AD", "#2C3E50")
imcellname=c("Naive CD4+T cell","Memory CD4+T cell","CD8+T cell","Treg","B cell","NK cell","Microglial","Myeloid DC","cDC","Inf DC","pDC","M2","M1", "OPC" , "APC" ,"G3_MB","G4_MB",  "SHH","WNT","Endothelial","Pericytes","Neuron", "Neuronal progenitor cell", "Cycling cells" )
names(color24) = imcellname
color = color24[names(color24) %in% unique(pbmc$cell)]
color = color[order(names(color))]
pbmc = NormalizeData(pbmc)

labels <- Idents(pbmc)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = pbmc, group.by = "cell")
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
gc()
cellchat <- identifyOverExpressedInteractions(cellchat)
gc()
options(future.globals.maxSize = 8000 * 1024^2)
cellchat <- computeCommunProb(cellchat, type = "triMean")
gc()
cellchat <- filterCommunication(cellchat, min.cells = 10)
gc()
cellchat <- computeCommunProbPathway(cellchat)
gc()
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

save(cellchat,file=paste(groupname,"cellchat.Rdata"))
load(paste(groupname,"cellchat.Rdata"))


pdf(paste(groupname,"cellchat_circle.pdf"),width = 8,height = 10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = color,alpha.edge=1)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = color,alpha.edge=1)
dev.off()

pdf(paste(groupname,"cellchat_single_circle.pdf"),width = 16,height = 20)
mat <- cellchat@net$weight
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i],color.use = color,alpha.edge=1)
}
dev.off()
pdf(paste(groupname,"cellchat_heatmap.pdf"),width = 7,height = 6.5)
netVisual_heatmap(cellchat, color.heatmap = c("#5371b3","#E31A1C"),color.use = color)
netVisual_heatmap(cellchat, color.heatmap =  c("#5371b3","#E31A1C"),measure = "weight",color.use = color)
dev.off()

cancercell= c(groupname)
immunecell = c("M1","M2", "Memory CD4+T cell" ,"Naive CD4+T cell" ,"CD8+T cell","Treg" , "NK cell","B cell", "Microglial" )
immunecell = immunecell[immunecell %in% unique(pbmc$cell)]
neuron = c("Neuron", "Astrocyte" ,"Oligodendrocyte",groupname)
neuron = neuron[neuron %in% unique(pbmc$cell)]
netVisual_bubble(cellchat, sources.use = cancercell, targets.use = immunecell, remove.isolate = FALSE,color.heatmap = "viridis")
ggsave(paste(groupname,"can-imm_dotplot.png"),width = 6,height = 8)
netVisual_bubble(cellchat, sources.use = immunecell, targets.use = cancercell, remove.isolate = FALSE,color.heatmap = "viridis")
ggsave(paste(groupname,"imm-can_dotplot.png"),width = 5,height = 6)
netVisual_bubble(cellchat, sources.use = immunecell, targets.use = immunecell, remove.isolate = FALSE,color.heatmap = "viridis")
ggsave(paste(groupname,"imm-imm_dotplot.png"),width = 12,height = 12)
netVisual_bubble(cellchat, sources.use = neuron, targets.use = neuron, remove.isolate = FALSE)
ggsave(paste(groupname,"neuron-can_dotplot.png"),width = 6,height = 9)

###############对比分析######
color24 = c("#E69F00", "#56B4E9", "#34495E", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#8E44AD", "#1ABC9C", "#F39C12", "#2ECC71", "#E74C3C", "#3498DB", "#9B59B6", "#16A085", "#F1C40F", "#56B4E9", "#E67E22", "#C0392B","#2980B9", "#27AE60","#8E44AD", "#2C3E50")
imcellname=c("Naive CD4+T cell","Memory CD4+T cell","CD8+T cell","Treg","B cell","NK cell","Microglial","Myeloid DC","cDC","Inf DC","pDC","M2","M1", "OPC" , "APC" ,"G3_MB","G4_MB",  "SHH","WNT","Endothelial","Pericytes","Neuron", "Neuronal progenitor cell", "Cycling cells" )
names(color24) = imcellname


load("cancer.Rdata")
pbmc = subset(pbmc,orig.ident %in% c("X966","X966.2"))
pbmc$cancer = pbmc$orig.ident
pbmc$cancer[pbmc$cancer == "X966"] = "Original"
pbmc$cancer[pbmc$cancer == "X966.2"] = "Recurrent"

check = table(pbmc$cancer,pbmc$cell)
pdf("compare_cellnumber.pdf",width = 5.5,height = 2.2)
pheatmap::pheatmap(log2(check+1),cluster_rows = F,cluster_cols = F,display_numbers = check,legend = F)
dev.off()
check = as.data.frame(check)
move = check[check$Freq == 0,] 
pbmc = subset(pbmc, cell %in% move$Var2,invert=T)

ori= subset(pbmc ,cancer == "Original")
ori =NormalizeData(ori)
cellchat <- createCellChat(object = ori, group.by = "cell")
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

rec = subset(pbmc ,cancer == "Recurrent")
rec =NormalizeData(rec)
cellchat <- createCellChat(object = rec, group.by = "cell")
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

color = color24[names(color24) %in% unique(cellchat@meta$cell)]
color = color[order(names(color))]
cellchatlist <- list(Original = cellchat1, Recurrent = cellchat2)
cellchat <- mergeCellChat(cellchatlist, add.names = names(cellchatlist))
g1<-compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = c("#5371b3","#E31A1C"))
g2<-compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use = c("#5371b3","#E31A1C"))
g1+g2
ggsave("sum_barplot.png",width = 4,height = 4)
cellname = c("GP4","CD8+T cell","M2")
color3 = color[names(color) %in% cellname]
pdf("compare_circle.pdf",width = 5,height = 5)
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color3,label.edge = T,remove.isolate = T)
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = color3,measure = "weight",label.edge = T,remove.isolate = T)
netVisual_heatmap(cellchat, color.heatmap = c("#5371b3","#E31A1C"),color.use = color)
netVisual_heatmap(cellchat, color.heatmap =  c("#5371b3","#E31A1C"),measure = "weight",color.use = color)
dev.off()

rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,color.use = c("#5371b3","#E31A1C"))
ggsave("compare_barplot.png",width = 4,height = 6)
p1<-netVisual_bubble(cellchat, sources.use = cellname, targets.use = cellname,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Recurrent", angle.x = 45, remove.isolate = T,color.heatmap = "viridis",color.text =c("#5371b3","#E31A1C") )
p2<-netVisual_bubble(cellchat, sources.use = cellname, targets.use = cellname,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Recurrent", angle.x = 45, remove.isolate = T,color.heatmap = "viridis",color.text =c("#5371b3","#E31A1C") )
p1+p2
ggsave("compare_In_De_pathway.png",width = 10,height = 8)
pos.dataset = "Recurrent"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Recurrent",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Original",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
p1<-netVisual_bubble(cellchat, angle.x = 45,pairLR.use = pairLR.use.up, sources.use = cellname, targets.use = cellname, comparison = c(1, 2),   remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(cellchatlist)[2]),color.heatmap = "viridis",color.text =c("#5371b3","#E31A1C"))
pairLR.use.down = net.down[, "interaction_name", drop = F]
p2<-netVisual_bubble(cellchat, angle.x = 45,pairLR.use = pairLR.use.down, sources.use = cellname, targets.use = cellname, comparison = c(1, 2),   remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(cellchatlist)[2]),color.heatmap = "viridis",color.text =c("#5371b3","#E31A1C"))
p1+p2
ggsave("compare_Up_Down_pathway.png",width = 10,height = 6)
cellchat@meta$datasets = factor(cellchat@meta$cancer, levels = c("Original", "Recurrent")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin",color.use =  c("#5371b3","#E31A1C"))
ggsave("compare_CXCL_vlnplot.png",width=5,height = 5)
plotGeneExpression(cellchat, signaling = "MHC-I", split.by = "datasets", colors.ggplot = T, type = "violin",color.use =  c("#5371b3","#E31A1C"))
ggsave("compare_MHC-I_vlnplot.png",width=5,height = 8)
save(cellchat,cellchatlist,file="Recurrent_compare_cellchat.Rdata")
load("Recurrent_compare_cellchat.Rdata")
