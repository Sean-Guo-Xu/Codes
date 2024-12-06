library("Seurat")
library("harmony")
setwd("D:\\生信\\target数据分析\\scRNA")
dir="GSE162454\\1"
list.files(dir)
counts <- Read10X(data.dir = dir)
pbmc <- CreateSeuratObject(counts = counts)

sample = rep("OS1",length(pbmc$orig.ident))


for(i in 2:6){
  dir=paste("GSE162454\\",i,sep="")
  counts = Read10X(data.dir = dir)
  pbmc=merge(pbmc,CreateSeuratObject(counts = counts))
  sample=c(sample,rep(paste("OS",i,sep = ""),length(colnames(counts))))
}
pbmc$sample=sample
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 300& nFeature_RNA <4500 & percent.mt < 10)  
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基???
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 14
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.15)                  #??????????,??????????????
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)   


DimPlot(pbmc,reduction="umap",group.by = "sample")

pbmc = pbmc %>% RunHarmony("sample", plot_convergence = TRUE)#耗时1min
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:pcSelect) %>%
  FindNeighbors(reduction = "harmony", dims = 1:pcSelect) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DimPlot(pbmc,reduction="umap",group.by = "seurat_clusters",label=T)
dev.off()
save(pbmc,file="OS.Rdata")

logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter
)
write.table(pbmc.markers,file="15Markers.txt",sep="\t",row.names=T,quote=F) 
############自动注释###############
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref, labels = ref$label.fine)
pbmc@meta.data$auto=cellpred$labels
DimPlot(pbmc,reduction="umap",group.by = "seurat_clusters",label=T)
###############手动注释###########
l=cbind(pbmc$orig.ident,pbmc$seurat_clusters)
l[which(l[,2] %in% c("1","2","6","7","19")),1]="myeloid cells"
l[which(l[,2] %in% c("4","8","9","17","20","21","12")),1]="OS cells"
l[which(l[,2] %in% c("3","5")),1]="NK/T cells"
l[which(l[,2] %in% c("13","11")),1]="OCs"
l[which(l[,2] %in% c("10")),1]="plasma cells"
l[which(l[,2] %in% c("14")),1]="endothelial cells"
l[which(l[,2] %in% c("15")),1]="B cells"
l[which(l[,2] %in% c("16","18")),1]="CAFs"
pbmc$cell = l[,1]
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DimPlot(pbmc,reduction="umap",group.by = "cell",label=T)
dev.off()
####################OS细胞###################
load("OS.Rdata")
pbmc=subset(pbmc,cell == "OS cells")
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 15
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.06)                  #??????????,??????????????
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)   
DimPlot(pbmc,reduction="umap",group.by = "seurat_clusters",label=T)

logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter
)
write.table(pbmc.markers,file="OS Markers.txt",sep="\t",row.names=T,quote=F) 

l=cbind(pbmc$cell,pbmc$seurat_clusters)
l[which(l[,2] %in% "4"),1] = "High Risk OS"
l[which(l[,2] %in% "5"),1] = "Low Risk OS"
pbmc$os = l[,1]
save(pbmc,file="onlyos.Rdata")
###################高低Os轨迹######################
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(monocle)
load("onlyos.Rdata")
pbmc=subset(pbmc,os %in% c("High OS","Low OS" ))
logFCfilter=0.5
pbmc.markers <- FindAllMarkers(pbmc,ide ,group.by="os",
                               test.use = "DESeq2" ,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter
)
pbmc.markers = FindMarkers(pbmc, ident.1 = "High OS",ident.2 = "Low OS",group.by = "os")
pbmc.markers=FindMarkers(pbmc,ident.1="Normal",ident.2="AML",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0.1)

pbmc.markers=FindMarkers(i,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0)

df_id<-bitr(rownames(pbmc.markers), #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
marker=cbind(rownames(pbmc.markers),pbmc.markers)
marker=marker[which(marker$`rownames(pbmc.markers)` %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker$avg_log2FC
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
GO = gseGO(gene_fc,ont= "ALL",OrgDb = org.Hs.eg.db)

kk_gse <- GO
kk_gse_cut = kk_gse[grep("collagen",kk_gse$Description)]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 1,]
pdf("细胞差异gsea.pdf",width=12,height=15.5)
gseaplot2(kk_gse,
          kk_gse_cut_up$ID,#瀵岄泦鐨処D缂栧???
          color = "red",#GSEA绾挎潯棰滆壊
          base_size = 20,#鍩虹瀛椾綋澶у???
          rel_heights = c(1.5, 0.5, 1),#鍓浘鐨勭浉瀵归珮搴?
          subplots = 1:3, #瑕佹樉绀哄摢浜涘壇鍥? 濡俿ubplots=c(1,3) #鍙绗竴鍜岀涓変釜????
          ES_geom = "line",#enrichment score鐢ㄧ嚎杩樻槸鐢ㄧ???"dot"
) #鏄剧ずpvalue绛変???
dev.off()
write.table(pbmc.markers,file="HL Markers.txt",sep="\t",row.names=T,quote=F) 


monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)


monocle.sample=read.table("07.monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("07.monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("OS Markers.txt",sep="\t",header=T,check.names=F)
marker=marker[marker$cluster %in% c("4","5"),]
library("dplyr")

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表???
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)

cds <- estimateSizeFactors(cds)
memory.limit(size=300000000000)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')

cds <- orderCells(cds)
load(file="hlcds.Rdata")

pdf(file="pseudotime.pdf",width=5,height=5)
plot_cell_trajectory(cds,color_by = "Pseudotime",show_backbone = T)
dev.off()
pdf(file="check.trajectory.pdf",width=5,height=5)
plot_complex_cell_trajectory(cds, x = 1, y = 2,
                             color_by = "os",show_branch_points =T)+
  theme(legend.title = element_blank()) 

plot_cell_trajectory(cds,color_by = "os")
dev.off()
library(ggplot2)
library(scales)
pdf(file="exp trajectory.pdf",width=5,height=5)
pData(cds)$COL5A2 = log2( exprs(cds)['COL5A2',]+1)
plot_cell_trajectory(cds, color_by = "COL5A2")  +scale_colour_gradient2(low = muted("cyan3"),
                                                                        
                                                                        mid = "cyan3",
                                                                        
                                                                        high = muted("brown3"))
pData(cds)$IGF1R = log2( exprs(cds)['IGF1R',]+1)
plot_cell_trajectory(cds, color_by = "IGF1R")  +scale_colour_gradient2(
  
  mid = "cyan3",
  
  high = muted("brown3"))
dev.off()

test=cds
expressed_genes=see[see$use_for_ordering,2] 
pseudotime_de <- differentialGeneTest(test[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(test[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]


saveRDS(test, file = "test_monocle.rds")
write.table(pseudotime_de, file = "pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(states_de, file = "states_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

test_genes=c("COL5A2","IGF1R")
pdf(file="gene time.pdf",width=6,height=4)
plot_genes_in_pseudotime(cds[test_genes,],color_by = "os")
dev.off()


time = differentialGeneTest(cds,cores = 6,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(time, "Time_diff_all.csv", row.names = F)

p=plot_pseudotime_heatmap(cds[time$gene_short_name,], num_clusters=2, show_rownames=T, return_heatmap=T)
see = p$tree_row
clusters <- cutree(p$tree_row, k = 2)
clustering <- data.frame(clusters)
clustering = cbind(clustering,p$tree_row$order)
clustering = clustering[order(rownames(clustering)),]
time = time[which(rownames(time) %in% rownames(clustering)),]
time = time[order(rownames(time)),]
time = cbind(time,clustering)

time1 = time[time$clusters %in% "1",]
time1 = time1[order(time1$qval),]
time1$qval = 1-time1$qval
time2 = time[time$clusters %in% "2",]
time2 = time2[order(time2$qval,decreasing = T),]
time2$qval = time2$qval - 1
time = rbind(time1,time2)

ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

#GSEA分析
library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
df_id<-bitr(rownames(time), #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
marker=cbind(rownames(time),time)
marker=marker[which(marker$`rownames(time)` %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker$qval
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
KEGG <- gseKEGG(gene_fc, organism = "hsa")
GO = gseGO(gene_fc,ont= "ALL",OrgDb = org.Hs.eg.db)
KEGG
kk_gse <- GO
kk_gse_entrez <- KEGG_kk_entrez

###鏉′欢绛涢???? 
#涓€鑸涓簗NES|>1锛孨OM pvalue<0.05锛孎DR锛坧adj????<0.25鐨勯€氳矾鏄樉钁楀瘜闆嗙殑
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.05 & abs(kk_gse$NES)>1 ]
kk_gse_cut = kk_gse[grep("collagen",kk_gse$Description)]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

#閫夋嫨灞曠幇NES鍓嶅嚑涓€氳矾 
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
pdf(file="时间相关 gsea.pdf",width=12,height=15.5)
gseaplot2(kk_gse,
          up_gsea$ID,#瀵岄泦鐨処D缂栧???
          color = "red",#GSEA绾挎潯棰滆壊
          base_size = 20,#鍩虹瀛椾綋澶у???
          rel_heights = c(1.5, 0.5, 1),#鍓浘鐨勭浉瀵归珮搴?
          subplots = 1:3, #瑕佹樉绀哄摢浜涘壇鍥? 濡俿ubplots=c(1,3) #鍙绗竴鍜岀涓変釜????
          ES_geom = "line",#enrichment score鐢ㄧ嚎杩樻槸鐢ㄧ???"dot"
) #鏄剧ずpvalue绛変???

dev.off()

diff = read.table("D:\\生信\\target数据分析\\新mi-m\\all.xls",header = T,sep = "\t")
rownames(diff) = diff$gene
df_id<-bitr(rownames(diff), #??????????df??????????SYMBOL??
            fromType = "SYMBOL",#????????ID????
            toType = "ENTREZID",#????????ID????
            OrgDb = "org.Hs.eg.db")
marker=cbind(rownames(diff),diff)
marker=marker[which(marker$`rownames(diff)` %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker$logFC
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
KEGG2 <- gseKEGG(gene_fc, organism = "hsa")
GO2 = gseGO(gene_fc,ont= "ALL",OrgDb = org.Hs.eg.db)
kk_gse <- GO2
kk_gse_cut = kk_gse[grep("collagen",kk_gse$Description)]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

#閫夋嫨灞曠幇NES鍓嶅嚑涓€氳矾 
up_gsea <- kk_gse_cut_up
pdf(file="转录组gsea.pdf",width=12,height=15.5)
gseaplot2(kk_gse,
          up_gsea$ID,#瀵岄泦鐨処D缂栧???
          color = "red",#GSEA绾挎潯棰滆壊
          base_size = 20,#鍩虹瀛椾綋澶у???
          rel_heights = c(1.5, 0.5, 1),#鍓浘鐨勭浉瀵归珮搴?
          subplots = 1:3, #瑕佹樉绀哄摢浜涘壇鍥? 濡俿ubplots=c(1,3) #鍙绗竴鍜岀涓変釜????
          ES_geom = "line",#enrichment score鐢ㄧ嚎杩樻槸鐢ㄧ???"dot"
) #鏄剧ずpvalue绛変???

dev.off()

##############################高低细胞通讯##########################
load("OS.Rdata")
load("onlyos.Rdata")
meta=pbmc@meta.data
osmeta = pbmc@meta.data
osh = osmeta[osmeta$os %in% "High OS",]
osl = osmeta[osmeta$os %in% "Low OS",]
meta[rownames(meta) %in% rownames(osh), "cell"] = "High Risk OS"
meta[rownames(meta) %in% rownames(osl), "cell"] = "Low Risk OS"

pbmc$newcell = meta$cell
DimPlot(pbmc,reduction="umap",group.by = "newcell")

load("OS.Rdata")
pbmc = subset(pbmc,newcell != "OS cells")
cellchat <- createCellChat(object=pbmc,group.by = "newcell")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 8)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可???
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

########
load("cellchat.Rdata")
par(mfrow = c(1,2), xpd=TRUE)

pdf(file="cellchat.pdf",width=6,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")

dev.off()
pdf(file="chat heatmap.pdf",width=12,height=6)
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
dev.off()

pathways.show <- c("COLLAGEN")  



cellchat@netP$pathways
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


levels(cellchat@idents)    # show all celltype
vertex.receiver = c(4,5) # define a numeric vector （淋系细胞）giving the index of the celltype as targets

pdf(file="hierarchy.pdf",width=10,height=10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")
dev.off()
pdf(file="chord.pdf",width=10,height=10)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
pdf(file="expression.pdf",width=10,height=10)
plotGeneExpression(cellchat, signaling = pathways.show)
dev.off()
pdf(file="bubble.pdf",width=5,height=8)
netVisual_bubble(cellchat, sources.use = c(4,5),signaling = pathways.show)
netVisual_bubble(cellchat, targets.use = c(4,5),sources.use = c(1,2,3,6,7,8,9),signaling = pathways.show)
dev.off()

######################?��?##############
load("OS.Rdata")

library(paletteer)
library(ggplot2)
pal librggsci)

DimPlot(pbmc, label = T, pF.size = 1,group.by = "cellsamp")+l_color_npg()+
  NoLe(x = "UMAP1", y = "UMAP2",title = "") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("all cell.tiff",width = 5,height = 5)ad("
lo
loodata")
pal <- paletteer_d("ggsci::nrc_npg")[c(1,2,5)]
pdf("os cell.pdf",width = 5,height = 5)
DimPlot(pbmc, label = T, pt.size = 1,cols = pal,group.by = "os")+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()
