library("Seurat")
setwd("D:\\AAA\\scRNA")
load("7AAA.Rdata")
load("all.Rdata")
testAB.integrated = subset(testAB.integrated,sample == "AAA")


dir="D:\\AAA\\scRNA\\AAA\\1"
list.files(dir)
counts <- Read10X(data.dir = dir)
pbmc <- CreateSeuratObject(counts = counts)

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <2500 & percent.mt < 25)  
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
sample = rep("AAA1",length(pbmc$orig.ident))

all=list()
pbmc$sample=sample
all[[1]] = pbmc

for(i in 2:7){
  dir=paste("D:\\AAA\\scRNA\\AAA\\",i,sep="")
  counts = Read10X(data.dir = dir)
  pbmc=CreateSeuratObject(counts = counts)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
  pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <2500 & percent.mt < 25)  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  #提取那些在细胞间变异系数较大的基因
  pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
  
  pbmc$sample=rep(paste("AAA",i,sep = ""),length(pbmc$orig.ident))
  all[[i]]=pbmc
  
}
testAB.anchors <- FindIntegrationAnchors(object.list = all, dims = 1:20)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)


###########质量控制############
#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <2500 & percent.mt < 25)  


#测序深度的相关性绘图
pdf(file="04.featureCor.pdf",width=10,height=6)              #保存基因特征相关性图
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)              #保存基因特征方差图
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





###################################05.PCA主成分分析###################################
##PCA分析
pbmc=ScaleData(pbmc)                     #PCA降维之前的标准预处理步骤
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc)) 
#PCA分析
library(harmony)
pbmc = pbmc %>% RunHarmony("sample", plot_convergence = TRUE)#耗时1min

#细胞聚类
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:17) %>%
  FindNeighbors(reduction = "harmony", dims = 1:17) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
#绘制每个PCA成分的相关基因
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#主成分分析图形
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#主成分分析热图
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#每个PC的p值分布和均匀分布
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

ElbowPlot(pbmc, ndims=20, reduction="pca") 


###################################06.TSNE聚类分析和marker基因###################################
##TSNE聚类分析
pcSelect=12
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)                      #TSNE聚类


pdf(file="AAA cluste tsne.pdf",width=15,height=6)
DimPlot(nor,reduction="umap",group.by = c("cell"),label=T,repel=T) #TSNE可视化
DimPlot(pbmc,reduction="tsne",group.by = "deleteano",label=F,repel=T) +ggtitle ("Unknown")+ theme(legend.position = "none")#TSNE可视化
ggsave("1.tiff",height = 8,width = 9.5)
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)

##寻找差异表达的特征
logFCfilter=0.5
adjPvalFilter=0.05
pbmc=ifnb.list$Normal
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter
)
pbmc.markers=FindMarkers(pbmc,ident.1="Normal",ident.2="AML",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0.1)
write.table(pbmc.markers,file="07.normal-aml monocleMarkers.txt",sep="\t",row.names=T,quote=F)                 

sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="AAA.markers.xls",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = cluster10Marker) + NoLegend()
dev.off()

#绘制marker的小提琴图
pdf(file="06.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = c("IGLL5", "MBOAT1"))
dev.off()

#绘制marker在各个cluster的散点图
pdf(file="normal.markerScatter.pdf",width=10,height=8)
FeaturePlot(object = pbmc, features = c("HOXA9", "SORT1","SH3BP5"),cols = c("white", "brown3"),keep.scale = "all" )
dev.off()

#绘制marker在各个cluster的气泡图
pdf(file="normal markerBubble.pdf",width=12,height=6)
cluster10Marker=c("HOXA9", "SORT1","SH3BP5")
DotPlot(object = ref, features = cluster10Marker,group.by = c("labels"),cols = c("cyan3","brown3"),dot.scale = 6)
dev.off()

######手动注释######·
#肿瘤的
library(plyr)
label=as.character(nor$seurat_clusters)
label=revalue(label, c("0"="B cell", "1"="B cell"))
label=revalue(label, c("2"="NK/T cell", "3"="NK/T cell","4"="NK/T cell","11"="NK/T cell"))
label=revalue(label, c("9"="Fibroblast", "12"="Neutrophil","8"="Neutrophil","7"="Monocyte"))
label=revalue(label, c("6"="Smooth muscle cell", "17"="Smooth muscle cell"))
label=revalue(label, c("5"="Macrophage", "14"="Dendritic cell","13"="Endothelial cell"))
label=revalue(label, c("15"="Mast cell", "10"="Plasma cell","16"="Plasma cell"))
pbmc$cell=label
DimPlot(pbmc,reduction="umap",group.by = "cell",label=T,repel=T) #TSNE可视化

#肿瘤的
library(plyr)
label=nor$seurat_clusters
label=revalue(label, c("21"="B cell", "15"="Neutrophil"))
label=revalue(label, c("3"="Endothelial cell", "5"="Endothelial cell","12"="Endothelial cell","22"="Endothelial cell"))
label=revalue(label, c("0"="Fibroblast", "1"="Fibroblast","9"="Monocyte","14"="Monocyte"))
label=revalue(label, c("25"="Mast cell", "11"="Mast cell"))
label=revalue(label, c("23"="Dendritic cell", "13"="Plasma cell","2"="Smooth muscle cell","4"="Smooth muscle cell","7"="Smooth muscle cell"))
label=revalue(label, c("18"="Macrophage", "16"="Macrophage","6"="Macrophage"))
label=revalue(label, c("10"="NK/T cell", "20"="NK/T cell","27"="NK/T cell","17"="NK/T cell","19"="NK/T cell"))
label=revalue(label, c("8"="Erythrocyte", "26"="Erythrocyte","24"="Erythrocyte"))
nor@meta.data$cell=label
DimPlot(pbmc,reduction="tsne",group.by = "cell",label=T,repel=T) #TSNE可视化
#正常的
load("Normal 5.Rdata")
label=nor$seurat_clusters
label=revalue(label, c("0"="Monocyte", "1"="Macrophage"))
label=revalue(label, c("2"="T cell", "3"="Neutrophils","4"=" ","5"="Erythrocyte"))
label=revalue(label, c("6"="Smooth muscle cell", "7"="T cell","8"="Basophil","9"="Fibroblast"))
label=revalue(label, c("10"="Endothelial cell", "11"="NK cell"))
save(nor,file="Normal 5.Rdata")
DimPlot(nor,reduction="tsne",group.by = "cell",label=T,repel=T) #TSNE可视化
###################################07.注释细胞类型###################################
library(SingleR)
library(celldex)
reff=celldex::BlueprintEncodeData()
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
nor=pbmc
load("AAA.Rdata")
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = testAB.integrated, labels = testAB.integrated$cell)
nor@meta.data$sutolabels=cellpred$labels
ref$labels=cellpred$labels
singler = CreateSinglerObject(counts, annot = ann, "pbmc", min.genes = 0,
                              species = "Human", citation = "",
                              ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
                              reduce.file.size = T, numCores = 1)
singler$seurat = pbmc
singler$meta.data$xy = pbmc@reductions$tsne@cell.embeddings
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
write.table(clusterAnn,file="07.clusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="07.cellAnn.txt",quote=F,sep="\t",col.names=F)
GetAssayData(pbmc, "data")
hpca.se <- HumanPrimaryCellAtlasData()
#准备monocle分析需要的文件
pbmc=ifnb.list$AML
rm(ifnb.list)
rm(merged.ifnb)
memory.limit(size=300000)
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(pbmc$labels,file="07.monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(pbmc.markers,file="07.normal-aml monocleMarkers.txt",sep="\t",row.names=F,quote=F)

###################################08.样本合并分析###################################
library(ggplot2)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggsci)
library(scales)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(colorspace)
setwd("D:\\生信\\单细胞\\pic")  
load("D:/生信/单细胞/normal/GSM5936942_readCounts_TME_NL.rda")
pbmc = CreateSeuratObject(counts = readCounts,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
l=read.table("D:\\生信\\单细胞\\normal\\GSM5936942_meta_data_TME_NL.txt",sep="\t",header=T)
pbmc@meta.data$labels =l$CELLTYPE
new=l$ident_nl[-which(l$ident_nl %in% "Normal")]
new=c(l$ident_nl[which(l$ident_nl %in% "Normal")],rep("AML",length(new)))
pbmc@meta.data$sample = new 
merged.ifnb<-pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
merged.ifnb <- NormalizeData(merged.ifnb, normalization.method ="LogNormalize", scale.factor = 10000)
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #对数据进行过滤
pbmc$orig.ident=substr(pbmc$orig.ident,14,16)
merged.ifnb <- FindVariableFeatures(merged.ifnb, selection.method = "vst", nfeatures = 2000)
merged.ifnb <- ScaleData(merged.ifnb)
merged.ifnb <- RunPCA(merged.ifnb)
merged.ifnb <- FindNeighbors(merged.ifnb, dims = 1:30)
merged.ifnb <- FindClusters(merged.ifnb)
merged.ifnb <- RunUMAP(merged.ifnb, dims = 1:30)
merged.ifnb$orig.ident=substr(merged.ifnb$orig.ident,14,16)

p=DimPlot(merged.ifnb, reduction = "umap", group.by = "sample", repel = TRUE)
ggsave("D:\\normal$aml.tiff",p,width =7,height = 6.5)
p=DimPlot(merged.ifnb, reduction = "umap", group.by = "labels",label = T)
ggsave("D:\\cells.tiff",p,width =10.5,height = 7.5)
p=VlnPlot(merged.ifnb , features = c("HOXA9", "SORT1","SH3BP5"), split.by = "sample", group.by = "labels", pt.size = 0,ncol=1,cols = c("cyan3", "brown3"))
ggsave("D:\\markersexp.tiff",p,width =7.5,height = 7.5)
p=FeaturePlot(object =merged.ifnb , features = c("HOXA9", "SORT1","SH3BP5"),cols = c("gray", "blue"),keep.scale = "feature",split.by = "sample" )
ggsave("D:\\aml vs normal.tiff",p,width =6,height = 8)
cell.prop<-as.data.frame(prop.table(table(merged.ifnb$labels, merged.ifnb$orig.ident)))
colnames(cell.prop)<-c("celltype","sample","proportion")
p=ggplot(cell.prop,aes(sample,proportion,fill=celltype))+geom_bar(stat="identity",position="fill")+ggtitle("")+theme_bw()+guides(fill=guide_legend(title=NULL))+coord_flip()+scale_color_manual(colorRampPalette((pal_npg( "nrc")(9)))(24))
ggsave("D:\\cell proportion.tiff",p,width =10,height = 3.5)
ifnb.list <- SplitObject (merged.ifnb, split.by = "sample")
#比例图绘制
cell.prop$sample=substr(cell.prop$sample,1,2)
data=c("8","9",1.1)
for (i in unique(cell.prop$celltype)) {
  pp=cell.prop[which(cell.prop$sample %in% "NL"),]
  pp=pp[which(pp$celltype %in% i),3]
  pp=sum(pp)
  data=rbind(data,c(i,"NL",pp/length(merged.ifnb$sample[which(merged.ifnb$sample %in% "Normal")])*length(merged.ifnb$sample)))
}
data=data[-1,]

for (i in unique(cell.prop$celltype)) {
  pp=cell.prop[which(cell.prop$sample %in% "PT"),]
  
  pp=pp[which(pp$celltype %in% i),3]
  pp=sum(pp)
  data=rbind(data,c(i,"PT",pp/length(merged.ifnb$sample[which(merged.ifnb$sample %in% "AML")])*length(merged.ifnb$sample)))
}
data=as.data.frame(data)
data[,3]=as.numeric(data[,3])
p=ggplot() + geom_bar(data = data, aes(x = V1,y=V3,fill=V2),width = 0.6, stat = "identity", position = 'dodge')+theme_bw()+theme(legend.position = "bottom",legend.box = "horizontal",axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid()+xlab("")+ylab("")+scale_fill_manual(values = rep(c("brown3","cyan3"),24))
ggsave("D:\\PT&NL proportion.tiff",p,width =5,height = 3)

cluster10Marker=c("HOXA9", "SORT1","SH3BP5")
#kegg气泡图绘制
t=read.table("D:\\生信\\单细胞\\KEGG_Symbols.txt",sep=" ",fill=T,row.names=1)
logfc=read.table("C:\\Users\\56988\\Documents\\Tencent Files\\569886166\\FileRecv\\KEGG_Symbols_LogFC.txt",header=T)
name=0
for(i in unique(logfc$KEGG)){
  name=name+1
  kegg=logfc[which(logfc$KEGG %in% i),]
  kegg=kegg[order(abs(kegg$LogFC),decreasing = T),]
  kegg=kegg[1:5,]
  t[name,]=t[name,which(c(t[name,]) %in% kegg$Symbols)]
  
}
t=t[,1:5]
for (i in 1:7) {
  gene=c( t[i,which(t[i,] %in% merged.ifnb@assays$RNA@data@Dimnames[[1]])])
  gene=as.vector(unlist(gene))
  
  data=c("x","y",1,1)
  p=VlnPlot(ifnb.list$AML,features = gene,group.by = "labels")
  for (j in 1:length(gene)) {
    allcell=p[[j]]$data
    data=rbind(data,c(gene[j],paste("AML",rownames(t)[i]),sum(allcell[,1])/length(allcell[,1]),length(allcell[allcell[,1]>0,2])/length(allcell[,2])))
  }
  data=data[-1,]
  p=VlnPlot(ifnb.list$Normal,features = gene,group.by = "labels")
  for (j in 1:length(gene)) {
    allcell=p[[j]]$data
    data=rbind(data,c(gene[j],paste("Normal",rownames(t)[i]),sum(allcell[,1])/length(allcell[,1]),length(allcell[allcell[,1]>0,2])/length(allcell[,2])))
  }
  data=as.data.frame(data)
  data$V4=as.numeric(data$V4)*100
  data$V3=as.numeric(data$V3)
  p=ggplot(data,aes(x=V1,y=V2,fill=V3))+geom_point(aes(size=V4,color=V3))+theme_bw()+theme(legend.box = "horizontal",axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid()+xlab("")+ylab("")+scale_color_gradient(low="cyan3",high="brown3")
  ggsave(paste("D:\\AMLfunction",i,".tiff",sep=""),p,width =6,height = 2)
}


marker=read.table("07.normal-aml monocleMarkers.txt",sep="\t",header=T,check.names=F)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(snow)
R.utils::setOption("clusterProfiler.download.method",'auto') 
df_id<-bitr(rownames(marker), #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")
marker=cbind(rownames(marker),marker)
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker$avg_log2FC
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
KEGG <- gseKEGG(gene_fc, organism = "hsa")
paths <- c("hsa04940", "hsa04914", "hsa04114", "hsa04672","hsa05332","hsa04110","hsa05330")#选取你需要展示的通路ID
p=gseaplot2(KEGG,paths,ES_geom = "dot")
ggsave(paste("D:\\Normal_AML_GSEA.tiff",sep=""),p,width =7,height = 7)
###############合并###########
setwd("D:\\AAA\\scRNA")
load("5Normal.Rdata")
load("7AAA.Rdata")
pbmc = merge(nor,pbmc)
label = c(rep("Normal",37207),rep("AAA",75946-37207))
pbmc$allsample = label

DefaultAssay(pbmc) = "RNA"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #??????????,???????????榛?
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)   
library(harmony)
pbmc = pbmc %>% RunHarmony("allsample", plot_convergence = TRUE)#鑰楁椂1min
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:pcSelect) %>%
  FindNeighbors(reduction = "harmony", dims = 1:pcSelect) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
DimPlot(pbmc,reduction="umap",group.by = "allsample",split.by = "cell")
##################AUCELL#####
library(Seurat)
library(clusterProfiler) 
library(AUCell)
library("EnrichmentBrowser")
library("KEGGREST")
setwd("D:\\AAA\\scRNA")
load("all Macro.Rdata")
path <- keggGet("hsa00600")
gene<- path[[1]]$GENE
gene <- unlist(lapply(gene,function(x) strsplit(x,";")))
gene <- gene[1:length(gene)%%3 == 2]
geneset = list(sphingolipid=gene)
memory.limit(size=100000)
cells_rankings <- AUCell_buildRankings(pbmc@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- as.numeric(getAUC(cells_AUC))
pbmc$sphingolipid  <- aucs
library(ggraph)
library(paletteer)
library(ggsci)

ggplot(data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=apoptosis)
) + geom_point( size=1.5
) + scale_color_viridis(option="C")  + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave("AUC apoptosis.tiff",width = 6,height = 5)
ggplot(data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=sphingolipid)
) + geom_point( size=1.5
) + scale_color_viridis(option="C")  + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave("AUC sphingolipid.tiff",width = 6,height = 5)
pal <- paletteer_d("ggsci::nrc_npg")[c(1,8,3,2,5,6,7,10)]

DimPlot(pbmc, label = T, pt.size = 1,group.by = "EC")+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "") +scale_color_npg()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("EC cell.tiff",width = 5,height = 5)
save(pbmc,file="5Normal.Rdata")
#####################GSEA#######################
library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
R.utils::setOption("clusterProfiler.download.method",'auto') 
setwd("D:\\AAA\\scRNA")
load("5Normal.Rdata")
nor = pbmc
load("7AAA.Rdata")
pbmc = merge(nor,pbmc)
label = c(rep("Normal",37207),rep("AAA",75946-37207))
pbmc$allsample = label
DefaultAssay(pbmc) = "RNA"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
memory.limit(size = 1000000)
pbmc = SplitObject (pbmc, split.by = "cell")
pbmc["Erythrocyte"] = NULL
for (i in pbmc){
  pbmc.markers=FindMarkers(i,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "allsample",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
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
  KEGG <- gseKEGG(gene_fc, organism = "hsa")
  
  
  re = KEGG@result
  write.table(re,paste(unique(i$cell),"gsea.txt"),sep = "\t",quote = F,row.names = F)
}
########################EC#############
library(Seurat)
setwd("D:\\AAA\\scRNA")
load("5Normal.Rdata")
nor = pbmc
load("7AAA.Rdata")
pbmc = SplitObject (pbmc, split.by = "cell")
nor = SplitObject (nor, split.by = "cell")
pbmc = merge(pbmc$`Endothelial cell`,nor$`Endothelial cell`)
rm(nor)
DefaultAssay(pbmc) = "RNA"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 10
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #??????????,???????????姒??
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)
pbmc$sample = c(rep("AAA",919),rep("Normal",6919-919))
data1 =cbind(pbmc$sample,pbmc$sphingolipid,rep("Sphingolipid",length(pbmc$sample)))
data2 =cbind(pbmc$sample,pbmc$apoptosis,rep("Apoptosis",length(pbmc$sample)))
data= rbind(data1,data2)

colnames(data) = c("Sample","AUC","Pathway")
data =as.data.frame(data)
data$AUC = as.numeric(data$AUC)
library(ggplot2)
library(ggpubr)
ggplot(data, aes(x=Pathway, y=AUC,fill=Sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                     method = "t.test")+
  geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9))+ #缁樺埗绠辩嚎鍥?
  scale_fill_npg()+ #璁剧疆濉厖鐨勯鑹?
  theme_bw()+ #鑳屾櫙鍙樹负鐧借壊
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #璁剧疆x杞村埢搴︽爣绛剧殑瀛椾綋鏄剧ず鍊炬枩瑙掑害涓?15搴︼紝骞跺悜涓嬭皟鏁?1(hjust = 1)锛屽瓧浣撶皣涓篢imes澶у皬涓?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #璁剧疆y杞村埢搴︽爣绛剧殑瀛椾綋绨囷紝瀛椾綋澶у皬锛屽瓧浣撴牱寮忎负plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #璁剧疆y杞存爣棰樼殑瀛椾綋灞炴�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #鍘婚櫎榛樿濉厖鐨勭伆鑹诧紝骞跺皢x=0杞村拰y=0杞村姞绮楁樉绀?(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #璁剧疆鍥句緥鐨勫瓙鏍囬鐨勫瓧浣撳睘鎬?
                                 size=12),
        legend.title=element_text(face="italic", family="Times", colour="black", #璁剧疆鍥句緥鐨勬�绘爣棰樼殑瀛椾綋灞炴�?
                                  size=12),
        panel.grid.major = element_blank(),   #涓嶆樉绀虹綉鏍肩嚎
        panel.grid.minor = element_blank())+  #涓嶆樉绀虹綉鏍肩嚎
  ylab("AUC")+xlab("")
ggsave("AUC_AAA_Nor.tiff",width = 4,height = 4)
data = cbind(pbmc$sphingolipid,pbmc$apoptosis,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Sphingolipid","Apoptosis","Sample")
data$Sphingolipid = as.numeric(data$Sphingolipid)
data$Apoptosis = as.numeric(data$Apoptosis)
ggplot(data, aes(x=Sphingolipid, y=Apoptosis, color=Sample)) + scale_color_manual(values = c("brown3","cyan3"))+
  geom_point()+theme_bw()+ stat_smooth(method=lm,formula=y~x)+stat_cor(aes(x =Sphingolipid, y =Apoptosis))
ggsave("col_AAA_Nor.tiff",width = 6,height = 5)
########################轨迹#####################
library(Seurat)
library(monocle)
library(ggpubr)
library("ggsci")
setwd("D:\\AAA\\scRNA")
load("7AAA.Rdata")
pbmc = subset(pbmc,cell =="Endothelial cell" )
DefaultAssay(pbmc) = "RNA"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 15
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #??????????,???????????姒??
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)  
DimPlot(pbmc,reduction="umap",group.by = "seurat_clusters",split.by = "cell")
library(ggraph)
ggplot(data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=sphingolipid)
) + geom_point( size=1.5
) + scale_color_viridis(option="C")  + theme_light(base_size = 15)+labs(title = "TNFA_SIGNALING_VIA_NFKB")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(pbmc,features = c("apoptosis"),  pt.size = 0, adjust = 2,group.by = "EC",combine = T)+scale_fill_npg()+stat_summary(fun.data = "mean_sdl",
             fun.args = list(mult = 1), geom = "pointrange", color = "white")+theme(legend.position = "none")
ggsave("EC_ap_vol.tiff",width = 6,height = 4)
VlnPlot(pbmc,features = c('sphingolipid'),  pt.size = 0, adjust = 2,group.by = "EC",combine = T)+scale_fill_npg()+theme(legend.position = "none")+stat_summary(fun.data = "mean_sdl",
                                                                                                                                            fun.args = list(mult = 1), geom = "pointrange", color = "white")
ggsave("EC_sp_vol.tiff",width = 6,height = 4)

logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter
)
write.table(pbmc.markers,file="EC Markers.txt",sep="\t",row.names=T,quote=F) 

monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="EC.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="EC.monocleGene.txt",quote=F,sep="\t",row.names=F)


monocle.sample=read.table("EC.monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("EC.monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("EC Markers.txt",sep="\t",header=T,check.names=F)

library("dplyr")

#灏哠eurat缁撴灉杞崲涓簃onocle闇�瑕佺殑缁嗚優鐭╅樀锛岀粏鑳炴敞閲婅〃鍜屽熀鍥犳敞閲婅〃琛?
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)

cds <- estimateSizeFactors(cds)
memory.limit(size=300000000000)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree',auto_param_selection = F)
cds <- orderCells(cds)
cds <- orderCells(cds,reverse = T)
data = cbind(cds$sphingolipid,cds$Pseudotime,cds$seurat_clusters)
data = as.data.frame(data)
colnames(data) = c("sphingolipid","Pseudotime","Cluster")
data$Cluster = as.factor(data$Cluster)
pbmc$seurat_clusters = as.factor(pbmc$seurat_clusters)

p1=ggplot(data, aes(x=Pseudotime, y=sphingolipid,color= Cluster))+geom_point()+scale_color_npg()+theme_bw()+ 
  stat_smooth(aes(x=Pseudotime, y=sphingolipid,color=c("Line")),data,method=lm,formula=y~x,color="black")+stat_cor(aes(x =Pseudotime, y =sphingolipid,color =c("Line")),color = "black")

data = cbind(cds$apoptosis,cds$Pseudotime,cds$seurat_clusters)
data = as.data.frame(data)
colnames(data) = c("apoptosis","Pseudotime","Cluster")
data$Cluster = as.factor(data$Cluster)
pbmc$seurat_clusters = as.factor(pbmc$seurat_clusters)
p2=ggplot(data, aes(x=Pseudotime, y=apoptosis,color= Cluster))+geom_point()+scale_color_npg()+theme_bw()+ 
  stat_smooth(aes(x=Pseudotime, y=apoptosis,color=c("Line")),data,method=lm,formula=y~x,color="black")+stat_cor(aes(x =Pseudotime, y =apoptosis,color =c("Line")),color = "black")
p1+p2
ggsave("EC_ap_dot.tiff",width = 10,height = 4)
load(file="ECcds.Rdata")
################################trajectory################
p1=plot_cell_trajectory(cds,color_by = "apoptosis",show_branch_points = F)+ scale_color_material("teal")
p2=plot_cell_trajectory(cds,color_by = "sphingolipid",show_branch_points = F)+ scale_color_material("cyan")

p3=plot_cell_trajectory(cds,color_by = "EC",show_branch_points = F)+scale_color_npg()
p4=plot_cell_trajectory(cds,color_by = "Pseudotime",show_branch_points = F)
p1+p2+p3+p4
ggsave("pseudotime.tiff",width=10,height=10)
pdf(file="check.trajectory.pdf",width=5,height=5)
plot_complex_cell_trajectory(cds, x = 1, y = 2,
                             color_by = "apoptosis",show_branch_points =T)+
  theme(legend.title = element_blank()) 

plot_cell_trajectory(cds,color_by = "os")
dev.off()
library(ggplot2)
library(scales)
test_genes=c("UGCG","SGSM2")
pdf(file="gene time.pdf",width=6,height=4)
plot_genes_in_pseudotime(cds[test_genes,],color_by = "seurat_clusters")
dev.off()

pData(cds)$UGCG = log2( exprs(cds)['UGCG',]+1)
plot_cell_trajectory(cds, color_by = "UGCG",show_branch_points =F)  +scale_colour_gradient2(low = muted("cyan3"),
                                                                        
                                                                        mid = "cyan3",
                                                                        
                                                                        high = muted("brown3"))
ggsave(file="UGCG trajectory.tiff",width=5,height=5)
pData(cds)$PSAP = log2( exprs(cds)['PSAP',]+1)
plot_cell_trajectory(cds, color_by = "PSAP",show_branch_points =F)  +scale_colour_gradient2(low = muted("cyan3"),
                                                                                            
                                                                                            mid = "cyan3",
                                                                                            
                                                                                            high = muted("brown3"))
ggsave(file="PSAP trajectory.tiff",width=5,height=5)
pData(cds)$SGSM2 = log2( exprs(cds)['SGSM2',]+1)
plot_cell_trajectory(cds, color_by = "SGSM2",show_branch_points =F)  +scale_colour_gradient2(
  
  mid = "cyan3",
  
  high = muted("brown3"))
ggsave(file="SGSM2 trajectory.tiff",width=5,height=5)
ordergene = marker$gene
##################分支点################
BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
save(BEAM_res,file="ECbean.Rdata")
pdf("point heatmap.pdf",width = 8,height = 15)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 3, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 6,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T)
dev.off()
clusters <- cutree(p$ph$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
p$ph$tree_row$order
gene=BEAM_res[BEAM_res$gene_short_name %in% rownames(clustering),]  
gene = cbind(gene,clustering)
badfate = gene[gene$Gene_Clusters %in% c("1","2"),]

goodfate = gene[gene$Gene_Clusters %in% c("3","4"),]

save(badfate,goodfate,file="fate DEgenes.Rdata")
##cell fate 2 is bad and 1 is good### cluster 1,2 is bad 3,4 is good####
############取交集############
library("EnrichmentBrowser")
library("KEGGREST")
setwd("D:\\AAA\\scRNA")
load("fate DEgenes.Rdata")
path <- keggGet("hsa04210") #hsa04210
gene<- path[[1]]$GENE
gene <- unlist(lapply(gene,function(x) strsplit(x,";")))
gene <- gene[1:length(gene)%%3 == 2]
goodfate[goodfate$gene_short_name %in% gene,]

badfate[badfate$gene_short_name %in% gene,]

library(Seurat)
load("ECseurat.Rdata")
marker = FindMarkers(pbmc,ident.1 = 2,ident.2 = 3,group.by = "seurat_clusters")
marker = marker[marker$p_val_adj <0.05,]
see = marker[rownames(marker) %in% gene,]

###################通讯################3
library(CellChat)
setwd("D:\\AAA\\scRNA")
load("ECseurat.Rdata")
cell = pbmc$seurat_clusters
cell = cell[cell %in% c(2,3)]
load("7AAA.Rdata")
meta=pbmc@meta.data

meta[rownames(meta) %in% names(cell)[cell %in% 2] , "cell"] = "High EC"
meta[rownames(meta) %in% names(cell)[cell %in% 3], "cell"] = "Low EC"

pbmc$newcell = meta$cell
DimPlot(pbmc,reduction="umap",group.by = "newcell")
save(pbmc,file="7AAA.Rdata")

pbmc = subset(pbmc,newcell != "Endothelial cell")
DefaultAssay(pbmc) = "RNA"
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
save(cellchat,file="cellchat.Rdata")

########
library(CellChat)
library(ggpubr)
library("ggsci")
library(scales)
show_col(pal_npg("nrc")(10))
setwd("D:\\AAA\\scRNA")
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
ap = c("IL4","TNF","TRAIL","FASLG","NGF")
sg = c("PSAP")
pathways.show <- ap[5]



cellchat@netP$pathways
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


levels(cellchat@idents)    # show all celltype
vertex.receiver = c(4,5) # define a numeric vector （淋系细胞）giving the index of the celltype as targets

pdf(file=paste(pathways.show,"hierarchy.pdf"),width=14,height=7)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")
dev.off()


plotGeneExpression(cellchat, signaling = pathways.show)
ggsave(file=paste(pathways.show,"expr.tiff"),width=10,height=2.5)
w=8
h=2
netVisual_bubble(cellchat, sources.use = c(4,5),signaling = pathways.show)
ggsave(file=paste(pathways.show,"bubble 1.tiff"),width=w,height=h)
netVisual_bubble(cellchat, targets.use = c(4,5),sources.use = c(1,2,3,6,7,8,9,10,11,12),signaling = pathways.show)
ggsave(file=paste(pathways.show,"bubble 2.tiff"),width=w,height=h)
dev.off()


library(NMF)
library(ggalluvial)
nPatterns = 5 
# 同样在这里遇到了bug，难道说是我没有安装好吗，de了它。
# cellchat <- myidentifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)  
pdf(file="pattrern.pdf",width=10,height=12)
cellchat <- myidentifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()

###############周期##########
library("tricycle")
library("Seurat")
library(ggplot2)
library(scattermore)
library(scater)
library("ggsci")
library(scales)
setwd("D:\\AAA\\scRNA")
load("ECseurat.Rdata")
neurosphere_example<- as.SingleCellExperiment(pbmc,assay = "RNA")
gocc_sce.o <- run_pca_cc_genes(neurosphere_example,
                               exprs_values = "logcounts", species = "human",gname.type = "SYMBOL")

new.ref <- attr(reducedDim(gocc_sce.o, 'PCA'), 'rotation')[, seq_len(2)]
neurosphere_example <- project_cycle_space(neurosphere_example,ref.m = new.ref,gname.type = "SYMBOL",species = "human")
neurosphere_example <- estimate_cycle_position(neurosphere_example, ref.m  = new.ref,
                                   dimred = 'tricycleEmbedding2')

neurosphere_example <- estimate_cycle_position(neurosphere_example,g)
neurosphere_example <- estimate_Schwabe_stage(neurosphere_example,
                                              gname.type = 'SYMBOL',
                                              species = 'human')

fit.l <- fit_periodic_loess(neurosphere_example$tricyclePosition,
                            assay(neurosphere_example, 'logcounts')["CDC16",],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(neurosphere_example), ")"))
fit.l$fig + theme_bw(base_size = 14)
see = pbmc$EC
see[see %in% c("EC1","EC2","EC5","EC6","EC7")] = "Others"
neurosphere_example$draw = see
plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$draw, 'EC type',
                    bw = 10, fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)+scale_color_npg()


plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$draw, 'EC type', type = "circular",
                    bw = 10,  fig.title = "Kernel density of \u03b8",line.size=1) +
  theme_bw(base_size = 14)
ggsave("cycle.tiff",width = 6,height = 5)
library(cowplot)
p <- plot_emb_circle_scale(neurosphere_example, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave("EC cycle.tiff",width = 7,height = 5)

EC1 = neurosphere_example[,neurosphere_example$EC =="EC4"]
p <- plot_emb_circle_scale(EC1, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave("EC4 cycle.tiff",width = 7,height = 5)

EC1 = neurosphere_example[,neurosphere_example$EC =="EC3"]
p <- plot_emb_circle_scale(EC1, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave("EC3 cycle.tiff",width = 7,height = 5)

library(ggpubr)

data = cbind(pbmc$sphingolipid,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Sphingolipid","Cycle","Sample")
data$Sphingolipid = as.numeric(data$Sphingolipid)
data$Cycle = as.numeric(data$Cycle)
ggplot(data, aes(y=Sphingolipid, x=Cycle)) +
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")
ggsave(paste("Sphingolipid","cycle.tiff"),width=5,height = 5)

data = cbind(pbmc$apoptosis,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Apoptosis","Cycle","Sample")
data$Apoptosis = as.numeric(data$Apoptosis)
data$Cycle = as.numeric(data$Cycle)
ggplot(data, aes(y=Apoptosis, x=Cycle)) +
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")+ scale_x_continuous(breaks =seq(0, 2*pi, pi),labels=c('0', expression(pi), expression(2*pi)))
ggsave(paste("Apoptosis","cycle.tiff"),width=5,height = 5)

data = cbind(cds$Pseudotime,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Apoptosis","Cycle","Sample")
data$Apoptosis = as.numeric(data$Apoptosis)
data$Cycle = as.numeric(data$Cycle)
ggplot(data, aes(y=Apoptosis, x=Cycle)) +
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")+ scale_x_continuous(breaks =seq(0, 2*pi, pi),labels=c('0', expression(pi), expression(2*pi)))



load(file="ECcds.Rdata")
library(monocle)
for(gene in c("PSAP","SGSM2","UGCG")){
pbmc$SGSM2 = log2( exprs(cds)[gene,]+1)
data = cbind(pbmc$SGSM2,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("SGSM2","Position","Sample")
data$SGSM2 = as.numeric(data$SGSM2)
data$Position = as.numeric(data$Position)
ggplot(data, aes(y=SGSM2, x=Position)) +ylab(gene)+xlab("Cycle")+
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")+ scale_x_continuous(breaks =seq(0, 2*pi, pi),labels=c('0', expression(pi), expression(2*pi)))
ggsave(paste(gene,"cycle.tiff"),width=5,height = 5)}
