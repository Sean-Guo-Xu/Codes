library(Seurat)
setwd("D:\\AAA\\scRNA")
load("5Normal.Rdata")
normal = pbmc
load("7AAA.Rdata")
pbmc = merge(pbmc,normal)
pbmc1 = subset(pbmc,sample=="AAA")
pbmc2 = subset(pbmc,sample=="Normal")
pbmc$sample = substr(pbmc$sample,1,nchar(pbmc$sample)-1)
pbmc = list(pbmc1,pbmc2)
testAB.anchors <- FindIntegrationAnchors(object.list = pbmc, dims = 1:20)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
pbmc = testAB.integrated
control = subset(CONTROL,celltype %in% c("M0 Macrophages","M1 Macrophages","M2 Macrophages"))
pbmc$macrotype[pbmc$macrotype %in% "M0-like Macrophages"]="M0 Macrophages"
pbmc$macrotype[pbmc$macrotype %in% "M1-like Macrophages"]="M1 Macrophages"
pbmc$macrotype[pbmc$macrotype %in% "M2-like Macrophages"]="M2 Macrophages"
pbmc$celltype = pbmc$macrotype
sample=c(rep("AAA",length(pbmc$orig.ident)),rep("Normal",length(control$celltype)))
pbmc=merge(pbmc,control)
pbmc$sample =sample
save(pbmc,file = "all macro.Rdata")

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)                     #PCA??????????????????????????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????
ElbowPlot(pbmc)
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #????????????
pbmc <- FindClusters(object = pbmc, resolution = 0.6)                  #??????????,??????????????
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)   
DimPlot(pbmc,label = T)
DimPlot(pbmc,label = T,group.by = "state")

ggplot(data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=sphingolipid)
) + geom_point( size=1.5
) + scale_color_viridis(option="C")  + theme_light(base_size = 15)+labs(title = "")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

VlnPlot(pbmc,features = c('resident'),  pt.size = 0, adjust = 2,group.by = "seurat_clusters",combine = T)+scale_fill_npg()+theme(legend.position = "none")+stat_summary(fun.data = "mean_sdl",
                                                                                                                                                                        fun.args = list(mult = 1), geom = "pointrange", color = "white")

VlnPlot(pbmc,features = c('TREM2'),  pt.size = 0, adjust = 2,group.by = "seurat_clusters",combine = T)+scale_fill_npg()+theme(legend.position = "none")+stat_summary(fun.data = "mean_sdl",
                                                                                                                                                                     fun.args = list(mult = 1), geom = "pointrange", color = "white")
VlnPlot(pbmc,features = c('inflammatory'),  pt.size = 0, adjust = 2,group.by = "seurat_clusters",combine = T)+scale_fill_npg()+theme(legend.position = "none")+stat_summary(fun.data = "mean_sdl",
                                                                                                                                                                            fun.args = list(mult = 1), geom = "pointrange", color = "white")
pbmc$sample = c(rep("AAA",919),rep("Normal",6919-919))


load("all macro.Rdata")
pbmc$sample6 = paste(pbmc$sample,pbmc$celltype)
pbmc = ScaleData(pbmc)
DefaultAssay(pbmc) = "counts"
data1 =cbind(pbmc$sample,FetchData(pbmc,vars = "rna_SPHK1"),pbmc$celltype)
data2 =cbind(pbmc$sample,pbmc$apoptosis,rep("Apoptosis",length(pbmc$sample)))
data= data1
colnames(data) = c("Sample","Expression","Pathway")
data =as.data.frame(data)
data$AUC = as.numeric(data$Expression)
library(ggplot2)
library(ggpubr)
library(ggsci)
ggplot(data, aes(x=Pathway, y=Expression,fill=Sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                     method = "t.test")+
  geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_npg()+ #设置填充的颜�?
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显�?(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属�?
                                 size=12),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属�?
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("AUC")+xlab("")
pbmc$diffcell = paste(pbmc$celltype,pbmc$sample)
VlnPlot(pbmc,group.by = "diffcell",features = "EP300")
exp = pbmc@assays$RNA@data
exp= as.matrix(exp)
meta = pbmc@meta.data
meta = cbind(meta,FetchData(pbmc,vars="rna_SPHK1"))
per = NULL
for (i in unique(pbmc$sample6)){
  count = nrow(meta[(meta$sample6 %in% i) & (meta$rna_SPHK1 >2) ,]) / nrow(meta[(meta$sample6 %in% i),])
  per = rbind(per,c(i,count))
}
per=as.data.frame(per)
per$Celltype=c("M0","M1","M2","M0","M1","M2")
per$group=c(rep("AAA",3),rep("Normal",3))
per$V2=as.numeric(per$V2)
per1=per
#####
per = NULL
for (i in unique(pbmc$sample6)){
  count = nrow(meta[(meta$sample6 %in% i) & (meta$rna_SPHK1 > 0) &(meta$rna_SPHK1<2),]) / nrow(meta[(meta$sample6 %in% i),])
  per = rbind(per,c(i,count))
}
per=as.data.frame(per)
per$Celltype=c("M0","M1","M2","M0","M1","M2")
per$V2=as.numeric(per$V2)
per$group=c(rep("AAA",3),rep("Normal",3))
per2=per
#####
per = NULL
for (i in unique(pbmc$sample6)){
  count = nrow(meta[(meta$sample6 %in% i) & (meta$rna_SPHK1 <= 0) ,]) / nrow(meta[(meta$sample6 %in% i),])
  per = rbind(per,c(i,count))
}
per=as.data.frame(per)
per$Celltype=c("M0","M1","M2","M0","M1","M2")
per$V2=as.numeric(per$V2)
per$group=c(rep("AAA",3),rep("Normal",3))
per3=per
###
per= rbind(per3,per2,per1)
per$Expression =factor( c(rep("0",6),rep("1~2",6),rep(">2",6)),levels =c(">2","1~2","0" ))


ggplot(per,aes(x=group,y=V2,
               fill=Expression))+scale_fill_npg()+
  geom_bar(stat="identity",
                                    position = "stack")+facet_wrap(vars(Celltype))+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab(label="")+ylab(label="Percentage")

ggsave("SPHK1_percent.tiff",height = 5,width = 5)
#####�����л#####
library(Seurat)
library(clusterProfiler) 
library(AUCell)
library("KEGGREST")
setwd("D:\\AAA\\scRNA")
load("all macro.Rdata")
R.utils::setOption("clusterProfiler.download.method",'auto')
marker = FindMarkers(pbmc,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
marker[,1]=rownames(marker)
df_id<-bitr(rownames(marker),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
marker=marker[which(marker$p_val %in% df_id$SYMBOL),]
df_id=df_id[!duplicated(df_id$SYMBOL),]
colnames(marker)[1]="SYMBOL"
marker=merge(marker, df_id,by="SYMBOL")
gene_fc=marker[,2]
names(gene_fc)=marker$ENTREZID
gene_fc=gene_fc[order(gene_fc,decreasing = T)]
KEGG<-gseKEGG(gene_fc,organism   = 'hsa')
setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

geneset = list(LacticAcid=gene)
memory.limit(size=100000)
cells_rankings <- AUCell_buildRankings(pbmc@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- as.numeric(getAUC(cells_AUC))
pbmc$sphingolipid  <- aucs
