library(Seurat)
setwd("D:\\AAA\\scRNA")
load("MMN.Rdata")
mmn=pbmc
load("7AAA.Rdata")
pbmc$cell[pbmc$cell %in% c("Macrophage","Monocyte","Neutrophil")]=mmn$newcell
save(pbmc,file = "7AAA.Rdata")

library(ggplot2)
library(ggpubr)
library(ggsci)
library(scales)
library(RColorBrewer)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
"#1F77B4FF" "#FF7F0EFF" "#2CA02CFF" "#D62728FF" "#9467BDFF" "#8C564BFF"
"#E377C2FF" "#7F7F7FFF" "#BCBD22FF" "#17BECFFF" "pink"    
pal = ggsci::pal_d3()(10)
pal = c(pal,"pink")
VlnPlot(pbmc,features = "PLEK",group.by = "cell",pt.size = 0.5,cols = pal)
ggsave("PLEK_vilin.tiff",width = 8,height = 4)
FeaturePlot(pbmc,features = "PLEK",order = T,cols = c("#1F77B4FF","#D62728FF"))+theme_void()
ggsave("PLEK_feature.tiff",width = 5,height = 5)
DimPlot(pbmc, label = T, pt.size = 1,cols = pal,group.by = "cell")+
  theme_void()+ NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "") 
ggsave("PLEK_UMAP.tiff",width = 5,height = 5)
dev.off()
load("N_A_MM.Rdata")
VlnPlot(pbmc,features = "PLEK",group.by = "cell",cols = c("#D62728FF","#1F77B4FF"),split.by = "sample")+stat_compare_means( aes(label = ..p.signif..), method = "t.test")+
  xlab("")
ggsave("PLEK_Macro_vln.tiff",width = 6,height = 4)

load("Macro.Rdata")
macro = pbmc
load("7AAA.Rdata")
pbmc=subset(pbmc,cell=="Monocyte")
macro$cell=macro$macrotype
mm=merge(pbmc,macro)
mm$sample = rep("AAA",length(mm$orig.ident))
unique(mm$cell)
load("all Macro.Rdata")
nor = subset(pbmc,sample=="Normal")
load("MMN normal.Rdata")
normo = subset(CONTROL,celltype=="Monocytes")
normo$cell = "Monocyte"
nor$cell = nor$celltype
nor$cell[nor$cell %in% "M0 Macrophages"]= "M0-like Macrophages"
nor$cell[nor$cell %in% "M1 Macrophages"]= "M1-like Macrophages"
nor$cell[nor$cell %in% "M2 Macrophages"]= "M2-like Macrophages"
nor = merge(nor,normo)
nor$sample = "Normal"
allmm = merge(mm,nor)
pbmc=allmm
save(pbmc,file="MM.Rdata")
DefaultAssay(allmm)="RNA"
FeaturePlot(pbmc,features = "PLEK")
Seurat::VlnPlot(pbmc,features = "PLEK",group.by = "cell")
ggsave("PLEK_feature.tiff",width = 5,height = 5)
load("7AAA.Rdata")
DefaultAssay(pbmc) = "RNA"
FeaturePlot(pbmc,features = "PLEK",order = T,cols = c("gray","#4DBBD5FF"))+theme_void()
FeaturePlot(aaa,features = "UGCG",order = T,cols = c("#4DBBD5FF","#F39B7FFF"))+theme_void()

ggsave("violin.tiff",width = 8,height = 5)
load("MM.Rdata")
#########����ͼ#############
exp = pbmc@assays$RNA@data
exp= as.matrix(exp)
meta = pbmc@meta.data
meta = cbind(meta,FetchData(pbmc,vars="rna_PLEK"))
meta$sample_cell=paste(meta$sample,meta$cell)
per = NULL
for (i in unique(meta$sample_cell)){
  count = nrow(meta[(meta$sample_cell %in% i) & (meta$rna_PLEK >1) ,]) / nrow(meta[(meta$sample_cell %in% i),])
  i=strsplit(i," ")
    per = rbind(per,c(i[[1]][1],i[[1]][2],count))
}
per=as.data.frame(per)
per[,3] = as.numeric(per[,3])
colnames(per)=c("sample","cell","percentage")
ggplot(per,aes(x=cell,y=percentage,
               fill=sample))+scale_color_npg()+geom_bar(stat="identity",position = "dodge",color = "black")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab(label="")+ylab(label="Percentage")
ggsave("ratio.tiff",width = 6,height = 5)

#########��ʱ��###############
library(monocle)
library(ggsci)
library(Seurat)
setwd("D:\\AAA\\scRNA")
load("MM.Rdata")
pbmc=subset(pbmc,sample=="AAA")
pbmc=subset(pbmc,cell != "Monocyte")
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
p_data$celltype <- pbmc$cell
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#��p_data��f_data��data.frameת��AnnotatedDataFrame����
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
pbmc@active.ident = as.factor(pbmc$cell)

marker=FindAllMarkers(pbmc)
cds <- setOrderingFilter(cds, unique(marker$gene))
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree',auto_param_selection = F)
cds <- orderCells(cds)
cds<-orderCells(cds,root_state = 11)
plot_ordering_genes(cds)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points = F,markers="AKT",markers_linear = T) 
load("Macro cds.Rdata")

 plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points = F,markers="UGCG",markers_linear = T,show_tree = T)+ scale_color_gradient(low = "#4DBBD5FF", high = "#E64B35FF")+
    theme(axis.ticks = element_blank(),legend.position = "right", axis.line.y = element_blank(),axis.line.x = element_blank(),axis.text = element_blank())+xlab("")+ylab("")
ggsave("UGCG_Time_trajectory.tiff",width = 6,height = 5)
 plot_cell_trajectory(cds,color_by="cell",show_branch_points = F,markers="UGCG",markers_linear = T,show_tree = T)+ scale_color_npg()+
   theme(axis.ticks = element_blank(),legend.position = "right", axis.line.y = element_blank(),axis.line.x = element_blank(),axis.text = element_blank())+xlab("")+ylab("")
ggsave("UGCG_Type_trajectory.tiff",width = 6,height = 5)
 
plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points = F,markers="PLEK",markers_linear = T,show_tree = T)+ scale_color_gradient(low = "#1F77B4FF", high = "#D62728FF")+
  theme(axis.ticks = element_blank(),legend.position = "right", axis.line.y = element_blank(),axis.line.x = element_blank(),axis.text = element_blank())+xlab("")+ylab("")
ggsave("PLEK_Time_trajectory.tiff",width = 6,height = 5)
plot_cell_trajectory(cds,color_by="cell",show_branch_points = F,markers="PLEK",markers_linear = T,show_tree = T)+ scale_color_manual(values=c("#2CA02CFF","#D62728FF","#1F77B4FF"))+
  theme(axis.ticks = element_blank(),legend.position = "right", axis.line.y = element_blank(),axis.line.x = element_blank(),axis.text = element_blank())+xlab("")+ylab("")
ggsave("PLEK_Type_trajectory.tiff",width = 6,height = 5)

plot_cell_trajectory(cds,color_by="cell",show_branch_points = F,show_tree = T)+  
  scale_size_continuous(range = c(1, 5))+geom_point(aes(size = sphingolipid,color=cell))+scale_color_npg()+theme(axis.ticks = element_blank(),legend.position = "right", axis.line.y = element_blank(),axis.line.x = element_blank(),axis.text = element_blank())+xlab("")+ylab("")
ggsave("sphingolipid_Type_trajectory.tiff",width = 6,height = 5)
plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points = F,show_tree = T)+  
  scale_size_continuous(range = c(1, 5))+geom_point(aes(size = sphingolipid,color=Pseudotime))+scale_color_gradient(low = "#4DBBD5FF", high = "#E64B35FF")+theme(axis.ticks = element_blank(),legend.position = "right", axis.line.y = element_blank(),axis.line.x = element_blank(),axis.text = element_blank())+xlab("")+ylab("")
ggsave("sphingolipid_time_trajectory.tiff",width = 6,height = 5)

############����ͼ############
library(Seurat)
setwd("D:\\AAA\\scRNA")
load("7AAA.Rdata")
aaa=pbmc
load("5Normal.Rdata")
pbmc=merge(aaa,pbmc)
pbmc$sample = c(rep("AAA",length(aaa$orig.ident)),rep("Normal",length(pbmc$orig.ident)-length(aaa$orig.ident)))
pbmc = subset(pbmc,cell != "Erythrocyte" )

genelist = c("")

for (gene in genelist) {
  

DefaultAssay(pbmc)="RNA"
data = cbind(pbmc$sample,FetchData(pbmc,paste("rna_",gene,sep="")),pbmc$cell,rep("Cell",length(pbmc$orig.ident)))
colnames(data)=c("sample","Expression","Celltype","cell")
library(ggpubr)
library(ggplot2)
library(ggsci)

ggplot(data, aes(x=Celltype, y=Expression,fill=sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                     method = "t.test")+
  geom_violin(trim=FALSE) + geom_boxplot(width=0.3,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_manual(values = c("#62B197","#E18E6D"))+ 
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显�?(size=1)
       
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  xlab("")
data$cell
ggplot(data, aes(x=cell,y=Expression,fill=sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                             method = "t.test")+
  geom_violin(trim=FALSE) + geom_boxplot(width=0.3,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_manual(values = c("#62B197","#E18E6D"))+ 
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=12), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显�?(size=1)
        
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  xlab("")
ggsave(paste(gene,"_diff_vilin.tiff"),width=3,height = 4)}
aaa$test[aaa$test %in% "1"] = "Macrophage"
aaa$test[aaa$test %in% "3"] = "DCs"
aaa$test[aaa$test %in% "2"] = "Neutrophil"
aaa$test[aaa$test %in% "4"] = "Monocyte"
################��ͼ����#############
library(ggplot2)
library(ggpubr)
library(ggsci)
library(scales)
library(RColorBrewer)

p = read.table("scmetabolism_pvalue.txt",sep = "\t")
path=read.table("scmetabolism_data.txt",sep="\t",header=T)
path=t(path)
path=as.data.frame(path)
path$sample = c(rep("AAA",924),rep("Normal",2801-924))
data=NULL
for (i in 1:85) {
  data=rbind(data,cbind(path[,i],path$sample,rep(colnames(path)[i],nrow(path))))
}
data=as.data.frame(data)
data$V1 = as.numeric(data$V1)
data=data[data$V3 %in% p$V1,]
colnames(data)=c("AUCscore","Sample","Pathway")
ggplot(data, aes(x=Pathway, y=AUCscore,fill=Sample)) +stat_compare_means( aes(label = ..p.signif..),
                                                                               method = "wilcox.test")+
  geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.8,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_npg()+ #设置填充的颜�?
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=10), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #设置y轴标题的字体属�?
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显�?(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属�?
                                 size=12),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属�?
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),legend.position = "left")+  #不显示网格线
  ylab("")+xlab("")
