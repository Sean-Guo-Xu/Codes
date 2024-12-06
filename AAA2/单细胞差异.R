library(Seurat)
setwd("E:\\AAA_scRNA")
load("8AAA.Rdata")
aaa=pbmc
load("6Normal.Rdata")
pbmc=merge(aaa,pbmc)

pbmc$sample = c(rep("AAA",length(aaa$orig.ident)),rep("Normal",length(pbmc$orig.ident)-length(aaa$orig.ident)))
pbmc = subset(pbmc,cell != "Erythrocyte" )

genelist = c("EIF2B2","CETP","MLH3","MRC2","LDAH","GDF7","NEK9","LRP1")
setwd("D:\\AAA\\qtl_scRNA")
for (gene in genelist) {
  if(gene %in% rownames(pbmc@assays$RNA)){
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
  ggsave(paste(gene,"_diff_vilin.tiff",sep=""),width=7,height = 4)
  
  data$cell
  p=ggplot(data, aes(x=cell,y=Expression,fill=sample)) +stat_compare_means( aes(label = ..p.signif..),
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
  ggsave(paste(gene,"all_diff_vilin.tiff",sep=""),p,width=3,height = 4)
  
  }}

setwd("E:\\AAA_scRNA\\Blood")
load("A_blood.Rdata")
aaa=pbmc
load("N_blood.Rdata")
pbmc=merge(aaa,pbmc)

pbmc$sample = c(rep("AAA",length(aaa$orig.ident)),rep("Normal",length(pbmc$orig.ident)-length(aaa$orig.ident)))
pbmc = subset(pbmc,cell != "Erythrocyte" )

genelist = c("EIF2B2","CETP","MLH3","NEK9","HMGCR","CERT1","LRP1")
setwd("D:\\AAA\\qtl_scRNA")
for (gene in genelist) {
  if(gene %in% rownames(pbmc@assays$RNA)){
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
    ggsave(paste(gene,"Blood_diff_vilin.tiff",sep=""),width=7,height = 4)
    
    data$cell
    p=ggplot(data, aes(x=cell,y=Expression,fill=sample)) +stat_compare_means( aes(label = ..p.signif..),
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
    ggsave(paste(gene,"Blood_all_diff_vilin.tiff",sep=""),p,width=3,height = 4)
    
  }}

setwd("E:\\AAA_scRNA\\Blood")
load("A_blood.Rdata")
  FeaturePlot(pbmc,features=genelist)

  
  genelist = c("EIF2B2","CETP","MLH3","MRC2","LDAH","GDF7","NEK9","LRP1")
  setwd("E:\\AAA_scRNA")
  load("8AAA.Rdata")
  FeaturePlot(pbmc,features=genelist)
  