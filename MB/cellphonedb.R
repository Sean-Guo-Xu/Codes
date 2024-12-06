setwd("E:\\brain\\cellphonedb")
load("E:\\brain\\cancer.Rdata")#这里直接读入我处理的pbmc3k的示例数据
pbmc = subset(pbmc, sample =="SHH")
pbmc = subset(pbmc,cell %in% c("Neuron", "Neuronal progenitor cell","Endothelial","Cycling cells","Pericytes"  ),invert = T)
counts <- as.matrix(pbmc@assays$RNA$counts)
write.table(counts,'cellphonedb_count.txt', row.names = T,col.names = T,sep='\t', quote=F)
meta_data <- cbind(rownames(pbmc@meta.data), pbmc@meta.data[,'cell', drop=F]) 
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #细胞类型不能为空
write.table(meta_data,'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

###############可视化##########
library(ktplots)
library(SingleCellExperiment)
library(Seurat)
library(pheatmap)


load("E:\\brain\\cancer.Rdata")
cancer = "SHH"
pbmc = subset(pbmc,sample ==cancer)
pbmc$celltype = as.factor(pbmc$cell)
setwd(paste("E:\\brain\\cellphonedb\\",cancer,sep=""))
filename = list.files()
pfile=filename[grep("statistical_analysis_pvalues",filename)]
mfile=filename[grep("statistical_analysis_means",filename)]
dfile=filename[grep("statistical_analysis_deconvoluted",filename)[1]]

pvals = read.table(pfile,header = T,check.names = F,sep = "\t")
means = read.table(mfile,header = T,check.names = F,sep = "\t")
decon = read.table(dfile,header = T,check.names = F,sep = "\t")
path = unique(pvals$classification)
path = path[path != ""]

p_values <- c(0.01, 0.02, 0.03, 0.04)

# 计算 Fisher's 统计量
integrate_p<-function(p_values){
fisher_stat <- -2 * sum(log(p_values))
combined_p_value <- 1 - pchisq(fisher_stat, df = 2 * length(p_values))
return(combined_p_value)
}
path_m = NULL
path_p = NULL
for (i in path) {
  partp = pvals[pvals$classification %in% i,]
  p = apply(partp[,14:ncol(partp)],2,integrate_p)
  partm = means[means$classification  %in% i,]
  m = apply(partm[,14:ncol(partm)],2,mean)
  path_m = rbind(path_m,m)
  path_p = rbind(path_p,p)
}
colnames(path_p) = names(p)
rownames(path_p) = path
colnames(path_m) = names(m)
rownames(path_m) = path
min = apply(path_p,1,min)
path_p = path_p[min<=0.05,]
path_m = path_m[min<=0.05,]
library(reshape2)
data_p = melt(path_p)
data_m = melt(path_m)
colnames(data_p) = c("Pathway","Celltype","P")
colnames(data_m) = c("Pathway","Celltype","Mean")
split_list <- strsplit(as.character(data_m$Celltype), "\\|")
celltype1 <- sapply(split_list, `[`, 1)
celltype2 <- sapply(split_list, `[`, 2)
data_m = cbind(data_m,celltype1,celltype2)
data_m$Celltype = gsub("\\|","->",data_m$Celltype)
data_m$P  = data_p$P
library(dplyr)

pathway_dotplot<-function(data,celltype1 = unique(data$celltype1),celltype2 = unique(data$celltype2)){
  library(ggplot2) 
  data  = data[data$celltype1 %in% celltype1,]
   data  = data[data$celltype2 %in% celltype2,]
   data <- data %>%
     group_by(Pathway) %>%
     mutate(scaled_Mean = (Mean - min(Mean)) / (max(Mean) - min(Mean)))
   
   ##wash
   for (i in unique(data$Celltype)) {
    if(min( data$P[data$Celltype %in% i] ) >= 0.05){
      data = data[data$Celltype != i, ] 
    }
    }
   for (i in unique(data$Pathway)) {
     if(min( data$P[data$Pathway %in% i] ) >= 0.05){
       data = data[data$Pathway != i, ] 
     }
   }
   data$Pathway = gsub("Adhesion by ","",data$Pathway)
   data$Pathway = gsub("Signaling by ","",data$Pathway)

   data$Celltype = as.factor(data$Celltype)
   data$size <- cut(data$P,
                    breaks = c(-Inf, 0.05, Inf),
                    labels = c( "p <= 0.05", "p > 0.05"))  # 点大小的相对值
 
   ggplot(data, aes(x = Celltype, y = Pathway)) +
     geom_point(aes(size = size, color = scaled_Mean)) +
     scale_size_manual(values = c("p <= 0.05" = 2,  "p > 0.05" = 1), 
                       guide = guide_legend(title = "p-value")) +  # 手动设置点大小
      scale_color_gradient2(low = "#5371b3",high= "#E31A1C", midpoint = 0.5) +  # 颜色渐变从冷色到暖色
     theme_minimal() +
     labs(x = "Celltype", y = "Pathway", title = "Scatter Plot of Celltype vs Pathway") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+ylab("")+ggtitle("") # 调整 x 轴文本角度
   
}
celltype= "SHH"
pathway_dotplot(data_m,celltype1 = celltype,celltype2 = unique(pbmc$cell))
ggsave(paste(celltype,"dotplot.png"),width = 6,height = 6)
celltype = c("APC","OPC")
pathway_dotplot(data_m,celltype1 =celltype ,celltype2 = unique(pbmc$cell))
ggsave(paste(celltype,"dotplot.png"),width = 10,height = 7)
celltype = c("cDC","pDC","Inf DC","Myeloid DC")
pathway_dotplot(data_m,celltype1 = celltype,celltype2 = unique(pbmc$cell))
ggsave(paste(celltype,"dotplot.png"),width = 12,height = 8)
celltype = c("NK cell","CD8+T cell","Naive CD4+T cell", "Memory CD4+T cell","Treg","B cell" )
pathway_dotplot(data_m,celltype1 =celltype ,celltype2 = unique(pbmc$cell))
ggsave(paste(celltype,"dotplot.png"),width = 18,height = 6)
celltype = c("M1"  ,"M2","Microglial")
pathway_dotplot(data_m,celltype1 =celltype ,celltype2 = unique(pbmc$cell))
ggsave(paste(celltype,"dotplot.png"),width = 10,height = 6)







pbmc = subset(pbmc, cell %in% c("Endothelial","Pericytes","Neuron","Neuronal progenitor cell","Cycling cells"),invert=T)

cancer = "SHH"
pdf(paste(cancer,"heatmap.pdf"),height = 6,width = 6)
plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10)
dev.off()

min = apply(pvals[,14:ncol(pvals)],1,min)
means = means[min <= 0.05,]
pvals = pvals[min <= 0.05,]
max = apply(means[,14:ncol(means)],1,max)
means = means[max >= 0.5,]
pvals = pvals[max >= 0.5,]
pbmc =  as.SingleCellExperiment(pbmc)
library(ggplot2)


for (i in unique(pbmc$cell)) {

p1=plot_cpdb(
  scdata=pbmc,
  cell_type1=cancer,
  cell_type2="CD8+T cell|Memory CD4+T cell|Naive CD4+T cell|NK cell|Treg",
  celltype_key="cell",
  means=means,
  pvals=pvals,
  #min_interaction_score = 2,
  scale_alpha_by_cellsign = T,
  max_size = 3,
  keep_significant_only=T
)+theme(legend.position = 'none',panel.grid.major =element_blank(),
        axis.text.x = element_text(angle = 90))
p2=plot_cpdb(
    scdata=pbmc,
    cell_type1=cancer,
    cell_type2="cDC|Myeloid DC|Inf DC|pDC",
    celltype_key="cell",
    means=means,
    pvals=pvals,
    #min_interaction_score = 2,
    scale_alpha_by_cellsign = T,
    max_size = 3,
    keep_significant_only=T
  )+theme(legend.position = 'none',panel.grid.major =element_blank(),
          axis.text.x = element_text(angle = 90))
p3=plot_cpdb(
    scdata=pbmc,
    cell_type1=cancer,
    cell_type2="M1|M2|Microglial",
    celltype_key="cell",
    means=means,
    max_size = 3,
    pvals=pvals,
    #min_interaction_score = 2,
    scale_alpha_by_cellsign = T,
    keep_significant_only=T
  )+theme(legend.position = 'none',panel.grid.major =element_blank(),
          axis.text.x = element_text(angle = 90))

  ggsave(paste(i,"dot.png"),height = 24,width = 20)
}

plot_cpdb(
  scdata=pbmc,
  cell_type1="M1|M2|Microglial",
  cell_type2="CD4+T cell|CD8+T cell|NK cells|Treg",
  celltype_key="cell",
  means=means,
  pvals=pvals,
 # min_interaction_score = 0.1
)
plot_cpdb(
  scdata=pbmc,
  cell_type1="CD4+T cell|CD8+T cell|NK cells|Treg",
  cell_type2="CD4+T cell|CD8+T cell|NK cells|Treg",
  celltype_key="cell",
  means=means,
  pvals=pvals,
  # min_interaction_score = 0.1
)
ggsave("T cell.png",height = 20,width = 18)

p1=plot_cpdb3(
  scdata = pbmc,
  cell_type1 = "B cell",
  cell_type2 = "SHH",
  celltype_key = "cell", # column name where the cell ids are located in the metadata
  means = means,
  pvals = pvals,
  deconvoluted = decon# new options from here on specific to plot_cpdb3
)

 receptor="B cell"
sendor ="Treg"
circlepic <- function(receptor,sendor,receptorcol="blue",sendorcol="red"){
plot_cpdb2(
  scdata = pbmc,
  cell_type1 = receptor,
  cell_type2 = ".",
  celltype_key = "cell", # column name where the cell ids are located in the metadata
  means = means,
  pvals = pvals,
  deconvoluted = decon, # new options from here on specific to plot_cpdb2
  desiredInteractions = list(
    c(sendor, receptor),
    c(receptor, sendor)
  ),
  interaction_grouping = interaction_annotation,
  edge_group_colors = c(
    "Activating" = "#e15759",
    "Chemotaxis" = "#59a14f",
    "Inhibitory" = "#4e79a7",
    "Intracellular trafficking" = "#9c755f",
    "DC_development" = "#B07aa1",
    "Unknown" = "#e7e7e7"
  ),
  node_group_colors = c(
  sendorcol,
  receptorcol
  ),
)
ggsave(paste(receptor,"circle.png"),width = 10,height = 10)
}
circlepic("B cell","Treg")
