rm(list=ls())

# Load
load(file='GSE7084_ID.Rdata')
newSet=raw_Set
dim(newSet)
colnames(newSet) <- paste(group_list,1:ncol(newSet),sep='_')

library(limma)
# design matrix
tmp=data.frame(case = c(rep(1,7), rep(0,8)), 
               control = c(rep(0,7), rep(1,8)))
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(newSet)
design
# contrast matrix
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                               levels = design)
contrast.matrix

# Check if it's necessary to normalize the data
library(RColorBrewer)
library(edgeR)
testing <- log2(newSet[,c(-7,-14,-15)]) # 只用了12个样本的数据，太多了不让用
dge <- DGEList(counts = testing)
col <- brewer.pal(ncol(dge$counts), "Paired")
par(mfrow=c(2,2))
# Un
boxplot(dge$counts,outline=F, col=col)
title(main="A. Unnormalised ",ylab="raw count")
# TMM
boxplot(calcNormFactors(dge, method = "TMM")$counts,outline=F,col=col)
title(main="B. TMM ",ylab="raw count")
# CPM
boxplot(cpm(dge$counts),outline=F, col=col)
title(main="C. CPM ",ylab="cpm")
# Log-CPM
boxplot(cpm(dge$counts,log=TRUE),outline=F, col=col)
title(main="D. Log-CPM ",ylab="log-cpm")
par(mfrow=c(1,1))
## 效果都不是很好，用Voom试一试

# using voom to normalize the data
dge <- DGEList(counts = newSet)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE) # 后续deg函数第一个参数用v
# 再试试刚刚的Log-CPM
# 画个完整的boxplot
## unnormalisied
par(mfrow=c(2,2))
dge <- DGEList(counts = newSet)
boxplot(dge$counts,outline=F, col=col)
title(main="Unnormalised ",ylab="log2count")
## normalized using Log-CPM
boxplot(cpm(dge$counts,log=TRUE),outline=F, col=col)
title(main="Log-CPM ",ylab="log-cpm")
cpm <- cpm(dge$counts,log=TRUE)

# define a DEG function
deg = function(newSet,design,contrast.matrix){
  #step1
  fit <- lmFit(newSet,design)
  #step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  
  fit2 <- eBayes(fit2)
  #eBayes() with trend=TRUE
  #step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}
re = deg(v,design,contrast.matrix) # 改第一个参数就行
nrDEG=re

## heatmap
library(pheatmap)
choose_gene=head(rownames(nrDEG),50) ## top 50 DEGs
choose_matrix=newSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename = 'DEG_top50_heatmap.png')

library(ggplot2)
## volcano plot
## logFC_cutoff is 1.10534
## p-value_cutoff is 0.05
colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))
DEG=nrDEG

logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs(logFC)) )

DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)

g = ggplot(data=DEG, 
           aes(x=logFC, y=-log10(P.Value), 
               color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
print(g)
ggsave(g,filename = 'volcano.png')

save(newSet,group_list,nrDEG,DEG, 
     file='GSE7084_DEG.Rdata')

load(file='GSE7084_DEG.Rdata')



# 得到上下调基因到DEGs
setwd("D:/学习/大创/韩彦槊老师/转录组学/R/AAA/GSE7084")
load(file='GSE7084_DEG.Rdata')
which(DEG$change!='NOT' ,arr.ind = T)
View(DEG[which(DEG$change!='NOT' ,arr.ind = T),])
DEGs_1 <- DEG[which(DEG$change!='NOT' ,arr.ind = T),]
write.csv(DEGs, file='voom_DEG.csv', quote = F)

DEG1 <- read.csv(file = 'voom_DEG.csv')
DEG1 <- DEG1[,1]

# GSE47472的:
setwd("D:/学习/大创/韩彦槊老师/转录组学/R/AAA/GSE47472")
load(file='GSE47472_DEG.Rdata')
which(DEG$change!='NOT' ,arr.ind = T)
View(DEG[which(DEG$change!='NOT' ,arr.ind = T),])
DEGs_2 <- DEG[which(DEG$change!='NOT' ,arr.ind = T),]
write.csv(DEGs, file='voom_DEG.csv', quote = F)

DEG2 <- read.csv(file = 'voom_DEG.csv')
DEG2 <- DEG2[,1]

# 得到各自up与down的degs
deg1up <- DEGs_1[which(DEGs_1$change=='UP' ,arr.ind = T),]
deg1down <- DEGs_1[which(DEGs_1$change=='DOWN' ,arr.ind = T),]
deg2up <- DEGs_2[which(DEGs_2$change=='UP' ,arr.ind = T),]
deg2down <- DEGs_2[which(DEGs_2$change=='DOWN' ,arr.ind = T),]
deg1up <- rownames(deg1up)
deg1down <- rownames(deg1down)
deg2up <- rownames(deg2up)
deg2down <- rownames(deg2down)


# 分别求up或down的交并集(保证两个数据集中都是up或都是down才行)
intered <- c(intersect(deg1up, deg2up), intersect(deg1down, deg2down))
uni <- union(DEG1, DEG2)
uni <- c(uni[!(uni %in% intersect(DEG1, DEG2))], intered)

# # 求并集
# uni <- union(DEG1, DEG2)
# 
# # 求交集
# intered <- intersect(DEG1, DEG2)

setwd("D:/学习/大创/韩彦槊老师/转录组学/R/AAA/GSE7084")
write.csv(intered, file = "intersected_DEGs.txt", quote = F, row.names = F)

# merge一下，得到uni和inter基因的完整数据
# inter 数据用的是deg后GSE7084的数据
inter_deg <- DEGs_1[which(rownames(DEGs_1) %in% intered),]
# uni 公共部分还是用的GSE7084，剩余不公共的，用的是各自的数据
uni_deg <- DEGs_1[which(rownames(DEGs_1) %in% uni),]
tmp <- uni[!(uni %in% rownames(uni_deg))]
tmp <- DEGs_2[which(rownames(DEGs_2) %in% tmp),]
uni_deg <- rbind(uni_deg, tmp)

save(uni_deg, inter_deg, file='uni_inter.Rdata')
load(file='uni_inter.Rdata')

# 画venn图
library(VennDiagram)
library(RColorBrewer)

venn.diagram(list(A=DEG1, B=DEG2),fill = c(brewer.pal(7,"Set1")[c(1,2)]),
             alpha = c(0.5, 0.5), cex = 2, 
             cat.cex=2, cat.fontface = 4,lty =2, fontfamily =3, resolution =300, 
             filename = "Venn.tiff")

DEG[rownames(DEG) %in% ,]
