library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(ggbiplot)
setwd("E:\\cervix")
anno=read.table("./GTEX/gencode.v23.annotation.gene.probemap",header = T)
gtex = read.table("./GTEX/gtex_RSEM_gene_tpm",header = T)
type=read.table("./GTEX/GTEX_phenotype",header = T,sep = "\t")
type = type[type$X_primary_site %in% c("Cervix Uteri","Uterus"),]
gtex[,2:ncol(gtex)] = 2^gtex[,2:ncol(gtex)]-0.001
gtex[,2:ncol(gtex)]=round(gtex[,2:ncol(gtex)],digits = 4)
gtex = gtex[order(gtex$sample),]
anno = anno[order(anno$id),]

library(limma)
gtex = gtex[,-1]
gtex=as.matrix(gtex)
rownames(gtex) = anno$gene
gtex = avereps(gtex)

load("TCGA-CESC.Rdata")
count <- assay(data,"tpm_unstrand")
see  =rowData(data)
rownames(count) =see$gene_name 
see = colData(data)
see = see$shortLetterCode
con = count[,see %in% "NT"]
tumor = count[,see %in% c("TP","TM")]
count = cbind(con,tumor)
count = count[order(rownames(count)),]
gtex = gtex[order(rownames(gtex)),]

count = count[rownames(count) %in% rownames(gtex), ]
gtex = gtex[rownames(gtex) %in% rownames(count), ]
colnames(gtex) = gsub("[.]","-",colnames(gtex))
count = avereps(count)
gtex = gtex[,colnames(gtex) %in% type$Sample]
group = cbind(c(colnames(gtex),colnames(count)),c(rep("GTEX",ncol(gtex)),rep("TCGA",ncol(count))))
exp=cbind(gtex,count)
exp = exp[rowSums(exp) != 0,]
exp = t(exp)
exp = log2(exp+1)
pca_result<- prcomp(exp,scale=T) 

ggbiplot(pca_result, 
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,         # 横纵比例 
         groups = group[,2], # 添加分组信息，将按指定的分组信息上色
         ellipse = T,           # 是否围绕分组画椭圆
         circle = F)
group=as.data.frame(group)
group$V3 = c(rep("Normal",ncol(con)+ncol(gtex)),rep("Tumor",ncol(tumor)))
design = model.matrix(~group$V3)
exp=round(exp,digits = 4)
rt = removeBatchEffect(t(exp),batch = group$V2,design = design)
rt=round(rt,digits = 4)

ggbiplot(prcomp(t(rt),scale=T) , 
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,         # 横纵比例 
         groups = group[,2], # 添加分组信息，将按指定的分组信息上色
         ellipse = T,           # 是否围绕分组画椭圆
         circle = F)
save(rt,file="combined.Rdata")
pdf("box.pdf")
boxplot(rt[,c(1:100,210:280)])
dev.off()

rt = normalizeBetweenArrays(rt)
fit=lmFit(rt,design)
fit=eBayes(fit)
options(digits = 4)
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
deg$gene = rownames(deg)
write.csv(deg,"deg.csv",row.names = F)
###########
library(dplyr)
library(ggrepel)
library(limma)
library(SummarizedExperiment)
library(ggsci)
mygene = "MICB"
DEG=deg
DEG<-DEG %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) 
ggplot(DEG,aes(logFC, -log10(P.Value)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.1, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#51B1B7", "#E1C855" ,"#E07B54"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG[DEG$gene %in% mygene,],color="black",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =mygene, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave("gtex_tcga_vol.tiff",width = 6,height = 5)
