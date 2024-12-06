setwd("E:\\bladder")
############TCGA##################
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
tcga <- GDCquery(project = "TCGA-BLCA" ,
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "STAR - Counts"
)
tcga = GDCdownload(tcga)
tcga <- GDCprepare(query = tcga,save = T,save.filename = "TCGA-BLCA.Rdata")
load("TCGA-BLCA.Rdata")
count <- assay(data,"unstranded")
see  =rowData(data)
rownames(count) = see$gene_name
see = colData(data)
see = see$shortLetterCode
N = count[,see %in% "NT"]
T = count[,see %in% "TP"]
count = cbind(N,T)
library(limma)
count = avereps(count)
save(count,file="TCGA_count.Rdata")
vmgene = read.table("vm.txt")
library(dplyr)
library(ggrepel)

library(ggsci)
list <- c(rep("Normal", 19), rep("BLCA",431-19)) %>% factor(., levels = c("BLCA", "Normal"), ordered = F)
list <- model.matrix(~factor(list)+0) 
colnames(list) = c("BLCA","Normal")

count<-voom(count,list,plot=TRUE)#
df.fit <- lmFit(count, list) 
df.matrix <- makeContrasts( BLCA - Normal , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
DEG <- topTable(fit,n = Inf, adjust = "fdr")
DEG$ID = rownames(DEG)
DEG = DEG[DEG$ID %in% vmgene$V1,]

DEG<-DEG %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) 
diffgene =list(TCGA = DEG[DEG$change %in% c("Down","Up"),])
ggplot(DEG,aes(logFC, -log10(P.Value )))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#0094CDFF", "#F39200FF" ,"#D51317FF"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =ID, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave("TCGA_Valcano.tiff",width = 6,height = 5)
#####################TCGA »»Õº###############
library(pheatmap)
count <- assay(data,"unstranded")
see  =rowData(data)
rownames(count) = see$gene_name
see = colData(data)
see = see$shortLetterCode
N = count[,see %in% "NT"]
T = count[,see %in% "TP"]
count = cbind(N,T)
count = avereps(count)
count = count[rownames(count) %in% vmgene$V1,]
count = count[order(rownames(count)),]
data=count
Type=c(rep("N",19),rep("T",431-19)) 
names(Type)=colnames(data)
Type=as.data.frame(Type)
data=log2(data+1)
library(ggsci)
pal_d3()(10)
typecol = c("#17BECFFF","#FF7F0EFF")
names(typecol)= c("N","T")
pdf(paste("TCGA_heatmap.pdf"),height=6,width=8)
pheatmap(data, 
         annotation=Type, 
         color = colorRampPalette(c("#0094CDFF", "white", "#D51317FF"))(50),
         labels_col = "",
         cluster_cols =F,
         fontsize = 8,
         fontsize_row=8,
         border_color = NA, 
         annotation_colors  = list(Type =  typecol),
         scale = "row",
         fontsize_col=5)
dev.off()

############GEO##################
geo = read.table("GSE13507_illumina_raw.txt",header=T,sep="\t")
colnames(geo) = gsub("GSM","",colnames(geo))

gene = read.table("GPL6102-11574.txt",header = T,sep = "\t",quote = "")
colnames(geo)[1] = "ID"
geo = merge(gene[,c(1,6)],geo,by="ID")
geo = as.matrix(geo)
rownames(geo) = geo[,2]
geo  = geo[,-c(1,2)]
geo = avereps(geo)

geo = rbind(colnames(geo),geo)
geo = geo[,order(geo[1,])]
geotype= c(rep("Control",10),rep("Surrounding",58),rep("Primary",165),rep("Recurrent",23))
geont = c(rep("Normal",68),rep("BLCA",188))
clinicgeo=cbind(geotype,geont)
geo=geo[-1,]


save(geo,clinicgeo,file="GEO-BLCA.Rdata")


library(dplyr)
library(ggrepel)
library(limma)
library(SummarizedExperiment)
library(ggsci)
setwd("E:\\bladder")
load("GEO-BLCA.Rdata") 
geo=as.data.frame(geo)
count = geo
save(count,clinicgeo,file="GEO-BLCA.Rdata")
count = count[,c(1:10,69:256)]
list <- c(rep("Normal", 10), rep("BLCA",188)) %>% factor(., levels = c("BLCA", "Normal"), ordered = F)
list <- model.matrix(~factor(list)+0) 
colnames(list) = c("BLCA","Normal")
name = rownames(count)
count = apply(count,2,as.numeric)
rownames(count) = name
for (i in rownames(count)) {
   if (min(count[i,]) < 0){
     count = count[!(rownames(count) %in% i),]
   }
}
count<-voom(count,list,plot=TRUE)#
df.fit <- lmFit(count, list) 
df.matrix <- makeContrasts( BLCA - Normal , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
DEG <- topTable(fit,n = Inf, adjust = "fdr")
DEG$gene = rownames(DEG)
DEG = DEG[DEG$gene %in% vmgene$V1,]

DEG<-DEG %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 0.25,
                                   ifelse(logFC > 0.25 ,'Up','Down'),'No change'))) 
diffgene$GEO =DEG[DEG$change %in% c("Down","Up"),]
ggplot(DEG,aes(logFC, -log10(P.Value)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#0094CDFF", "#F39200FF" ,"#D51317FF"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =gene, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave("GEO_Valcano.tiff",width = 6,height = 5)



######################GEO»»Õº##################
count = geo
count = apply(geo,2,as.numeric)
rownames(count) = rownames(geo)
for (i in rownames(count)) {
  if (min(count[i,]) < 0){
    count = count[!(rownames(count) %in% i),]
  }
}
count=count[rownames(count) %in% vmgene$V1,]

data=count
Type=c(rep("N",10),rep("S",58),rep("T",256-68)) 
names(Type)=colnames(data)
Type=as.data.frame(Type)
data=log2(data+1)
library(ggsci)
typecol = c("#17BECFFF","#2CA02CFF","#FF7F0EFF")
names(typecol)= c("N","S","T")
pdf(paste("GEO_heatmap.pdf"),height=6,width=6)
pheatmap(data, 
         annotation=Type, 
         color = colorRampPalette(c("#0094CDFF", "white", "#D51317FF"))(50),
         labels_col = "",
         cluster_cols =F,
         fontsize = 8,
         fontsize_row=8,
         border_color = NA, 
         annotation_colors  = list(Type =  typecol),
         scale = "row",
         fontsize_col=5)
dev.off()

save(diffgene,file="diffgene.Rdata")
