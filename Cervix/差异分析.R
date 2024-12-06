library(affy)
library(GEOquery)
library(affyPLM)
library(dplyr)
setwd("E:\\cervix")
mygene="HDAC10"
gse="GSE9750"
gset <- getGEO(gse,getGPL = F,destdir = ".")
data.raw <- ReadAffy(celfile.path = paste("./",gse,"_raw",sep = ""))
sampleNames(data.raw) <- sub(pattern = "\\.CEL.gz",replacement = "",sampleNames(data.raw))
data.rma <- rma(data.raw)
data.eset <- exprs(data.rma)
pd <- pData(gset[[1]])
pd <- pd %>% 
  select(geo_accession,title)
colnames(pd) <- c("id","sample")
pd <- pd[order(pd$sample),]
data.eset <- data.eset[,pd$id]

pd$sample[1:42]="CC"
pd$sample[43:66]="Normal"
group_list <- factor(pd$sample,levels = c("Normal","CC"))

names(group_list) <- pd$id
design <- model.matrix(~group_list)

library(limma)
fit <- lmFit(data.eset,design)
fit1 <- eBayes(fit)
options(options=4)
deg <- topTable(fit1,coef=2,n=Inf) 
cc=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "CC"]])
normal=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "Normal"]])
cc - normal  
deg = cbind(rownames(deg),deg)
colnames(deg)[1]="probe_id"
vm = read.table("vm.txt")
gpl <- gset[[1]]@annotation
library(AnnoProbe)
library(ggplot2)
library(ggrepel)
data.eset=data.eset[order(rownames(data.eset)),]
probe2id <- idmap(gpl)
probe2id=probe2id[order(probe2id$probe_id),]
data.eset = data.eset[rownames(data.eset) %in% probe2id$probe_id,]
probe2id = probe2id[probe2id$probe_id %in% rownames(data.eset),]
rownames(data.eset) = probe2id$symbol
data.eset = data.eset[rownames(data.eset) %in% vm$V1,]
data.eset = avereps(data.eset)
deg1 <- probe2id %>% 
  inner_join(deg,by="probe_id") %>% 
  select(-probe_id) %>% 
  arrange(desc(logFC)) %>% 
  distinct(symbol,.keep_all=T)
deg1$type <- case_when(deg1$P.Value<0.05&deg1$logFC>0.5~"up",
                       deg1$P.Value<0.05&deg1$logFC< -0.5~"down",
                       T~"stable")

DEG<-deg1 %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 0.5,
                                   ifelse(logFC > 0.5 ,'Up','Down'),'No change'))) 
ggplot(DEG,aes(logFC, -log10(P.Value )))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#51B1B7", "#E1C855" ,"#E07B54"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG[DEG$symbol %in% mygene,],color="black",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =symbol, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave(paste(gse,"volcano.tiff"),width = 6,height = 5)
gseexp=list()
gseexp[["GSE9750"]]=list(data.eset,pd,DEG)
#############?ڶ??????ݼ?#######
gse="GSE39001"
gset <- getGEO(gse,getGPL = F,destdir = ".",)
data.raw <- ReadAffy(celfile.path = paste("./",gse,"_raw_1",sep = ""))
sampleNames(data.raw) <- sub(pattern = "\\.CEL.gz",replacement = "",sampleNames(data.raw))
data.rma <- rma(data.raw)
data.eset <- exprs(data.rma)
colnames(data.eset) = substr(colnames(data.eset),1,9)
pd <- pData(gset[[1]])
pd <- pd %>% 
  select(geo_accession,title)
colnames(pd) <- c("id","sample")
pd <- pd[order(pd$sample),]
data.eset <- data.eset[,pd$id]

boxplot(data.eset,las=3,mian="after-RMA") 
pd$sample[1:43]="CC"
pd$sample[44:55]="Normal"
group_list <- factor(pd$sample,levels = c("Normal","CC"))

names(group_list) <- pd$id
design <- model.matrix(~group_list)

library(limma)
fit <- lmFit(data.eset,design)
fit1 <- eBayes(fit)
options(options=4)
deg <- topTable(fit1,coef=2,n=Inf) 
cc=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "CC"]])
normal=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "Normal"]])
cc - normal  
deg = cbind(rownames(deg),deg)
colnames(deg)[1]="probe_id"

gpl <- gset[[1]]@annotation

library(AnnoProbe)
library(ggplot2)
library(ggrepel)
data.eset=data.eset[order(rownames(data.eset)),]
probe2id <- idmap(gpl)
probe2id=probe2id[order(probe2id$probe_id),]
data.eset = data.eset[rownames(data.eset) %in% probe2id$probe_id,]
rownames(data.eset) = probe2id$symbol
data.eset = data.eset[rownames(data.eset) %in% vm$V1,]
deg1 <- probe2id %>% 
  inner_join(deg,by="probe_id") %>% 
  select(-probe_id) %>% 
  arrange(desc(logFC)) %>% 
  distinct(symbol,.keep_all=T)
deg1$type <- case_when(deg1$P.Value<0.05&deg1$logFC>0.5~"up",
                       deg1$P.Value<0.05&deg1$logFC< -0.5~"down",
                       T~"stable")

DEG<-deg1 %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 0.5,
                                   ifelse(logFC > 0.5 ,'Up','Down'),'No change'))) 
ggplot(DEG,aes(logFC, -log10(P.Value )))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#51B1B7", "#E1C855" ,"#E07B54"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG[DEG$symbol %in% mygene,],color="black",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =symbol, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave(paste("GSE39001_1 volcano.tiff"),width = 6,height = 5)
gseexp[["GSE39001_1"]] = list(data.eset,pd,DEG)

###########################3
gse="GSE39001"
gset <- getGEO(gse,getGPL = F,destdir = ".",)
data.raw <- ReadAffy(celfile.path = paste("./",gse,"_raw_2",sep = ""))
sampleNames(data.raw) <- sub(pattern = "\\.CEL.gz",replacement = "",sampleNames(data.raw))
data.rma <- rma(data.raw)
data.eset <- exprs(data.rma)
colnames(data.eset) = substr(colnames(data.eset),1,9)
pd <- pData(gset[[2]])
pd <- pd %>% 
  select(geo_accession,title)
colnames(pd) <- c("id","sample")
pd <- pd[order(pd$sample),]
data.eset <- data.eset[,pd$id]
pd$sample
pd$sample[1:19]="CC"
pd$sample[20:24]="Normal"
group_list <- factor(pd$sample,levels = c("Normal","CC"))

names(group_list) <- pd$id
design <- model.matrix(~group_list)

library(limma)
fit <- lmFit(data.eset,design)
fit1 <- eBayes(fit)
options(options=4)
deg <- topTable(fit1,coef=2,n=Inf) 
cc=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "CC"]])
normal=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "Normal"]])
cc - normal  
deg = cbind(rownames(deg),deg)
colnames(deg)[1]="probe_id"

gpl <- gset[[2]]@annotation

library(AnnoProbe)
library(ggplot2)
library(ggrepel)
data.eset=data.eset[order(rownames(data.eset)),]
probe2id <- idmap(gpl)
probe2id=probe2id[order(probe2id$probe_id),]
data.eset = data.eset[rownames(data.eset) %in% probe2id$probe_id,]
probe2id = probe2id[probe2id$probe_id %in% rownames(data.eset),]
rownames(data.eset) = probe2id$symbol
data.eset = data.eset[rownames(data.eset) %in% vm$V1,]
deg1 <- probe2id %>% 
  inner_join(deg,by="probe_id") %>% 
  select(-probe_id) %>% 
  arrange(desc(logFC)) %>% 
  distinct(symbol,.keep_all=T)
deg1$type <- case_when(deg1$P.Value<0.05&deg1$logFC>0.5~"up",
                       deg1$P.Value<0.05&deg1$logFC< -0.5~"down",
                       T~"stable")

DEG<-deg1 %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 0.5,
                                   ifelse(logFC > 0.5 ,'Up','Down'),'No change'))) 
ggplot(DEG,aes(logFC, -log10(P.Value )))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#51B1B7", "#E1C855" ,"#E07B54"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG[DEG$symbol %in% mygene,],color="black",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =symbol, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave(paste("GSE39001_2 volcano.tiff"),width = 6,height = 5)
gseexp[["GSE39001_2"]] = list(data.eset,pd,DEG)

###########?????????ݼ?##############33

gse="GSE63514"
gset <- getGEO(gse,getGPL = F,destdir = ".",)
data.raw <- ReadAffy(celfile.path = paste("./",gse,"_raw",sep = ""))
sampleNames(data.raw) <- sub(pattern = "\\.CEL.gz",replacement = "",sampleNames(data.raw))
data.rma <- rma(data.raw)
data.eset <- exprs(data.rma)
colnames(data.eset) = substr(colnames(data.eset),1,10)
pd <- pData(gset[[1]])
pd <- pd %>% 
  select(geo_accession,title)
colnames(pd) <- c("id","sample")
pd <- pd[order(pd$sample),]
data.eset <- data.eset[,pd$id]
pd$sample
pd$sample[1:28]="CC"
pd$sample[105:128]="Normal"
pd$sample[29:42]="CIN1"
pd$sample[43:64]="CIN2"
pd$sample[65:104]="CIN3"
group_list <- factor(pd$sample[pd$sample %in% c("Normal","CC")],levels = c("Normal","CC"))

names(group_list) <- pd$id[pd$sample %in% c("Normal","CC")]
design <- model.matrix(~group_list)
data = data.eset
data.eset=data.eset[,colnames(data.eset) %in% pd$id[pd$sample %in% c("CC","Normal")]]
library(limma)
fit <- lmFit(data.eset,design)
fit1 <- eBayes(fit)
options(options=4)
deg <- topTable(fit1,coef=2,n=Inf) 
cc=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "CC"]])
normal=mean(data.eset[rownames(data.eset) %in% rownames(deg)[1],colnames(data.eset) %in% pd$id[pd$sample %in% "Normal"]])
cc - normal  
deg = cbind(rownames(deg),deg)
colnames(deg)[1]="probe_id"

gpl <- gset[[1]]@annotation

library(AnnoProbe)
library(ggplot2)
library(ggrepel)
data.eset=data.eset[order(rownames(data.eset)),]
probe2id <- idmap(gpl)
probe2id=probe2id[order(probe2id$probe_id),]
data.eset = data.eset[rownames(data.eset) %in% probe2id$probe_id,]
probe2id = probe2id[probe2id$probe_id %in% rownames(data.eset),]
rownames(data.eset) = probe2id$symbol
data.eset = data.eset[rownames(data.eset) %in% vm$V1,]
deg1 <- probe2id %>% 
  inner_join(deg,by="probe_id") %>% 
  select(-probe_id) %>% 
  arrange(desc(logFC)) %>% 
  distinct(symbol,.keep_all=T)
deg1$type <- case_when(deg1$P.Value<0.05&deg1$logFC>0.5~"up",
                       deg1$P.Value<0.05&deg1$logFC< -0.5~"down",
                       T~"stable")

DEG<-deg1 %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 0.5,
                                   ifelse(logFC > 0.5 ,'Up','Down'),'No change'))) 
ggplot(DEG,aes(logFC, -log10(P.Value )))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#51B1B7", "#E1C855" ,"#E07B54"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG[DEG$symbol %in% mygene,],color="black",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =symbol, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave(paste(gse,"volcano.tiff"),width = 6,height = 5)
gseexp[[gse]] = list(data.eset,pd,DEG)

#########TCGA############3
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
load("TCGA-CESC.Rdata")

count <- assay(data,"unstranded")
see  =rowData(data)
rownames(count) = see$gene_name
see = colData(data)
see = see$shortLetterCode
N = count[,see %in% "NT"]
T = count[,see %in% c("TP","TM")]
count = cbind(N,T)
library(limma)
count = avereps(count)
save(count,file="TCGA_count.Rdata")
data.eset = count[rownames(count) %in% vm$V1,]
library(dplyr)
library(ggrepel)

library(ggsci)
list <- c(rep("Normal", 3), rep("CC",306)) %>% factor(., levels = c("CC", "Normal"), ordered = F)
list <- model.matrix(~factor(list)+0) 
colnames(list) = c("CC","Normal")

count<-voom(count,list)#
df.fit <- lmFit(count, list) 
df.matrix <- makeContrasts( CC - Normal , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
DEG <- topTable(fit,n = Inf, adjust = "fdr")
DEG$symbol = rownames(DEG)

DEG<-DEG %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 0.5,
                                   ifelse(logFC > 0.5 ,'Up','Down'),'No change'))) 
diffgene =list(TCGA = DEG[DEG$change %in% c("Down","Up"),])
ggplot(DEG,aes(logFC, -log10(P.Value )))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  scale_color_manual(values =c("#51B1B7", "#E1C855" ,"#E07B54"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = DEG[DEG$symbol %in% mygene,],color="black",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label =symbol, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(p-value)")
ggsave("TCGA_Valcano.tiff",width = 6,height = 5)
gseexp[["TCGA"]]=list(data.eset = count,pd = cbind(colnames(count),c(rep("Normal",3),rep("CC",306))),deg1=DEG)

save(gseexp,file="all_matirx.Rdata")