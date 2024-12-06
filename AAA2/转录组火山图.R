library(affy)
library(GEOquery)
library(affyPLM)
library(dplyr)
setwd("D:\\AAA\\qtl_scRNA")
mygene="HDAC10"
gse="GSE57691"
gset <- getGEO(gse,getGPL = F,destdir = ".")
data.raw <- ReadAffy(celfile.path = "./GSE57691_series_matrix.txt")
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