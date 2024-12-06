library(biomaRt)
library(limma)
library(dplyr)
setwd("D:\\AAA\\chipseq")
counts= read.table("exp.txt",header = T)
counts = counts[!duplicated(counts$Geneid),]
rownames(counts)=counts$Geneid
counts=counts[,-1]
list <-  factor(c(rep("EP300", 3), rep("NC",3)), levels = c("EP300", "NC"), ordered = F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("EP300", "NC")
df.fit <- lmFit(counts, list)  
df.matrix <- makeContrasts(EP300 - NC , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
id = getBM(attributes = c("mgi_symbol", "ensembl_gene_id", "refseq_mrna"), 
           filters = "mgi_symbol",
           values = tempOutput$ID,
           mart = ensembl)
out = merge(id,tempOutput,by.x = "hgnc_symbol",by.y="ID")
out = out[!duplicated(out$refseq_mrna),]
out = out[!(out$refseq_mrna %in% ""),]
colnames(out)[3]="#ID"
write.table(out[,3:9],"exp.txt",sep="\t",quote = F,row.names = F,col.names = T)
