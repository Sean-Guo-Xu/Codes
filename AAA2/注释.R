library(ChIPseeker)
library(ggplot2)
require("biomaRt")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
setwd("D:\\AAA\\chipseq\\H3K27")
namelist=c("EP300AAA","EP300Normal","104699_peaks","104697_peaks")
for (name in namelist) {
peak <- readPeakFile(paste(name,"bed",sep = "."))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(
  peak,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db")
pdf(paste(name,"anno.pdf"),width=7,height = 5)
plotAnnoPie(peakAnno)
dev.off()
see=as.data.frame(peakAnno)
see=see[grep("Promoter",see$annotation),]
see= see[order(see$V5,decreasing = T),]
see=see[!duplicated(see$SYMBOL),]
data=cbind(see$SYMBOL,see$V5)
colnames(data) =c("Name","Score")

write.table(data,paste(name,"chipseq_omics.txt"),sep = "\t",quote = F,row.names = F,col.names = T)


}

##########ส๓ิด
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")

namelist=c("88612_peaks","84703_peaks")
for (name in namelist) {
  peak <- readPeakFile(paste(name,"bed",sep = "."))
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  peakAnno <- annotatePeak(
    peak,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db")
  
  pdf(paste(name,"anno.pdf"),width=7,height = 5)
  plotAnnoPie(peakAnno)
  dev.off()
  see=as.data.frame(peakAnno)
  see=see[grep("Promoter",see$annotation),]
  see= see[order(see$V5,decreasing = T),]
  see=see[!duplicated(see$SYMBOL),]
  data=cbind(see$SYMBOL,see$V5)
  data=as.data.frame(data)
  colnames(data) =c("Name","Score")
  genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = data[,1], 
                 mart = mouse, 
                 attributesL = c("hgnc_symbol"), 
                 martL = human, uniqueRows=T)
  data=merge(data,genes,by.x="Name",by.y="MGI.symbol")
  data=data[,-1]
  data=cbind(data$HGNC.symbol,data$Score)
  colnames(data) =c("Name","Score")
  write.table(data,paste(name,"chipseq_omics.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
}

