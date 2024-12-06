library(ChIPseeker)
library(ggplot2)
library(ggsci)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
setwd("D:\\AAA\\ATACseq")
namelist=c("AAA1","AAA2","N")
for (name in namelist) {
  peak  = read.table(paste(name,"bed",sep = "."))
  peak$V1 = gsub("chrchr","",peak$V1)
  peak$V1 = gsub("chr","",peak$V1)
peak$V1 = paste("chr",peak$V1,sep = "")
write.table(peak,paste(name,"bed",sep = "."),row.names = F,col.names = F,quote = F,sep = "\t")
}
pal_d3()(10)
color= c(pal_d3()(10),"#6BB7CA")
for (name in namelist) {
  peak <- readPeakFile(paste(name,"bed",sep = "."))
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  peakAnno <- annotatePeak(
    peak,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db")
  pdf(paste(name,"anno.pdf"),width=7,height = 5)
  plotAnnoPie(peakAnno,col =color)
  dev.off()
  see=as.data.frame(peakAnno)
  see=see[grep("Promoter",see$annotation),]
  see= see[order(see$V5,decreasing = T),]
  see=see[!duplicated(see$SYMBOL),]
  data=cbind(see$SYMBOL,see$V5)
  colnames(data) =c("Name","Score")
  
  write.table(data,paste(name,"chipseq_omics.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
  
  
}

EIF2B2