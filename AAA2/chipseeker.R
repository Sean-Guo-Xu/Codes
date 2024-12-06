setwd("D:\\AAA\\chipseq\\H3K27")
library(ChIPseeker)
covplot("84703_peaks.bed")
peakAnno <- annotatePeak("84703_peaks.bed", 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
diff=read.table("EP300.txt",sep="\t",header = T)
name=c(1:22)
name=c(name,"X","Y")
diff=diff[diff$seqnames %in% name,]
diff$seqnames = paste("chr",diff$seqnames,sep="")
out = cbind(diff$seqnames,diff$start,diff$end,diff$peakid,diff$logFC)
write.table(out,"EP300diff.bed",sep="\t",quote=F,col.names = F,row.names = F)

over=read.table("overlapboth.bed")
over=over[over$V5 > 0,]
over=diff[diff$peakid %in% over$V4,]
over=over[!duplicated(over$GeneID),]
over=over[over$GeneID %in% upre$Name,]
covplot("overlapboth.bed")
