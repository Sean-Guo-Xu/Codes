setwd("E:\\AAA_gwas")
library(CMplot)
allgwas = read.table("E:\\AAA_gwas\\Finngen_R10_I9_ABAORTANEUR",sep="\t")

allgwas = allgwas[order(allgwas$V7),]
allgwas = allgwas[!duplicated(allgwas$V5),]
data = cbind(allgwas$V5,allgwas$V1,allgwas$V2,allgwas$V7)
colnames(data)=c("SNP","CHR","BP","pvalue")
data=as.data.frame(data)

snp=read.table("MR_SNP_res.txt",header = T,sep = "\t")
mygenes <- c("CSK","PLEKHJ1","AMH","LDAH","NEK9","PSRC1","SCAPER","HLA-DRB1","FAM66C")  # 示例基因

snp = snp[snp$Gene %in% mygenes,]
snp = snp[!duplicated(snp$Gene),]
col = rep("#A8D08D",length(mygenes))

CMplot(data,  plot.type="m",
       LOG10=TRUE, threshold=1e-04, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,col = c("#6BB7CA",  "#E07B54"),
       chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,15),chr.border = F,band = 0,highlight.text.font = 18,
       file="png", file.output=TRUE, width=15, height=9, verbose=TRUE,highlight =snp$SNP,highlight.col =col,highlight.text = snp$Gene)


