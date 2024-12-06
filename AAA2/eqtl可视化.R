library(VariantAnnotation)
library(gwasglue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Hs.eg.db)
library(CMplot)

allgwas = read.table("D:\\AAA\\MR\\Finngen_R9_I9_ABAORTANEUR",sep="\t")
allgwas = allgwas[order(allgwas$V7),]
allgwas = allgwas[!duplicated(allgwas$V5),]



data = cbind(allgwas$V5,allgwas$V1,allgwas$V2,allgwas$V7)
colnames(data)=c("SNP","CHR","BP","pvalue")
data=as.data.frame(data)
CMplot(data,  plot.type="m",
       LOG10=TRUE, threshold=1e-04, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,col = c("#6BB7CA",  "#E07B54"),
       chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,15),chr.border = F,band = 0,highlight.text.font = 18,
       file="pdf", file.output=TRUE, width=15, height=9, verbose=TRUE,highlight = c("rs4937515","rs12882111","rs10170771"),highlight.col = c("#A8D08D","#A8D08D","#A8D08D"),highlight.text = c("rs4937515","rs12882111","rs10170771"))
library(eQTpLot)
library(BiocFileCache)
library(biomaRt)
library(ggsci)
gwas_df$PHE =NULL
setwd("D:\\AAA\\eqtl\\smr")

snplist =  c("rs4937515","rs12882111","rs10170771")
allgene = c("ADAMTS8","NEK9","LDAH")
snp="rs10170771"
gene="LDAH"
tissue = "Artery"
gwas_df = cbind(allgwas$V1,allgwas$V2,allgwas$V5,allgwas$V7,allgwas$V9,rep("CC",nrow(allgwas)))
colnames(gwas_df) = c("CHR","BP","SNP","P","BETA","PHE") 
gwas_df = as.data.frame(gwas_df)
gwas_df$CHR = as.integer(gwas_df$CHR)
gwas_df$BP = as.integer(gwas_df$BP)
gwas_df$P =as.numeric(gwas_df$P)
gwas_df$BETA = as.numeric(gwas_df$BETA)
cen = gwas_df[gwas_df$SNP %in% snp,]
commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery_Aorta --snp",snp," --query 1 --snp-wind 200 --out eqtl_print")
system(commond,intern = T)
eqtl = read.table("eqtl_print.txt",header=T)
eqtl_df = cbind(eqtl$Probe,eqtl$SNP,eqtl$p,eqtl$b,rep(tissue,nrow(eqtl)))
eqtl_df = as.data.frame(eqtl_df)
colnames(eqtl_df) = c("Gene.Symbol","SNP.Id","P.Value","NES","Tissue")
ensembls <- mapIds(org.Hs.eg.db, keys = eqtl_df$Gene.Symbol, keytype = "ENSEMBL", column="SYMBOL")
eqtl_df$Gene.Symbol= ensembls
eqtl_df$P.Value =as.numeric(eqtl_df$P.Value)
eqtl_df$NES =as.numeric(eqtl_df$NES)
gwas_df = gwas_df[gwas_df$CHR %in% cen[,1],]
gwas_df =gwas_df[gwas_df$SNP %in% eqtl_df$SNP.Id,]
eqtl =eqtl[eqtl$SNP %in% gwas_df$SNP,]
eqtl=eqtl[!duplicated(eqtl$SNP),]
gwas_df = gwas_df[order(gwas_df$CHR),]
eqtl = eqtl[order(eqtl$SNP),]
gwas_df = gwas_df[order(gwas_df$SNP),]
gwas_df$PHE = rep("AAA",nrow(gwas_df))
eqtl_df = eqtl_df[order(eqtl_df$P.Value),]
all = merge(gwas_df,eqtl_df[eqtl_df$SNP.Id %in% gwas_df$SNP,],by.x="SNP",by.y="SNP.Id")
all = all[all$P<1e-04,]
all = all[all$P.Value<1e-08,]
write.table(all$SNP,paste("D:\\AAA\\qtl_scRNA\\artery",gene,"rs.txt",sep = "_"),row.names = F,col.names = F,quote = F)
source("E:\\cervix\\eqtlplot.R")
ploteqtl(GWAS.df = gwas_df, eQTL.df = eqtl_df, gene =  gene, 
         gbuild = "hg38",  trait = "AAA", tissue = tissue, CollapseMethod = "min",sigpvalue_GWAS = 1e-04,sigpvalue_eQTL = 0.0001,range = 200)

ggsave(paste(gene,tissue,"eqtplot.tiff"),width = 15,height = 18)

snp="rs10170771"
gene="LDAH"
tissue = "Blood"
gwas_df = cbind(allgwas$V1,allgwas$V2,allgwas$V5,allgwas$V7,allgwas$V9,rep("CC",nrow(allgwas)))
colnames(gwas_df) = c("CHR","BP","SNP","P","BETA","PHE") 
gwas_df = as.data.frame(gwas_df)
gwas_df$CHR = as.integer(gwas_df$CHR)
gwas_df$BP = as.integer(gwas_df$BP)
gwas_df$P =as.numeric(gwas_df$P)
gwas_df$BETA = as.numeric(gwas_df$BETA)
cen = gwas_df[gwas_df$SNP %in% snp,]
commond = paste("smr-1.3.1-win.exe --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --snp",snp," --query 1 --snp-wind 200 --out eqtl_print")
system(commond,intern = T)
eqtl = read.table("eqtl_print.txt",header=T)
eqtl_df = cbind(eqtl$Gene,eqtl$SNP,eqtl$p,eqtl$b,rep(tissue,nrow(eqtl)))
eqtl_df = as.data.frame(eqtl_df)
colnames(eqtl_df) = c("Gene.Symbol","SNP.Id","P.Value","NES","Tissue")
ensembls <- mapIds(org.Hs.eg.db, keys = eqtl_df$Gene.Symbol, keytype = "ENSEMBL", column="SYMBOL")
eqtl_df$Gene.Symbol= ensembls
eqtl_df$P.Value =as.numeric(eqtl_df$P.Value)
eqtl_df$NES =as.numeric(eqtl_df$NES)
gwas_df = gwas_df[gwas_df$CHR %in% cen[,1],]
gwas_df =gwas_df[gwas_df$SNP %in% eqtl_df$SNP.Id,]
eqtl =eqtl[eqtl$SNP %in% gwas_df$SNP,]
eqtl=eqtl[!duplicated(eqtl$SNP),]
gwas_df = gwas_df[order(gwas_df$CHR),]
eqtl = eqtl[order(eqtl$SNP),]
gwas_df = gwas_df[order(gwas_df$SNP),]
gwas_df$PHE = rep("AAA",nrow(gwas_df))
eqtl_df = eqtl_df[order(eqtl_df$P.Value),]
all = merge(gwas_df,eqtl_df[eqtl_df$SNP.Id %in% gwas_df$SNP,],by.x="SNP",by.y="SNP.Id")
all = all[all$P<1e-04,]
all = all[all$P.Value<1e-08,]
write.table(all$SNP,paste("D:\\AAA\\qtl_scRNA\\blood",gene,"rs.txt",sep = "_"),row.names = F,col.names = F,quote = F)
source("E:\\cervix\\eqtlplot.R")
ploteqtl(GWAS.df = gwas_df, eQTL.df = eqtl_df, gene =  gene, 
         gbuild = "hg38",  trait = "AAA", tissue = tissue, CollapseMethod = "min",sigpvalue_GWAS = 1e-04,sigpvalue_eQTL = 0.0001,range = 200)
ggsave(paste(gene,tissue,"eqtplot.tiff"),width = 15,height = 18)
