setwd("D:\\AAA\\eqtl")
library(VennDiagram)
asmr=read.table("D:\\AAA\\eqtl\\Artery\\Artery.smr",header = T)
bsmr = read.table("D:\\AAA\\eqtl\\Blood\\Blood.smr",header = T)
asmr = asmr[asmr$p_SMR <0.05,]
asmr = asmr[asmr$p_GWAS<0.0005,]
bsmr = bsmr[bsmr$p_SMR <0.05,]
bsmr = bsmr[bsmr$p_GWAS<0.0005,]
acoloc = read.table("D:\\AAA\\eqtl\\Artery\\Artery_coloc_result.txt",header = T)
bcoloc = read.table("D:\\AAA\\eqtl\\Blood\\Blood_coloc_result.txt",header = T)
phy = c("Blood_SMR","Artery_SMR","Blood_coloc","Artery_coloc")
col = c("#E64B35FF" ,"#4DBBD5FF" ,"#00A087FF","#3C5488FF")
phylist  = list(bsmr$Symbol,asmr$Symbol,bcoloc$Symbol,acoloc$Symbol)
venn= venn.diagram(phylist,compression = "lzw",filename = "venn.tiff",imagetype="tiff",fill=col,height = 2500,width = 3000, cat.default.pos = "outer",category.names = phy)


load("E:\\AAA_scRNA\\blood_diff_all.Rdata")
bdiff = diff
load("E:\\AAA_scRNA\\artery_diff_all.Rdata")
adiff=diff
library(pheatmap)
load("E:\\AAA_scRNA\\")
adiff=