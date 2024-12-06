library("coloc")
library(TwoSampleMR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(locuscomparer)
library(ggplot2)

setwd("E:\\SMR")
gwas = gwas[,c(4,6,5,8,11,12,14,10)]
colnames(gwas) = c("SNP","A1","A2","freq","b","se","p","n")
gwas = gwas[gwas$SNP != "", ]
gwas = gwas[gwas$A1 %in% c("A","T","G","C"),]
gwas = gwas[gwas$A2 %in% c("A","T","G","C"),]
gwas$SNP = gsub(",.+","",gwas$SNP) ##delete multiple SNPs
write.table(gwas,"Asian_AAA_gwas.txt",sep = "\t",row.names = F,col.names = T,quote = F)

load("E:\\AAA_gwas\\asian_gwas.Rdata")
system("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary Asian_AAA_gwas.txt --beqtl-summary Artery/Artery_Aorta --out Asian_Artery --thread-num 10 --diff-freq-prop 0.99")
system("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary Asian_AAA_gwas.txt --beqtl-summary Blood/Blood --out Asian_Blood --thread-num 10 --diff-freq-prop 0.99")
vector = c(1:22)
for (i in vector) {
  command = paste("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary Asian_AAA_gwas.txt --beqtl-summary mqtl/bl_mqtl_chr",i," --out Asian_mqtl_",i," --thread-num 10 --diff-freq-prop 0.99",sep = "")
  system(command)
}
allsmr = NULL                          
for (i in vector) {
  smr=read.table(paste("Asian_mqtl_",i,".smr",sep = ""),sep = "\t",header = T)
  allsmr = rbind(allsmr,smr)
}
allsmr = allsmr[allsmr$p_SMR<0.05 & allsmr$p_GWAS<0.0001,]

smrgene = allgwas[allgwas$V5 %in% allsmr$topSNP,5:6]
for (i in unique(allsmr$topSNP)) {
  allsmr[allsmr$topSNP %in% i, 3] = smrgene[smrgene$V5 %in% i,2]
}
write.table(allsmr,"E:\\AAA_gwas\\Asian_mqtl_smr.txt",,row.names = F,col.names = T,sep = "\t",quote = F)

allsmr =read.table("Asian_Blood.smr",sep = "\t",header = T)
allsmr = allsmr[allsmr$p_SMR<0.05 & allsmr$p_GWAS<0.0001,]
#smrgene = gwas[gwas$SNPID %in% allsmr$topSNP,5:6]
#for (i in unique(allsmr$topSNP)) {
  allsmr[allsmr$topSNP %in% i, 3] = smrgene[smrgene$V5 %in% i,2]
}
write.table(allsmr,"E:\\AAA_gwas\\Asian_Blood_smr.txt",,row.names = F,col.names = T,sep = "\t",quote = F)

allsmr =read.table("Asian_Artery.smr",sep = "\t",header = T)
allsmr = allsmr[allsmr$p_SMR<0.05 & allsmr$p_GWAS<0.0001,]
#smrgene = allgwas[allgwas$V5 %in% allsmr$topSNP,5:6]
#for (i in unique(allsmr$topSNP)) {
#  allsmr[allsmr$topSNP %in% i, 3] = smrgene[smrgene$V5 %in% i,2]
#}
write.table(allsmr,"E:\\AAA_gwas\\Asian_Artery_smr.txt",,row.names = F,col.names = T,sep = "\t",quote = F)

setwd("E:\\AAA_gwas")
as=read.table("Asian_Artery_smr.txt",header=T)
bs=read.table("Asian_Blood_smr.txt",header=T)
mq=read.table("Asian_mqtl_smr.txt",header = T)
allsnp = c(as$topSNP,bs$topSNP,mq$topSNP)
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas =  174756   #gwas样本量
case = 1155#得病数
neqtl = 31684 #eqtl样本量

allgwas = gwas[order(gwas$p.value),]
allgwas = allgwas[!duplicated(allgwas$SNPID),]
mygwas = allgwas[allgwas$p.value < 0.0001,]
result = list()
gwaslist =list()
setwd("E:\\SMR")

for (i in unique(allsnp)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Blood/Blood  --snp",i," --query 1 --snp-wind 100 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl = eqtl[order(eqtl$p),]
  eqtl = eqtl[!duplicated(eqtl$SNP),]
  gwas = mygwas[mygwas$SNPID %in% eqtl$SNP,]
  colnames(gwas)[4] = "SNP"
  input = merge(gwas,eqtl,by="SNP")
  input$p[input$p %in% 0] = 1e-300 
  if(nrow(input) ==0 ){
    next
  }
  re= coloc.abf(dataset1 = list(pvalues=input$p.value,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                dataset2 = list(pvalues=input$p,snp=input$SNP,type="quant",N=neqtl), MAF = input$Freq)
  if(re$summary[6]>0.75){
    
    input$H4 = re$summary[6]
    result[[i]] = input
  }
}

save(result,file="E:\\AAA_gwas\\Asian_Blood_coloc.Rdata")
#save(gwaslist,file="E:\\AAA_gwas\\Blood_gene.Rdata")
################GTEX##############
setwd("E:\\AAA_gwas")
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas =  174756   #gwas样本量
case = 1155#得病数
neqtl = 387 #eqtl样本量


alleqtl = read.table("E:\\AAA_gwas\\Artery_eqtl.txt",header=T)
colnames(alleqtl)[19] = "SNP"
alleqtl = alleqtl[order(alleqtl$pval_nominal),]
alleqtl = alleqtl[!duplicated(alleqtl$SNP),]

result = list()
gwaslist =list()
setwd("E:\\SMR")
allsnp = unique(c(as$topSNP,bs$topSNP,mq$topSNP))
for (i in allsnp) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery/Artery_Aorta --snp",i," --query 1 --snp-wind 100 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl = eqtl[order(eqtl$p),]
  eqtl = eqtl[!duplicated(eqtl$SNP),]
  eqtl = eqtl[eqtl$SNP %in% alleqtl$SNP,]
  eqtl = eqtl[order(eqtl$SNP),]
  partalleqtl = alleqtl[alleqtl$SNP %in% eqtl$SNP, ]
  eqtl$Freq = partalleqtl$maf
  gwas = mygwas[mygwas$SNPID %in% eqtl$SNP,]
  colnames(gwas)[4] = "SNP"
  input = merge(gwas,eqtl,by="SNP")
  input$p[input$p %in% 0] = 1e-300 
  if(nrow(input) ==0 ){
    next
  }
  re= coloc.abf(dataset1 = list(pvalues=input$p.value,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                dataset2 = list(pvalues=input$p,snp=input$SNP,type="quant",N=neqtl), MAF = input$Freq)
  if(re$summary[6]>0.75){
    
    input$H4 = re$summary[6]
    result[[i]] = input
  }
}
save(result,file="E:\\AAA_gwas\\Asian_coloc_result.Rdata")

##############基因名注释##################333

library(clusterProfiler)
library(org.Hs.eg.db)
load("E:\\AAA_gwas\\Asian_coloc_result.Rdata")
setwd("E:\\SMR")
egene = NULL
for (i in names(result)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery/Artery_Aorta --snp",i," --snp-wind 100, --query 0.0001 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl  =eqtl[,c(1,14,12,13,7)]
  ensembls <- mapIds(org.Hs.eg.db, keys = eqtl$Probe, keytype = "ENSEMBL", column="SYMBOL")
  eqtl$Probe= ensembls
  colnames(eqtl)[5] = "Gene"
  gwas = result[[i]]
  gwas = gwas[gwas$SNP %in% eqtl$SNP,]
  gwas = gwas[gwas$p.value < 0.0001,]
  gwas = gwas[,c(1,14,11,12,34)]
  colnames(gwas) = c("SNP","p","b","SE","H4")
  output= merge(gwas,eqtl,by="SNP")
  
  if(nrow(eqtl)<= 0){
    next
  }
  output = output[order(output$p.y),]
  output = output[!duplicated(output$Gene),]
  egene = rbind(egene,output)
}

egene = egene[!duplicated(egene$Gene) | !duplicated(egene$SNP),]

write.table(egene,"E:\\AAA_gwas\\Asian_Artery_coloc.txt",col.names = T,row.names = F,quote = F,sep="\t")

###############################


load("E:\\AAA_gwas\\Asian_Blood_coloc.Rdata")
setwd("E:\\SMR")
egene = NULL
for (i in names(result)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Blood/Blood --snp",i," --snp-wind 100, --query 0.0001 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  ensembls <- mapIds(org.Hs.eg.db, keys = eqtl$Probe, keytype = "ENSEMBL", column="SYMBOL")
  eqtl$Gene= ensembls
  eqtl  =eqtl[,c(1,14,12,13,10)]
  gwas = result[[i]]
  gwas = gwas[gwas$SNP %in% eqtl$SNP,]
  gwas = gwas[gwas$p.value < 0.0001,]
  gwas = gwas[,c(1,14,11,12,34)]
  colnames(gwas) = c("SNP","p","b","SE","H4")
  output= merge(gwas,eqtl,by="SNP")
  
  if(nrow(eqtl)<= 0){
    next
  }
  output = output[order(output$p.y),]
  output = output[!duplicated(output$Gene),]
  egene = rbind(egene,output)
  
}
egene = egene[!duplicated(egene$Gene) | !duplicated(egene$SNP),]
write.table(egene,"E:\\AAA_gwas\\Asian_Blood_coloc.txt",col.names = T,row.names = F,quote = F,sep="\t")

