library("coloc")
library(TwoSampleMR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

setwd("E:\\SMR")

system("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary Europe_AAA_gwas.txt --beqtl-summary Artery/Artery_Aorta --out Europe_Artery --thread-num 10 --diff-freq-prop 0.99")
system("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary Europe_AAA_gwas.txt --beqtl-summary Blood/Blood --out Europe_Blood --thread-num 10 --diff-freq-prop 0.99")
vector = c(1:22)
for (i in vector) {
  command = paste("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary Europe_AAA_gwas.txt --beqtl-summary mqtl/bl_mqtl_chr",i," --out Europe_mqtl_",i," --thread-num 10 --diff-freq-prop 0.99",sep = "")
  system(command)
}
allsmr = NULL
for (i in vector) {
  smr=read.table(paste("Europe_mqtl_",i,".smr",sep = ""),sep = "\t",header = T)
  allsmr = rbind(allsmr,smr)
}
allsmr = allsmr[allsmr$p_SMR<0.05 & allsmr$p_GWAS<0.0001,]
allgwas = read.table("E:\\AAA_gwas\\finngen_R10_I9_ABAORTANEUR",sep="\t")
smrgene = allgwas[allgwas$V5 %in% allsmr$topSNP,5:6]
for (i in unique(allsmr$topSNP)) {
  allsmr[allsmr$topSNP %in% i, 3] = smrgene[smrgene$V5 %in% i,2]
}
write.table(allsmr,"E:\\AAA_gwas\\Europe_mqtl_smr.txt",,row.names = F,col.names = T,sep = "\t",quote = F)

allsmr =read.table("Europe_Blood.smr",sep = "\t",header = T)
allsmr = allsmr[allsmr$p_SMR<0.05 & allsmr$p_GWAS<0.0001,]
smrgene = allgwas[allgwas$V5 %in% allsmr$topSNP,5:6]
for (i in unique(allsmr$topSNP)) {
  allsmr[allsmr$topSNP %in% i, 3] = smrgene[smrgene$V5 %in% i,2]
}
write.table(allsmr,"E:\\AAA_gwas\\Europe_Blood_smr.txt",,row.names = F,col.names = T,sep = "\t",quote = F)

allsmr =read.table("Europe_Artery.smr",sep = "\t",header = T)
allsmr = allsmr[allsmr$p_SMR<0.05 & allsmr$p_GWAS<0.0001,]
smrgene = allgwas[allgwas$V5 %in% allsmr$topSNP,5:6]
for (i in unique(allsmr$topSNP)) {
  allsmr[allsmr$topSNP %in% i, 3] = smrgene[smrgene$V5 %in% i,2]
}
write.table(allsmr,"E:\\AAA_gwas\\Europe_Artery_smr.txt",,row.names = F,col.names = T,sep = "\t",quote = F)

setwd("E:\\AAA_gwas")
load("Europe_gwas.Rdata")
as=read.table("Europe_Artery_smr.txt",header=T)
bs=read.table("Europe_Blood_smr.txt",header=T)
mq=read.table("Europe_mqtl_smr.txt",header = T)
allsnp = c(as$topSNP,bs$topSNP,mq$topSNP)
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas =429209  #gwas样本量
neqtl = 31684 #eqtl样本量
case = 3869#得病数
allgwas=gwas
allgwas = allgwas[order(allgwas$V7),]
allgwas = allgwas[!duplicated(allgwas$V5),]
mygwas = allgwas[allgwas$V7 < 0.0001,]
result = list()
gwaslist =list()
setwd("E:\\SMR")

for (i in unique(allsnp)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Blood/Blood  --snp",i," --query 1 --snp-wind 100 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl = eqtl[order(eqtl$p),]
  eqtl = eqtl[!duplicated(eqtl$SNP),]
  gwas = mygwas[mygwas$V5 %in% eqtl$SNP,]
  colnames(gwas)[5] = "SNP"
  input = merge(gwas,eqtl,by="SNP")
  input$p[input$p %in% 0] = 1e-300 
  if(nrow(input) ==0 ){
    next
  }
  re= coloc.abf(dataset1 = list(pvalues=input$V7,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                dataset2 = list(pvalues=input$p,snp=input$SNP,type="quant",N=neqtl), MAF = input$Freq)
  if(re$summary[6]>0.75){
    
    input$H4 = re$summary[6]
    result[[i]] = input
    gwas = allgwas[allgwas$V5 %in% eqtl$SNP,]
    gwas = cbind(gwas$V5,gwas$V7)
    colnames(gwas) = c("rsid","pval")
    eqtl = cbind(eqtl$SNP,eqtl$p)
    colnames(eqtl) = c("rsid","pval")
    eqtl = as.data.frame(eqtl)
    eqtl$pval[eqtl$pval %in% 0] = "1e-300"
    write.table(gwas,"gwas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    write.table(eqtl,"eqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    gwaslist[[i]] = unique(input$V6)
  #  if(!(i %in% eqtl$rsid) & !(i %in% eqtl$pval)){ input = input[order(input$V7),]
  #  locuscompare("eqtl.tsv","gwas.tsv",legend =F,snp = input$SNP[1]) } else {
  #    locuscompare("eqtl.tsv","gwas.tsv",legend =F,snp = i) } 
   # ggsave(paste(i,"locus.png"),width=10,height = 5)
  }
}

  save(result,file="E:\\AAA_gwas\\Europe_Blood_coloc.Rdata")
#save(gwaslist,file="E:\\AAA_gwas\\Blood_gene.Rdata")
################GTEX##############
setwd("E:\\AAA_gwas")
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas = 429209  #gwas样本量
neqtl = 387 #eqtl样本量
case = 3869 #得病数
allgwas = read.table("finngen_R10_I9_ABAORTANEUR",sep="\t")
allgwas = allgwas[order(allgwas$V7),]
allgwas = allgwas[!duplicated(allgwas$V5),]
mygwas = allgwas[allgwas$V7 < 0.0001,]

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
  gwas = mygwas[mygwas$V5 %in% eqtl$SNP,]
  colnames(gwas)[5] = "SNP"
  input = merge(gwas,eqtl,by="SNP")
  input$p[input$p %in% 0] = 1e-300 
  if(nrow(input) ==0 ){
    next
  }
  re= coloc.abf(dataset1 = list(pvalues=input$V7,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                dataset2 = list(pvalues=input$p,snp=input$SNP,type="quant",N=neqtl), MAF = input$Freq)
  if(re$summary[6]>0.5){
    
    input$H4 = re$summary[6]
    result[[i]] = input
    gwas = allgwas[allgwas$V5 %in% eqtl$SNP,]
    gwas = cbind(gwas$V5,gwas$V7)
    colnames(gwas) = c("rsid","pval")
    eqtl = cbind(eqtl$SNP,eqtl$p)
    colnames(eqtl) = c("rsid","pval")
    eqtl = as.data.frame(eqtl)
    eqtl$pval[eqtl$pval %in% 0] = "1e-300"
    write.table(gwas,"gwas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    write.table(eqtl,"eqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    gwaslist[[i]] = unique(input$V6)
}
save(result,file="E:\\AAA_gwas\\Artery_coloc_result.Rdata")

##############基因名注释##################333
load("E:\\AAA_gwas\\Artery_coloc_result.Rdata")
setwd("D:\\AAA\\eqtl\\smr")
library(clusterProfiler)
library(org.Hs.eg.db)
#load("D:\\AAA\\eqtl\\blood\\Blood_coloc_result.Rdata")
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
  gwas = gwas[gwas$V7 < 0.0001,]
  gwas = gwas[,c(1,7,9,10,27)]
  colnames(gwas) = c("SNP","p","b","SE","H4")
  output= merge(gwas,eqtl,by="SNP")
  
  if(nrow(eqtl)<= 0){
    next
  }
  output = output[order(output$p.y),]
  output = output[!duplicated(output$Gene),]
  egene = rbind(egene,output)
}
egene = egene[egene$Gene %in% c(as$gene,bs$gene,mq$Gene),]
egene = egene[!duplicated(egene$Gene) | !duplicated(egene$SNP),]

write.table(egene,"E:\\AAA_gwas\\Europe_Artery_coloc.txt",col.names = T,row.names = F,quote = F,sep="\t")

###############################


load("E:\\AAA_gwas\\Europe_Blood_coloc.Rdata")
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
  gwas = gwas[gwas$SNP %in%  eqtl$SNP,]
  gwas = gwas[gwas$V7 < 0.0001,]
  gwas = gwas[,c(1,7,9,10,27)]
  colnames(gwas) = c("SNP","p","b","SE","H4")
  output= merge(gwas,eqtl,by="SNP")
  
  if(nrow(eqtl)<= 0){
    next
  }
  output = output[order(output$p.y),]
  output = output[!duplicated(output$Gene),]
  egene = rbind(egene,output)
  
}
egene = egene[egene$Gene %in% c(as$gene,bs$gene,mq$Gene),]
egene = egene[!duplicated(egene$Gene) | !duplicated(egene$SNP),]
write.table(egene,"E:\\AAA_gwas\\Europe_Blood_coloc.txt",col.names = T,row.names = F,quote = F,sep="\t")

