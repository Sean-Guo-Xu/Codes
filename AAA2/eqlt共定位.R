library("coloc")
library(TwoSampleMR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(locuscomparer)
library(ggplot2)
setwd("D:\\AAA\\eqtl\\blood")
myrs=read.table("Blood_coloc_result.txt",header = T)
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas = 392423  #gwas样本量
neqtl = 31684 #eqtl样本量
case = 3548 #得病数
allgwas = read.table("D:\\AAA\\MR\\Finngen_R9_I9_ABAORTANEUR",sep="\t")
allgwas = allgwas[order(allgwas$V7),]
allgwas = allgwas[!duplicated(allgwas$V5),]
gwas = allgwas[allgwas$V7 < 0.0001,]
gwas = format_data(
  gwas,
  pos = "V2",
  chr = "V1",
  type='outcome',
  gene = "V6",
  snp_col = "V5",
  beta_col = "V9",
  se_col = "V10",
  effect_allele_col ="V4",
  other_allele_col = "V3",
  eaf_col = "V11",
  pval_col = "V7"
)

gwas_ld <- clump_data(gwas, clump_kb=10000, clump_r2=0.2)
result = list()
gwaslist =list()
setwd("D:\\AAA\\eqtl\\smr")


for (i in gwas_ld$SNP) {
  snp = gwas_ld[gwas_ld$SNP %in% i,]
  commond = paste("smr-1.3.1-win.exe --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --snp",i," --query 1 --snp-wind 100 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  if(i %in% eqtl$SNP){
  eqtl = eqtl[order(eqtl$p),]
  eqtl = eqtl[!duplicated(eqtl$SNP),]
  gwas = allgwas[allgwas$V1 %in% snp$chr.outcom,]
  gwas = gwas[gwas$V5 %in% eqtl$SNP,]
  colnames(gwas)[5] = "SNP"
  input = merge(gwas,eqtl,by="SNP")
  input$p[input$p %in% 0] = 1e-300 
  re= coloc.abf(dataset1 = list(pvalues=input$V7,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                dataset2 = list(pvalues=input$p,snp=input$SNP,type="quant",N=neqtl), MAF = input$Freq)
  if(re$summary[6]>0.5){
    input$H4 = re$summary[6]
    result[[i]] = input
    gwas = cbind(input$SNP,input$V7)
    colnames(gwas) = c("rsid","pval")
    eqtl = cbind(input$SNP,input$p)
    colnames(eqtl) = c("rsid","pval")
    write.table(gwas,"gwas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    write.table(eqtl,"eqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
    gwaslist[[i]] = unique(input$V6)
    locuscompare("eqtl.tsv","gwas.tsv",legend =F,snp = i)  
    ggsave(paste(i,"locus.png"),width=10,height = 5)
    }
}
}

save(result,file="D:\\AAA\\eqtl\\Blood_coloc_result.Rdata")
save(gwaslist,file="D:\\AAA\\eqtl\\Blood_coloc_gene.Rdata")
load("D:\\AAA\\eqtl\\blood\\Blood_coloc_result.Rdata")
setwd("D:\\AAA\\eqtl\\smr")
egene = NULL
for (i in names(result)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --snp",i," --query 0.0001 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl  =eqtl[,c(1,14,12,13,10)]
  gwas = result[[i]]
  gwas = gwas[gwas$SNP %in% i,]
  gwas = gwas[,c(1,7,9,10,27)]
  colnames(gwas) = c("SNP","p","b","SE","H4")
  output= merge(gwas,eqtl,by="SNP")
  if(nrow(eqtl)<= 0){
      next
    }
  egene = rbind(egene,output)
}

egene$Gene = mapIds(x = org.Hs.eg.db,#注释包
       keys =egene$Gene, #需要转换的基因ID
       keytype = "ENSEMBL", #需要转换的类型
       column = "SYMBOL")
write.table(egene,"D:\\AAA\\eqtl\\blood\\coloc_result.txt",col.names = T,row.names = F,quote = F,sep="\t")

################GTEX##############
setwd("D:\\AAA\\eqtl")
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas = 392423  #gwas样本量
neqtl = 387 #eqtl样本量
case = 3548 #得病数
allgwas = read.table("D:\\AAA\\MR\\Finngen_R9_I9_ABAORTANEUR",sep="\t")
allgwas = allgwas[order(allgwas$V7),]
allgwas = allgwas[!duplicated(allgwas$V5),]
gwas = allgwas[allgwas$V7 < 0.0001,]
gwas = format_data(
  gwas,
  pos = "V2",
  chr = "V1",
  type='outcome',
  gene = "V6",
  snp_col = "V5",
  beta_col = "V9",
  se_col = "V10",
  effect_allele_col ="V4",
  other_allele_col = "V3",
  eaf_col = "V11",
  pval_col = "V7"
)
alleqtl = read.table("D:\\AAA\\eqtl\\Artery_eqtl.txt",header=T)
colnames(alleqtl)[19] = "SNP"
alleqtl = alleqtl[order(alleqtl$pval_nominal),]
alleqtl = alleqtl[!duplicated(alleqtl$SNP),]
gwas_ld <- clump_data(gwas, clump_kb=10000, clump_r2=0.2)
result = list()
gwaslist =list()
setwd("D:\\AAA\\eqtl\\smr-1.3.1-win-x86_64")

for (i in gwas_ld$SNP[99:401]) {
  snp = gwas_ld[gwas_ld$SNP %in% i,]
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery_Aorta --snp",i," --query 1 --snp-wind 100 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  if(i %in% eqtl$SNP){
    eqtl = eqtl[order(eqtl$p),]
    eqtl = eqtl[!duplicated(eqtl$SNP),]
    smreqtl = eqtl 
    eqtl  = alleqtl[alleqtl$SNP %in% eqtl$SNP,]
    if(nrow(eqtl)!=0){
    gwas = allgwas[allgwas$V1 %in% snp$chr.outcom,]
    smrgwas = gwas
    gwas = gwas[gwas$V5 %in% eqtl$SNP,]
    if(nrow(gwas) != 0){
    colnames(gwas)[5] = "SNP"
    input = merge(gwas,eqtl,by="SNP")
    input$pval_nominal[input$pval_nominal %in% 0] = 1e-300 
    re= coloc.abf(dataset1 = list(pvalues=input$V7,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                  dataset2 = list(pvalues=input$pval_nominal,snp=input$SNP,type="quant",N=neqtl), MAF = input$maf)
    if(re$summary[6]>0.5){
      input$H4 = re$summary[6]
      result[[i]] = input
      gwas = cbind(smrgwas$V5,smrgwas$V7)
      colnames(gwas) = c("rsid","pval")
      eqtl = cbind(smreqtl$SNP,smreqtl$p)
      colnames(eqtl) = c("rsid","pval")
      eqtl = eqtl[eqtl[,1] %in% gwas[,1],]
      gwas= gwas[gwas[,1] %in% eqtl[,1],]
      write.table(gwas,"gwas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
      write.table(eqtl,"eqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
      gwaslist[[i]] = unique(input$V6)
      locuscompare("eqtl.tsv","gwas.tsv",legend =F,snp = i)  
      ggsave(paste(i,"locus.png"),width=10,height = 5)
    }}
  }}
}

save(result,file="D:\\AAA\\eqtl\\Artery_coloc_result.Rdata")
save(gwaslist,file="D:\\AAA\\eqtl\\Artery_coloc_gene.Rdata")

##############基因名注释##################333
load("D:\\AAA\\eqtl\\Artery\\Artery_coloc_gene.Rdata")
load("D:\\AAA\\eqtl\\Artery\\Artery_coloc_result.Rdata")
setwd("D:\\AAA\\eqtl\\smr")
library(clusterProfiler)
library(org.Hs.eg.db)
#load("D:\\AAA\\eqtl\\blood\\Blood_coloc_result.Rdata")
setwd("D:\\AAA\\eqtl\\smr")
egene = NULL
for (i in names(result)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery_Aorta --snp",i," --snp-wind 100, --query 0.0001 --out eqtl_print")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl  =eqtl[,c(1,14,12,13,10)]
  gwas = result[[i]]
  gwas = gwas[gwas$SNP %in% eqtl$SNP,]
  gwas = gwas[gwas$V7 < 0.0001,]
  gwas = gwas[,c(1,7,9,10,32)]
  colnames(gwas) = c("SNP","p","b","SE","H4")
  output= merge(gwas,eqtl,by="SNP")

  if(nrow(eqtl)<= 0){
    next
  }
  output = output[order(output$p.y),]
  output = output[!duplicated(output$Gene),]
  egene = rbind(egene,output)
}


write.table(egene,"D:\\AAA\\eqtl\\Artery_coloc_result.txt",col.names = T,row.names = F,quote = F,sep="\t")

 ###############################


load("D:\\AAA\\eqtl\\blood\\Blood_coloc_result.Rdata")
setwd("D:\\AAA\\eqtl\\smr")
egene = NULL
for (i in names(result)) {
  commond = paste("smr-1.3.1-win.exe --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --snp",i," --snp-wind 100, --query 0.0001 --out eqtl_print")
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
write.table(egene,"D:\\AAA\\qtl_scRNA\\Blood_coloc.txt",col.names = T,row.names = F,quote = F,sep="\t")

