library("coloc")
library(TwoSampleMR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(locuscomparer)
library(ggplot2)
################GTEX##############
setwd("D:\\AAA\\eqtl")
################Blood#############
memory.limit(size=100000)
options(scipen = 2)
ngwas = 392423  #gwas样本量
neqtl = 670 #eqtl样本量
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
alleqtl = read.table("D:\\AAA\\sqtl\\Blood_sqtl.txt",header=T)
colnames(alleqtl)[19] = "SNP"
alleqtl = alleqtl[order(alleqtl$pval_nominal),]
alleqtl = alleqtl[!duplicated(alleqtl$SNP),]
gwas_ld <- clump_data(gwas, clump_kb=10000, clump_r2=0.2)
result = list()
gwaslist =list()
setwd("D:\\AAA\\eqtl\\smr")

for (i in gwas_ld$SNP) {
  snp = gwas_ld[gwas_ld$SNP %in% i,]
  commond = paste("smr-1.3.1-win.exe --beqtl-summary ./sqtl/sQTL_Whole_Blood.query.chr",snp$chr.outcome," --snp ",i," --query 1 --snp-wind 100 --out eqtl_print",sep="")
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
      write.table(eqtl,"sqtl.tsv",col.names = T,row.names = F,sep="\t",quote = F)
      gwaslist[[i]] = unique(input$V6)
      locuscompare("sqtl.tsv","gwas.tsv",legend =F,snp = i)  
      ggsave(paste(i,"locus.png"),width=10,height = 5)
    }}
  }}
}

save(result,file="D:\\AAA\\sqtl\\Blood_coloc_result.Rdata")
save(gwaslist,file="D:\\AAA\\sqtl\\Blood_coloc_gene.Rdata")

##############基因名注释##################333
egene = NULL
for (i in names(gwaslist)) {
  snp = result[[i]]
  commond = paste("smr-1.3.1-win.exe --beqtl-summary ./sqtl/sQTL_Whole_Blood.query.chr",unique(result[[i]]$V1.x),"  --snp ",i," --query 0.0001 --out eqtl_print",sep="")
  system(commond,intern = T)
  eqtl = read.table("eqtl_print.txt",header=T)
  eqtl  =eqtl[,c(1,14,12,13,7,10)]
  gwas = result[[i]]
  gwas = gwas[gwas$SNP %in% i,]
  gwas = gwas[,c(1,7,9,10)]
  colnames(gwas) = c("SNP","p","b","SE")
  output= merge(gwas,eqtl,by="SNP")
  egene = rbind(egene,output)
}
write.table(egene,"D:\\AAA\\sqtl\\Blood_coloc_result.txt",col.names = T,row.names = F,quote = F,sep="\t")

write.table(egene,"D:\\AAA\\sqtl\\Artery_coloc_result.txt",col.names = T,row.names = F,quote = F,sep="\t")

 ###############################
