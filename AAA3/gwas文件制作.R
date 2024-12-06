#############finnen###############
gwas = read.table("D:\\AAA\\eqtl\\smr\\aaagwas.txt")
setwd("E:\\AAA_gwas")
gwas = read.table("finngen_R10_I9_ABAORTANEUR",sep="\t")
gwas = gwas[,c(5,4,3,11,9,10,7)]
gwas = cbind(gwas,rep(NA,nrow(gwas)))
colnames(gwas) = c("SNP","A1","A2","freq","b","se","p","n")
gwas = gwas[gwas$SNP != "", ]
gwas = gwas[gwas$A1 %in% c("A","T","G","C"),]
gwas = gwas[gwas$A2 %in% c("A","T","G","C"),]
gwas$SNP = gsub(",.+","",gwas$SNP) ##delete multiple SNPs
write.table(gwas,"Europe_AAA_gwas.txt",sep = "\t",row.names = F,col.names = T,quote = F)
