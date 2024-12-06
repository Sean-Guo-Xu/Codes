setwd("E:\\bulk-qtl_v8_single-tissue-cis-qtl_GTEx_Analysis_v8_eQTL_EUR")
library(stringr)
gtex = read.table("Muscle_Skeletal.v8.EUR.signif_pairs.txt",sep = "\t",header=T)
locate = as.data.frame(str_split(gtex$variant_id,"_",n=5,simplify = T))
gtex = cbind(locate,gtex)
gtex$V1 = gsub("chr","",gtex$V1)
gtex$pos = paste(gtex$V1,gtex$V2,sep = "_")
out =NULL
for (i in unique(gtex$V1)) {
  part = gtex[gtex$V1 %in% i,]
  file = paste("F:\\SNPG38\\chr",i,".txt",sep="")
  snp = read.table(file,sep="\t")
  part=merge(part,snp,by.x="pos",by.y = "V1")
  out = rbind(out,part) 
}
write.table(out,"Muscle_eqtl.txt",sep="\t",quote = F,col.names = T,row.names = F)
#######
eqtl=read.table("D:\\AAA\\eqtl\\Artery_eqtl.txt",sep = "\t",header = T)
snp=eqtl[eqtl$V3 %in% c("A","T","C","G"),]
snp=snp[snp$V4 %in% c("A","T","C","G"),]
esi = cbind(snp$V1,snp$V2.y,rep(0,nrow(snp)),snp$pos,snp$V3,snp$V4,snp$maf)
