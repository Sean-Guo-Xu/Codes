


setwd("D:\\AAA\\eqtl\\smr")
for(i in 1:22){
commond = paste("smr-1.3.1-win.exe --bfile g1000/g1000_eur  --gwas-summary aaagwas.txt --beqtl-summary ./sqtl/sQTL_Artery_Aorta.query.chr",i," --out A_sqtl_",i," --thread-num 10 --diff-freq-prop 0.99",sep="")
system(commond,intern = T)
}







setwd("D:\\AAA\\eqtl\\smr")
alldata=NULL
for (i in 1:22) {
  data=read.table(paste("sqtl_",i,".smr",sep=""),header = T,sep="\t")
  alldata = rbind(alldata,data)
}
alldata = alldata[alldata$p_SMR<0.05,]
champ.import(directory=system.file("extdata",package="ChAMPdata"))
write.table(alldata,"sqtl_Artery_smr.txt",quote=F,sep="\t",row.names = F,col.names = T)

cg=read.table("cg_gene.txt",sep="\t")
alldata = merge(alldata,cg,by.x="probeID",by.y = "V1")
alldata = alldata[!(alldata$V2 %in% ""),]
alldata = alldata[!duplicated(alldata$probeID),]
alldata$Gene = alldata$V2
alldata$V2 =NULL
write.table(alldata,"AAA_mqtl.txt",quote=F,sep="\t",row.names = F,col.names = T)
