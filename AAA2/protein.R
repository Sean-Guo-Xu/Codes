setwd("D:\\AAA\\protein")
pro = read.table("protein_Data.txt",sep="\t",header=F)
id = read.table("id.txt",header=T)
out = merge(id,pro,by.x="From",by.y="V1")
out = out[order(out$V12,decreasing = T),]
out = out[,c(2,13)]
write.table(out,"protein_omics.txt",sep="\t",quote = F,row.names = F,col.names = F)
