setwd("D:\\bigslice\\»¯ºÏÎï")
b_c=read.table("bgc-class.txt",header = F)
h=read.table("bgc_with mutiple class.txt",header=F,sep="\t",stringsAsFactors = F)
newclass=c(unlist(h[1,]),"a","b")
for (i in 1:length(h$V2)){
  bgc=b_c[which(b_c$V1 %in% h[i,2]),]
  cc=c(unlist(h[i,]),paste(bgc$V2,collapse = " "),paste(bgc$V5,collapse = " "))
  newclass=rbind(newclass,cc)
}
newclass=newclass[-1,]
write.table(newclass,"hybrid bgc.txt",sep = "\t",quote=F,row.names = F,col.names = F)
