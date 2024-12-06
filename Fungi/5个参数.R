for (j in c(3,4,6,7)) {
  

gcf=read.table(paste("D:\\bigslice\\bgc-gcf",j,"50.txt",sep=""))
family=read.table("D:\\bigslice\\gcf-rank\\Family.txt",header=T,sep="\t")
family=family[1:11,]
family=family[-3,]
genus=read.table("D:\\bigslice\\gcf-rank\\Genus.txt",header=T,sep="\t")
genus=genus[1:11,]
genus=genus[-6,]
classify=read.table("D:\\bigslice\\new11608.txt",header=F,sep = "\t")
colnames(gcf)[3]="genome"
colnames(classify)[1]="genome"
bgc=merge(gcf,classify,by="genome")
genustop=c(0)
genusname=c("1",genus[1:10,1])
for (i in 1:10){
  newg=bgc[which(bgc$V7 %in% genus[i,1]),]
  newg=newg[!duplicated(newg$V2.x),]
  genustop=c(genustop,nrow(newg))
}
genustop=cbind(genusname[2:11],genustop[2:11])
colnames(genustop)=c("Genus","GCFNumber")                                                                                                                                   
write.table(genustop,paste("D:\\bigslice\\",j,"50genus.txt",sep=""),sep="\t",quote = F,row.names = F)

genustop=c(0)
genusname=c("1",family[1:10,1])
for (i in 1:10){
  newg=bgc[which(bgc$V6 %in% family[i,1]),]
  newg=newg[!duplicated(newg$V2.x),]
  genustop=c(genustop,nrow(newg))
}
genustop=cbind(genusname[2:11],genustop[2:11])
colnames(genustop)=c("Family","GCFNumber")
write.table(genustop,paste("D:\\bigslice\\",j,"50family.txt",sep=""),sep="\t",quote = F,row.names = F)
}
