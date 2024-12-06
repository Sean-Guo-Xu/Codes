setwd("D:\\bigslice")
np=read.table("NPAltas Ù–≈œ¢.txt",header=T,sep="\t",fill=T)
ge=read.table("764result.txt",sep="\t",fill=T)
clustnum=c(0)
compdnum=c(0)
for (i in 1:nrow(ge)){
  genus=np[which(np$genus %in% ge[i,1]),]
  compdnum=c(compdnum,nrow(genus[!duplicated(genus$npaid),]))
  clustnum=c(clustnum,nrow(genus[!duplicated(genus$compound_cluster_id),]))
}
clustnum=clustnum[-1]
compdnum=compdnum[-1]
ge = cbind(ge,compdnum,clustnum)
phy=ge[!duplicated(ge$V9),9]
cla=ge[!duplicated(ge$V8),8]
cla=cla[-(cla== "")]
c1=c("Name")
c2=c("CompondNum")
c3=c("ClusterNum")
for(i in phy){
  phynum=ge[which(ge$V9 == i),]
  c1=c(c1,i)
  c2=c(c2,sum(phynum[!duplicated(phynum$compdnum),11]))
  c3=c(c3,sum(phynum[!duplicated(phynum$clustnum),12]))
}
phy=cbind(c1,c2,c3)

c1=c("Name")
c2=c("CompondNum")
c3=c("ClusterNum")
for(i in cla){
  phynum=ge[which(ge$V8 == i),]
  c1=c(c1,i)
  c2=c(c2,sum(phynum[!duplicated(phynum$compdnum),11]))
  c3=c(c3,sum(phynum[!duplicated(phynum$clustnum),12]))
}
cla=cbind(c1,c2,c3)
write.table(phy,"phyCompound.txt",row.names = F,col.names = F,sep = "\t",quote=F)

write.table(cla,"claCompound.txt",row.names = F,col.names = F,sep = "\t",quote=F)
