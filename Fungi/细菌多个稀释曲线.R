library(iNEXT)
setwd("D:\\bigslice")
b=read.table("bac.csv",header=T,sep=",")
class = read.table("bgc-gcf400.txt",sep="\t",header=F,fill=T)
memory.limit(size=80000)
bb=b
data=c(1185995)
for(i in c(bb[!duplicated(bb$gcf_id),2])){
  gcf=bb[which(bb$gcf_id == i),]
  data=c(data,nrow(gcf))
}
y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
write.table(y1[["iNextEst"]],"Bac118.txt",sep="\t",quote=F,row.names = F)

for (j in 1:30){
bb=b[sample(1185995,293926),]
data=c(293926)
for(i in c(bb[!duplicated(bb$gcf_id),2])){
  gcf=bb[which(bb$gcf_id == i),]
  data=c(data,nrow(gcf))
}
y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
write.table(y1[["iNextEst"]],paste("Bac29Number-",j,".txt",sep = ""),sep="\t",quote=F,row.names = F)
}

for (j in 1:30){
  bb=b[sample(1185995,100000),]
  data=c(100000)
  for(i in c(bb[!duplicated(bb$gcf_id),2])){
    gcf=bb[which(bb$gcf_id == i),]
    data=c(data,nrow(gcf))
  }
  y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
  write.table(y1[["iNextEst"]],paste("Bac10Number-",j,".txt",sep = ""),sep="\t",quote=F,row.names = F)
}
for (j in 1:30){
  bb=b[sample(1185995,200000),]
  data=c(200000)
  for(i in c(bb[!duplicated(bb$gcf_id),2])){
    gcf=bb[which(bb$gcf_id == i),]
    data=c(data,nrow(gcf))
  }
  y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
  write.table(y1[["iNextEst"]],paste("Bac20Number-",j,".txt",sep = ""),sep="\t",quote=F,row.names = F)
}
for (j in 1:30){
  bb=b[sample(1185995,600000),]
  data=c(600000)
  for(i in c(bb[!duplicated(bb$gcf_id),2])){
    gcf=bb[which(bb$gcf_id == i),]
    data=c(data,nrow(gcf))
  }
  y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
  write.table(y1[["iNextEst"]],paste("Bac60Number-",j,".txt",sep = ""),sep="\t",quote=F,row.names = F)
}


  

for(j in 1:30){
cc=class
data=c[sample(293926,100000),]
for(i in cc[!duplicated(cc$V2),2]){
  gcf=cc[which(cc$V2 == i),]
  data=c(data,nrow(gcf))
}
y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
write.table(y1[["iNextEst"]],paste("Fun10Number-",j,".txt",sep = ""),sep="\t",quote=F,row.names = F)}

for(j in 1:30){
  cc=class
  data=c[sample(293926,200000),]
  for(i in cc[!duplicated(cc$V2),2]){
    gcf=cc[which(cc$V2 == i),]
    data=c(data,nrow(gcf))
  }
  y1=iNEXT(data, datatype="incidence_freq", size=round(seq(1,10000000,2500)), se=FALSE)
  write.table(y1[["iNextEst"]],paste("Fun20Number-",j,".txt",sep = ""),sep="\t",quote=F,row.names = F)}
