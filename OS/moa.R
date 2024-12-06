setwd("D:\\生信\\target数据分析\\药学")
drug<-read.table("top20.txt",sep = "\t",header=T)
drug[,3]<-as.character(drug[,3])
drug[,2]<-as.character(drug[,2])

  c<-c("1","2","3")
  for (i in 1:length(drug[,1])) {
    c1<-strsplit(drug[i,3],",")
    for(j in c1){
    c2<-rep(drug[i,2],length(c1))
    c3<-rep(drug[i,1],length(c1))
    cc<-cbind(j,c2,c3)
    c<-rbind(c,cc)
    }}
  c<-c[-1,]
  write.table (c, file ="drugfuc.txt", sep ="\t", row.names =F, col.names =F, quote =FALSE)
  