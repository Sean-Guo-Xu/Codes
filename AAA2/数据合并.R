setwd("D:\\AAA\\omics")
library(RobustRankAggreg)
set4=read.table("47472.txt",sep='\t',header=T)
set5=read.table("57691.txt",sep='\t',header=T)
set7=read.table("7084.txt",sep='\t',header=T)
setlist=list("set4"=set4,"set5"=set5,"set6"=set7)
uplist=list()
downlist=list()
for (i in names(setlist) ) {
  set=setlist[[i]]
  set_up=set[set$logFC>0,]
  set_down=set[which(set$logFC<(0)),]
  set_up=set_up[order(set_up$logFC,decreasing = T),]
  set_down=set_down[order(set_down$logFC),]
  uplist[[i]]=set_up$Gene
  downlist[[i]]=set_down$Gene
}
downre=RobustRankAggreg::aggregateRanks(downlist)
upre=RobustRankAggreg::aggregateRanks(uplist)
downre[,3] = downre$Score
upre[,3] = upre$Score
downre[,2] = rep(-1,nrow(downre))
upre[,2] = rep(1,nrow(upre))
downre=downre[downre$Score<0.05,]
upre = upre[upre$Score<0.05,]

relog=rbind(upre,downre)
write.table(relog,"3GEOforBETA.txt",row.names = F,col.names = F,sep = "\t",quote = F)
