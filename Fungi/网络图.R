setwd("D:\\bigslice\\pictures")
load("GCF_O_distance.Rdata")
class = read.table("D:\\bigslice\\bgc_withclass.txt",sep="\t")
gcfclass = table(class$V4,class$V14)

allclass = NULL
for(i in 1:26825){
  gcf = colnames(gcfclass)[gcfclass[i,] %in% max(gcfclass[i,])]
  gcf = gcf[1]
 allclass=c(allclass,gcf)
  }
gcfclass= cbind(rownames(gcfclass),allclass)
gcfclass=as.data.frame(gcfclass)
out = NULL
for(i in 1:36){
  node = allout[[i]]
  node = node[,1]
  node = class[class$V4 %in% node,]
    count=as.data.frame(table(node$V4))
    count$Var1=as.character(count$Var1)
    count=merge(count,gcfclass,by.x="Var1",by.y="V1")
    count=cbind(count,rep(i,nrow(count)))
    out = rbind(out,count)
}
colnames(out) = c("GCF","Size","Class","Group")
write.table(out,"node.txt",col.names = T,row.names = F,sep="\t",quote = F)

out = NULL
for (i in 1:36) {
  edge= allout[[i]]
  edge=cbind(edge,rep(i,nrow(edge)))
  out = rbind(out,edge)
  }
colnames(out) = c("GCF_1","GCF_2","Distance","Group")
write.table(out,"edge.txt",col.names = T,row.names = F,sep="\t",quote = F)
