setwd("D:\\bigslice\\化合物")
b_c=read.table("bgc-class.txt",header = F)
class=c()
for(i in 1:293926){
  bgc=b_c[which(b_c$V1 %in% i),]
  if (nrow(bgc) == 1){
    if(bgc$V2 %in% c(169,144,145))
    {class=c(class,bgc$V5)}
    else if(bgc$V2 %in% c(172)){class=c(class,"Terpene")}
    else{
      class=c(class,"others")
    }
  }
  else {
    if((169 %in% bgc$V2) & (144 %in% bgc$V2) )
    {class=c(class,"T1PKS&NRPS")}
    else if((169 %in% bgc$V2) & (145 %in% bgc$V2)){
      class = c(class,"T1PKS&NRPS-like")
    }
    else{
      class=c(class,"others")
    }
    }
  }
bgcclass = read.table("D:\\bigslice\\bioinfo_withgcfcluster.txt",sep="\t")

bgcclass=cbind(bgcclass,class)
withclass = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
withclass=withclass[order(withclass$V2),]
withclass$V3=bgcclass[,2]
write.table(bgcclass,"D:\\bigslice\\化合物\\bgc_withclass.txt",row.names = F,col.names=F,quote = F,sep="\t")
data=c("GCF","total","172","144","145","169","others","hybrid","129","112","167","116")
for(j in withclass$V4){
  bgc=bgc[withclass$V4 %in% j,]
  class=c(j,nrow(class))
  for (i in unique(withclass$V3)){
  bgc=bgc[bgc$V3 %in% i, ]
  class=c(class,nrow(bgc))
  }
  data=rbind(data,class)
}

withclass=read.table("bgc_withclass.txt",sep="\t")
withclass=withclass[which(withclass$V3 %in% "hybrid"),]
write.table(withclass,"D:\\bigslice\\化合物\\bgc_with mutiple class.txt",row.names = F,col.names=F,quote = F,sep="\t")
