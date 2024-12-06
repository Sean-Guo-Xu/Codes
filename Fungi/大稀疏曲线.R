library("iNEXT")
library(ggplot2)
library(ggthemes)
library(ggsci)
class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
data=c(11097)
for(i in 1:26825){
  gcf=class[which(class$V4 == i),]
  gcf=gcf[!duplicated(gcf$V1),]
  data=c(data,nrow(gcf))
}

y=iNEXT(data, datatype="incidence_freq", size=round(seq(1,150000,500)), se=T)
ggiNEXT(y,type = 1,se=T)+theme_stata()+xlab("Count of Genomes")+ylab("Count of GCFs")+scale_color_npg()
ggsave("D:\\bigslice\\all GCFs curve.tiff",width = 4,height = 6)
write.table(y$iNextEst,"D:\\bigslice\\data\\curve.txt")
