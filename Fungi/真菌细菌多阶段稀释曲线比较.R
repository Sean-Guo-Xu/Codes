path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac10Number-",1,".txt",sep = "")
library(ggplot2)
bac=read.table(path,sep = "\t")
data=bac$V1
for (i in 1:30){
  path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac10Number-",i,".txt",sep = "")
  bac=read.table(path,sep = "\t")
  bac=bac$V2
  bac=as.numeric(bac)
  data=cbind(data,bac)
}
data=data[-1,]
var=c()
ave=c()
for (i in 1:length(data[,1])) {
  
  s=mean(as.numeric(data[i,2:31]))
  v=sd(as.numeric(data[i,2:31]))
  ave=c(ave,s)
  var=c(var,v)
}
data=cbind(data[,1],ave,var)
colnames(data)[1]="x"
data=as.data.frame(data)
data$x=as.numeric(data$x)
data$ave=as.numeric(data$ave)
data$var=as.numeric(data$var)
data1=data
###################################

path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac20Number-",1,".txt",sep = "")
library(ggplot2)
bac=read.table(path,sep = "\t")
data=bac$V1
for (i in 1:30){
  path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac20Number-",i,".txt",sep = "")
  bac=read.table(path,sep = "\t")
  bac=bac$V2
  bac=as.numeric(bac)
  data=cbind(data,bac)
}
data=data[-1,]
var=c()
ave=c()
for (i in 1:length(data[,1])) {
  
  s=mean(as.numeric(data[i,2:31]))
  v=sd(as.numeric(data[i,2:31]))
  ave=c(ave,s)
  var=c(var,v)
}
data=cbind(data[,1],ave,var)
colnames(data)[1]="x"
data=as.data.frame(data)
data$x=as.numeric(data$x)
data$ave=as.numeric(data$ave)
data$var=as.numeric(data$var)
data2=data

#########################################

path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac20Number-",1,".txt",sep = "")
library(ggplot2)
bac=read.table(path,sep = "\t")
data=bac$V1
for (i in 1:30){
  path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac20Number-",i,".txt",sep = "")
  bac=read.table(path,sep = "\t")
  bac=bac$V2
  bac=as.numeric(bac)
  data=cbind(data,bac)
}
data=data[-1,]
var=c()
ave=c()
for (i in 1:length(data[,1])) {
  
  s=mean(as.numeric(data[i,2:31]))
  v=sd(as.numeric(data[i,2:31]))
  ave=c(ave,s)
  var=c(var,v)
}
data=cbind(data[,1],ave,var)
colnames(data)[1]="x"
data=as.data.frame(data)
data$x=as.numeric(data$x)
data$ave=as.numeric(data$ave)
data$var=as.numeric(data$var)
data2=data
########################################
path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Fun10Number-",1,".txt",sep = "")
library(ggplot2)
bac=read.table(path,sep = "\t")
data=bac$V1
for (i in 1:30){
  path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Fun10Number-",i,".txt",sep = "")
  bac=read.table(path,sep = "\t")
  bac=bac$V2
  bac=as.numeric(bac)
  data=cbind(data,bac)
}
data=data[-1,]
var=c()
ave=c()
for (i in 1:length(data[,1])) {
  
  s=mean(as.numeric(data[i,2:31]))
  v=sd(as.numeric(data[i,2:31]))
  ave=c(ave,s)
  var=c(var,v)
}
data=cbind(data[,1],ave,var)
colnames(data)[1]="x"
data=as.data.frame(data)
data$x=as.numeric(data$x)
data$ave=as.numeric(data$ave)
data$var=as.numeric(data$var)
data4=data
##############################################33
path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Fun20Number-",1,".txt",sep = "")
library(ggplot2)
bac=read.table(path,sep = "\t")
data=bac$V1
for (i in 1:30){
  path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Fun20Number-",i,".txt",sep = "")
  bac=read.table(path,sep = "\t")
  bac=bac$V2
  bac=as.numeric(bac)
  data=cbind(data,bac)
}
data=data[-1,]
var=c()
ave=c()
for (i in 1:length(data[,1])) {
  
  s=mean(as.numeric(data[i,2:31]))
  v=sd(as.numeric(data[i,2:31]))
  ave=c(ave,s)
  var=c(var,v)
}
data=cbind(data[,1],ave,var)
colnames(data)[1]="x"
data=as.data.frame(data)
data$x=as.numeric(data$x)
data$ave=as.numeric(data$ave)
data$var=as.numeric(data$var)
data5=data
#################################33
path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac29Number-",1,".txt",sep = "")
library(ggplot2)
bac=read.table(path,sep = "\t")
data=bac$V1
for (i in 1:30){
  path=paste("D:\\bigslice\\化合物稀释曲线\\bigslice\\Bac29Number-",i,".txt",sep = "")
  bac=read.table(path,sep = "\t")
  bac=bac$V2
  bac=as.numeric(bac)
  data=cbind(data,bac)
}
data=data[-1,]
var=c()
ave=c()
for (i in 1:length(data[,1])) {
  
  s=mean(as.numeric(data[i,2:31]))
  v=sd(as.numeric(data[i,2:31]))
  ave=c(ave,s)
  var=c(var,v)
}
data=cbind(data[,1],ave,var)
colnames(data)[1]="x"
data=as.data.frame(data)
data$x=as.numeric(data$x)
data$ave=as.numeric(data$ave)
data$var=as.numeric(data$var)
data6=data

data3 = read.table("D:\\bigslice\\化合物稀释曲线\\bigslice\\Fun29.txt",header=T,sep = "\t")

data=list(data1,data2,data3,data4,data5,data6)
for(i in 1:6){ 
  p=data[[i]]
  p = p[p[,1] < 3000000,]
  data[[i]] = p 
}
######################################
library(ComplexHeatmap)
library(ggsci)
library(ggthemes)

col = c("#FF8888FF","#EE3F3FFF" ,"#D60C00FF","#C6C0FFFF","#6B58EEFF","#4500ACFF" )
names(col)=c("Fungi 100k","Fungi 200k","Fungi 293k","Bacteria 100k","Bacteria 200k","Bacteria 293k")
p=ggplot() +
  geom_line(data=data[[2]],aes(x=x, y=ave),lwd=1,color=col[5]) +geom_ribbon(data=data[[2]],aes(ymin=ave-var, ymax=ave+var, x=x), fill = "#D4D4FFFF", alpha = 0.3)+  
  geom_line(data=data[[1]],aes(x=x, y=ave),lwd=1,color=col[4])+geom_ribbon(data=data[[1]],aes(ymin=ave-var, ymax=ave+var, x=x), fill = "#D4D4FFFF", alpha = 0.3)+    
  geom_line(data=data[[5]],aes(x=x, y=ave),lwd=1,color=col[1]) +geom_ribbon(data=data[[5]],aes(ymin=ave-var, ymax=ave+var, x=x), fill = "#FFBFE5FF", alpha = 0.3)+  
  geom_line(data=data[[6]],aes(x=x, y=ave),lwd=1,color=col[6]) +geom_ribbon(data=data[[6]],aes(ymin=ave-var, ymax=ave+var, x=x), fill = "#D4D4FFFF", alpha = 0.3)+  
  geom_line(data=data[[4]],aes(x=x, y=ave),lwd=1,color=col[2]) +geom_ribbon(data=data[[4]],aes(ymin=ave-var, ymax=ave+var, x=x), fill = "#FFBFE5FF", alpha = 0.3)+
  geom_line(data=data[[3]],aes(x=t, y=qD),lwd=1,color=col[3])+theme_stata()+xlab("BGC/CompoundNumber")+ylab("ClusterNumber")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
ggsave("D:\\bigslice\\compare.tiff",height = 8,width = 6)

leg = Legend(title = "Curve", at = names(col), 
                  legend_gp = gpar(fill = col))
pdf("D:\\bigslice\\legend.pdf",width=5,height=5)
draw(leg)
dev.off()
###1,2,6细菌###4,5,3