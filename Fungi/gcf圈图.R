setwd("D:\\bigslice")
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
class=read.table("D:\\bigslice\\bgc_withclass.txt",sep="\t")
numdata=matrix(0,nrow = 26825,ncol = 10)
colnames(numdata)=c("GCF","GenomesNum","BGCNum","PhylumNum",
                    "ClassNum","OrderNum","FamilyNum","GenusNum","SpeciesNum","StrainNum")
for(i in unique(class$V4)){
  subclass = class[class$V4 %in% i,]
  numdata[i,] =c(i,
                 length(unique(subclass$V1)),
                 length(unique(subclass$V2)),
                 length(unique(subclass$V6)),
                 length(unique(subclass$V7)),
                 length(unique(subclass$V8)),
                 length(unique(subclass$V9)),
                 length(unique(subclass$V10)),
                 length(unique(subclass$V11)),
                 length(unique(subclass$V12)))
}

##化合物类型数据准备########
classdata = table(class$V4,class$V14)
classdata = as.matrix(classdata)
sum = apply(classdata, 1, sum)
for (i in 1:7) {
classdata[,i] =   classdata[,i]/sum
}
#####门数据准备########
phydata = table(class$V4,class$V6)
phydata = as.matrix(phydata)
sum = apply(phydata, 1, sum)
for (i in 1:10) {
  phydata[,i] =   phydata[,i]/sum
}
#####order######
order= read.table("圈图顺序.txt",sep="\t",header = T)
level  = unique(order$type)
order  =order[order(order$GCF_id),]

data = cbind(numdata,classdata,phydata,order$type)
data =as.data.frame(data)
colnames(data)[28]="Order"
data$Order = factor(data$Order,levels = level)
save(data,file="circros data.Rdata")

###############绘图####################

data = data[order(data$BGCNum,decreasing=T),]
data = data[order(data$Order),]
phycol=read.table("color_phy.txt",sep = "\t", comment.char = "")
name=phycol$V1
phycol = phycol$V2
names(phycol) = name

clacol=read.table("color_class.txt",sep = "\t", comment.char = "")
name=clacol$V1
clacol = clacol$V2
names(clacol) = name
data[,1:27]  =apply(data[,1:27],2,as.numeric)

data$BGCNum = log10(data$BGCNum+1)
data$GenomesNum = log10(data$GenomesNum+1)
data$StrainNum = log10(data$StrainNum+1)
data$SpeciesNum = log2(data$SpeciesNum+1)
data$GenusNum = log2(data$GenusNum+1)
data$FamilyNum = log2(data$FamilyNum+1)
data$OrderNum = log2(data$OrderNum+1)

numcol= pal_frontiers()(9)
names(numcol) =  c("PhylumNum","ClassNum","OrderNum","FamilyNum","GenusNum","SpeciesNum","StrainNum","GenomesNum","BGCNum")
circos.clear()
creatwhite = function(){
  white = "white"
  names(white) = "empty"
  circos.heatmap(rep("empty",26825), col = white, track.height = 0.05,track.margin=c(0,0))
}



circos.clear()
pdf("circos.pdf",width = 10,height = 10)
circos.par(gap.after=0)
circos.heatmap.initialize(data,cluster = F)

for(i in names(phycol)){
  col = colorRamp2(c(0,1),c("white",phycol[i]))
  circos.heatmap(data[,i], col = col, track.height = 0.01,track.margin=c(0,0.005))
}
creatwhite()
for(i in names(clacol)){
  col = colorRamp2(c(0,1),c("white",clacol[i]))
  circos.heatmap(data[,i], col = col, track.height = 0.01,track.margin=c(0,0.005))
}
creatwhite()

for (i in names(numcol)) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  circos.heatmap(data[,i], col = col, track.height = 0.02,track.margin=c(0,0.005))
}

dev.off()

phyleg=Legend(title = "Phylum", at = names(phycol), 
              legend_gp = gpar(fill = phycol))
classleg = Legend(title = "Type", at = names(clacol), 
                  legend_gp = gpar(fill = clacol))
listleg = list()
circleorder = names(numcol)
for (i in circleorder[1:2]) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  listleg[[i]] =Legend(title = i, col_fun = col )
}
for (i in circleorder[3:6]) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  listleg[[i]] =Legend(title = paste("log2(",i,")",sep = ""), col_fun = col )
}
for (i in circleorder[7:9]) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  listleg[[i]] =Legend(title = paste("log10(",i,")",sep = ""), col_fun = col )
}
h = dev.size()[2]
lgd_list = packLegend(phyleg,classleg,listleg[[1]],listleg[[2]],listleg[[3]]
                      ,listleg[[4]],listleg[[5]],listleg[[6]],listleg[[7]],listleg[[8]],listleg[[9]],max_height = unit(0.9*h, "inch"))
pdf("circle legend.pdf",width = 5,height = 10)
draw(lgd_list)
dev.off()
