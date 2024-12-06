################数据准备#############

##最内部的热图数据##
class=read.table("D:\\bigslice\\bgc_withclass.txt",sep="\t")
numdata=matrix(0,nrow = 26825,ncol = 11)
colnames(numdata)=c("GCF","GCFgroup","GenomesNum","BGCNum","PhylumNum",
                    "ClassNum","OrderNum","FamilyNum","GenusNum","SpeciesNum","StrainNum")
for(i in unique(class$V4)){
  subclass = class[class$V4 %in% i,]
  numdata[i,] =c(i,subclass[1,13],
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
#######数据合并######
phydata = cbind(rownames(phydata),phydata)
classdata = cbind(rownames(classdata),classdata)
data=merge(phydata,classdata,by="V1")
colnames(data)[1]="GCF"
data=merge(numdata,data,by="GCF")
ano = data$GCFgroup
ano[ano %in% c(10,16,17,24,29,30,31,32,34,36,4)] = 0
ano[!(ano %in% 0)] = 1
data = cbind(ano,data)
colnames(data)[1]="Known&Unknown"
data = data[order(data$GCFgroup),]
data  =data[order(data$`Known&Unknown`,decreasing = T),]
write.table(data,"D:\\bigslice\\pictures\\circos data.txt",sep="\t",quote = F,row.names = F,col.names = T)
