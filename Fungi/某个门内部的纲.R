library(circlize)
name=read.table("D:\\bigslice\\gcf-rank\\Class.txt",header=T)
classify=read.table("D:\\bigslice\\new11608.txt",sep ="\t")
classify=classify[which(classify$V3 %in% "Ascomycota" ),]
name=name[which(name$Name %in% classify$V4),]
class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
gcf=list()
for (i in name$Name){
  e=class[which(class$V7 %in% i),]
  e=e[!duplicated(e$V4),]
  gcf[[i]]=e$V4
}
mat=c(rep(0,length(name$Name)))
for (i in 1:(length(name$Name)-1)){
  mat=cbind(mat,rep(0,length(name$Name)))}

colnames(mat)=c(1:nrow(mat))
rownames(mat)=c(1:nrow(mat))


for (i in 1:nrow(mat)){
  other=c(0)
  for (j in 1:nrow(mat)){
    t=gcf[[i]][which(gcf[[i]] %in% gcf[[j]])]
    if(i !=j ){
      other=c(other,t)}
    mat[i,j]=length(t)
  }
  other=other[!duplicated(other)]
  mat[i,i]=mat[i,i]-length(other)+1
}

for(i in 1:nrow(mat)){
  for(j in 1:nrow(mat)){
    if(i<j){
      mat[i,j]=0
    }
  }
}
tiff("D:\\bigslice\\Basidiomycotacircle.tiff",width = 3000,height = 3000,units = "px")
chordDiagram(mat, self.link = 1,
             annotationTrack = c( "grid"),
             # 指定显示标签和网格轨道
             # 指定标签和网格轨道的高度
)   
dev.off()

