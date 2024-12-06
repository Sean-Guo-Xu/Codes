library(circlize)
target1="Aspergillaceae"
target2="Nectriaceae"
classify=read.table("D:\\bigslice\\new11608.txt",sep ="\t")
small=classify[which(classify$V6 %in% target1 ),]
cutphy=small[1,3]
class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
name=class[which(class$V6 %in% cutphy),9]
name=name[!duplicated(name)]
gcf=list()
g=c(0)
for (i in name[-which(name %in% c(target1,target2))] ){
  e=class[which(class$V9 %in% i),]
  e=e[!duplicated(e$V4),]
  g=c(g,e$V4)
}
g=g[!duplicated(g)]
g=g[-1]
gcf[["Others"]]=g

g=c(0)
for (i in target1){
  e=class[which(class$V9 %in% i),]
  e=e[!duplicated(e$V4),]
  g=c(g,e$V4)
}
g=g[!duplicated(gcf)]
g=g[-1]
gcf[[target1]]=g

g=c(0)
for (i in target2){
  e=class[which(class$V9 %in% i),]
  e=e[!duplicated(e$V4),]
  g=c(g,e$V4)
}
g=g[!duplicated(gcf)]
g=g[-1]
gcf[[target2]]=g

mat=c(rep(0,3))

mat=cbind(mat,rep(0,3))
mat=cbind(mat,rep(0,3))
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
tiff("D:\\bigslice\\Aspergillaceae&Nectriaceae&others.tiff",width = 3000,height = 3000,units = "px")
chordDiagram(mat, self.link = 1,
             annotationTrack = c( "name","grid"),
             # 指定显示标签和网格轨道
             # 指定标签和网格轨道的高度
)   
dev.off()

