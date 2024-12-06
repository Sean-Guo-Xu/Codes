library(circlize)
name=read.table("D:\\bigslice\\gcf-rank\\Class.txt",header=T)
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
save(mat,file = "D:\\bigslice\\data\\class_circos_bigslice.Rdata")
save(class,file = "D:\\bigslice\\data\\bgc_all_information.Rdata")
tiff("D:\\bigslice\\phycircle.tiff",width = 3000,height = 3000,units = "px")
chordDiagram(mat, self.link = 1,
             annotationTrack = c( "grid"),
)   
dev.off()
