library(circlize)
library(ggsci)
library(ComplexHeatmap)
setwd("D:\\bigslice\\pictures")
genus  = read.table("D:\\bigslice\\gcf-rank\\genus.txt",sep="\t",header = T)
genus = genus[1:11,1]
genus = genus[-6]


class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
gcf = NULL
for (i in unique(class$V10)) {
  part = class[class$V10 %in% i,]
  part = part[!duplicated(part$V4),]
  gcf = rbind(gcf,c(part[1,6],i,nrow(part)))
}
group = gcf
group=as.data.frame(group)
group = group[group$V2 %in% genus,]
gcf=list()
for (i in group$V2){
  e=class[which(class$V10 %in% i),]
  e=e[!duplicated(e$V4),]
  gcf[[i]]=e$V4
}
othergcf = class[!(class$V10 %in% genus),]
othergcf =unique(othergcf$V4)
gcf[["Others"]]=othergcf
col = pal_frontiers()(10)
col = c(col,"gray")
group = rbind(group,c("oo","Others",length(othergcf)))
genus=c(genus,"Others")

mat=c(rep(0,length(genus)))
for (i in 1:(length(group$V2)-1)){
  mat=cbind(mat,rep(0,length(genus)))}

colnames(mat)=genus
rownames(mat)=genus


for (i in genus){
  other=c(0)
  for (j in genus){
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

pdf("10 genus circle.pdf",width = 8,height = 8)
chordDiagram(mat,
             grid.col=col,
             annotationTrack = c("grid"),
             self.link = 1,
             preAllocateTracks = list(
               track.height = mm_h(4),
               track.margin = c(mm_h(4), 0)
             ))
dev.off()
leg=Legend(title = "Genus", at = genus, 
              legend_gp = gpar(fill = col))
pdf("10 genus legend.pdf",width = 8,height = 5)
draw(leg)
dev.off()
