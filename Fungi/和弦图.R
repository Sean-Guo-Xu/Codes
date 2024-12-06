library(circlize)
library(ggsci)
setwd("D:\\bigslice")
genus  =c("Aspergillus","Fusarium")
  class = read.table("bgc_withclass.txt",sep="\t")
gcf = list()

phy = unique(class$V6)
phy = phy[!(phy %in% class[class$V10 %in% genus,6])]
phy = phy[!(phy %in% "Others")]
for (i in genus ){
  gcf[[i]] = unique(class[class$V10 %in% i,4])  
}
for (i in phy) {
  gcf[[i]] = unique(class[class$V6 %in% i,4])  
}
mat = matrix(0,nrow = length(gcf),ncol = length(gcf))

phycol=read.table("color_phy.txt",sep = "\t", comment.char = "")
phycol = phycol[phycol$V1 %in% names(gcf),]
gencol=pal_d3()(10)
names( gencol)=c("Aspergillus","Fusarium","Xylaria","Hypoxylon","Penicillium","Colletotrichum","Talaromyces","Diaporthe","Nemania","Calonectria")
gencol = gencol[names(gencol) %in% genus]
col = c(gencol,c(phycol$V2))
names(col) = c(names(gencol),phycol$V1)
colnames(mat) = names(col)
rownames(mat) = names(col)
for (i in names(col)) {
   g = NULL
  for (j in names(col)) {
    mat[i,j] = length(gcf[[i]][gcf[[i]] %in% gcf[[j]]])
    if(i != j){
    g  = c(g,gcf[[i]][gcf[[i]] %in% gcf[[j]]])}
  }
   mat[i,i] = length(gcf[[i]][! (gcf[[i]] %in% g)])
}

for(i in 1:nrow(mat)){
  for(j in 1:nrow(mat)){
    if(i<j){
      mat[i,j]=0
    }
  }
}


colmat = matrix(0,nrow = length(gcf),ncol = length(gcf))
colnames(colmat) = names(gcf)
rownames(colmat) = names(gcf)
for (i in names(gcf)) {
  colmat[,i] = col
}
pdf(paste(length(genus),"个属门外比较.pdf"),width = 5,height=5)
chordDiagram(mat, self.link = 1,grid.col =col ,
             annotationTrack = c( "grid"),
)   
dev.off()
names(gcf)
#######################################################
library(circlize)
library(ggsci)
setwd("D:\\bigslice")
gcf = list()

phy = unique(class$V6)
phy = phy[(phy %in% class[class$V10 %in% genus,6])]

for (i in genus ){
  gcf[[i]] = unique(class[class$V10 %in% i,4])  
}
for (i in phy) {
  subclass = class[!(class$V10 %in% genus), ]
  gcf[[i]] = unique(subclass[subclass$V6 %in% i,4])  
}
mat = matrix(0,nrow = length(gcf),ncol = length(gcf))

phycol=read.table("color_phy.txt",sep = "\t", comment.char = "")
phycol = phycol[phycol$V1 %in% names(gcf),]
gencol=pal_d3()(10)
names( gencol)=c("Aspergillus","Fusarium","Xylaria","Hypoxylon","Penicillium","Colletotrichum","Talaromyces","Diaporthe","Nemania","Calonectria")
gencol = gencol[names(gencol) %in% genus]
col = c(gencol,c(phycol$V2))
names(col) = c(names(gencol),phycol$V1)
colnames(mat) = names(col)
rownames(mat) = names(col)
for (i in names(col)) {
  g = NULL
  for (j in names(col)) {
    mat[i,j] = length(gcf[[i]][gcf[[i]] %in% gcf[[j]]])
    if(i != j){
      g  = c(g,gcf[[i]][gcf[[i]] %in% gcf[[j]]])}
  }
  mat[i,i] = length(gcf[[i]][! (gcf[[i]] %in% g)])
}

for(i in 1:nrow(mat)){
  for(j in 1:nrow(mat)){
    if(i<j){
      mat[i,j]=0
    }
  }
}


colmat = matrix(0,nrow = length(gcf),ncol = length(gcf))
colnames(colmat) = names(gcf)
rownames(colmat) = names(gcf)
for (i in names(gcf)) {
  colmat[,i] = col
}
pdf(paste(length(genus),"个属门内比较.pdf"),width = 5,height=5)
chordDiagram(mat, self.link = 1,grid.col =col ,
             annotationTrack = c( "grid"),
)   
dev.off()
names(gcf)


#################十个属#############333

library(circlize)
library(ggsci)
setwd("D:\\bigslice")
gcf = list()
genus=c("Aspergillus","Fusarium","Xylaria","Hypoxylon","Penicillium","Colletotrichum","Talaromyces","Diaporthe","Nemania","Calonectria")
for (i in genus ){
  gcf[[i]] = unique(class[class$V10 %in% i,4])  
}

mat = matrix(0,nrow = length(gcf),ncol = length(gcf))


gencol=pal_d3()(10)
names( gencol)=c("Aspergillus","Fusarium","Xylaria","Hypoxylon","Penicillium","Colletotrichum","Talaromyces","Diaporthe","Nemania","Calonectria")
gencol = gencol[names(gencol) %in% genus]
col = c(gencol)
names(col) = c(names(gencol))
colnames(mat) = names(col)
rownames(mat) = names(col)
for (i in names(col)) {
  g = NULL
  for (j in names(col)) {
    mat[i,j] = length(gcf[[i]][gcf[[i]] %in% gcf[[j]]])
    if(i != j){
      g  = c(g,gcf[[i]][gcf[[i]] %in% gcf[[j]]])}
  }
  mat[i,i] = length(gcf[[i]][! (gcf[[i]] %in% g)])
}

for(i in 1:nrow(mat)){
  for(j in 1:nrow(mat)){
    if(i<j){
      mat[i,j]=0
    }
  }
}


colmat = matrix(0,nrow = length(gcf),ncol = length(gcf))
colnames(colmat) = names(gcf)
rownames(colmat) = names(gcf)
for (i in names(gcf)) {
  colmat[,i] = col
}
pdf(paste(length(genus),"个属比较.pdf"),width = 5,height=5)
chordDiagram(mat, self.link = 1,grid.col =col ,
             annotationTrack = c( "grid"),
)   
dev.off()
names(gcf)
