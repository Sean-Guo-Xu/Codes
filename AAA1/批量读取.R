setwd("E:\\AAA_scRNA\\N")
fs=list.files('./','^GSM')
library(stringr)
samples=str_split(fs,'_',simplify = T)[,1]
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste(str_split(y[1],'_',simplify = T)[,2:3],collapse = '')
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"features.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})
folders=list.files('./')
library(Seurat)

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = gsub("barcodes.tsv.gz","",folder) )
})
genelist = sceList[[i]]@assays$RNA@counts@Dimnames[[1]]
for (i in 1:4) {
  sceList[[i]]$sample = rep(paste(folders[i],sep = ""),length(sceList[[i]]$orig.ident))
  genelist = genelist[genelist %in% sceList[[i]]@assays$RNA@counts@Dimnames[[1]]] 
}


pbmc <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]]))

pbmc = subset(pbmc,features = genelist )
save(pbmc,file="mouse_aaa_all.Rdata")

setwd("E:\\AAA_scRNA\\GSE166676")
fs=list.files('./','^GSM')
library(stringr)
samples=str_split(fs,'_',simplify = T)[,1]
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste(str_split(y[1],'[.]',simplify = T)[,1],collapse = '')
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"features.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})
folders=list.files('./')
library(Seurat)

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = gsub("barcodes.tsv.gz","",folder) )
})
for (i in 1:5) {
  sceList[[i]]$sample = rep(paste("Normal_",i,sep = ""),length(sceList[[i]]$orig.ident)) 
}
pbmc <- merge(sceList[[1]], 
              y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]]))
save(pbmc,file="M_Normal.Rdata")


setwd("E:\\AAA_scRNA\\GSE226492")
fs=list.files('./','^GSM')
library(stringr)

folders=list.files('./')
library(Seurat)

sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = paste("GSE226492",folder,sep = ""))
})
for (i in 1:6) {
  sceList[[i]]$sample = rep(paste("GSE226492_",i,sep = ""),length(sceList[[i]]$orig.ident)) 
}
pbmc <- merge(sceList[[4]], 
              y = c(sceList[[5]],sceList[[6]]))
save(pbmc,file="GSE226492_N.Rdata")
load("E:\\AAA_scRNA\\6Normal.Rdata")
pbmc$orig.ident[pbmc$orig.ident %in% "Normal1"]="Control4"
pbmc$orig.ident[pbmc$orig.ident %in% "Normal2"]="Control5"
pbmc$orig.ident[pbmc$orig.ident %in% "Normal3"]="Control6"
aaa1=pbmc
setwd("E:\\AAA_scRNA")
load("GSE166676_N.Rdata")
pbmc$orig.ident[pbmc$orig.ident %in% "GSM5077731_AAA20190703Tissue"]="GSE166676_1"
pbmc$orig.ident[pbmc$orig.ident %in% "GSM5077732_AAA20191127"]="GSE166676_2"
aaa3=pbmc
load("GSE226492_N.Rdata")
pbmc$orig.ident[pbmc$orig.ident %in% "GSE2264924"]="GSE226492_1"
pbmc$orig.ident[pbmc$orig.ident %in% "GSE2264925"]="GSE226492_2"
pbmc$orig.ident[pbmc$orig.ident %in% "GSE2264926"]="GSE226492_3"
aaa4=pbmc

pbmc = merge(aaa1,aaa3)
pbmc=merge(pbmc,aaa4)
unique(pbmc$orig.ident)
save(pbmc,file = "E:\\AAA_scRNA\\all_Normal.Rdata")
