library(Seurat)
library(org.Hs.eg.db)
setwd("E:\\AAA_scRNA\\AAA")
dir = list.dirs()
pbmc=Read10X(dir[2],gene.column = 1)
gene = rownames(pbmc)
gene = gsub("[.].*","",gene)
gene = mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL" , keytype ="ENSEMBL")
rownames(pbmc) = gene
pbmc = pbmc[!is.na(rownames(pbmc)),]
pbmc=limma::avereps(pbmc)
pbmc=CreateSeuratObject(pbmc)
pbmc$sample = gsub("./", "", dir[2])
allgene = gene
for (i in 3:23) {
  part = Read10X(dir[i],gene.column = 1)
  gene = rownames(part)
  gene = gsub("[.].*","",gene)
  gene = mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL" , keytype ="ENSEMBL")
  rownames(part) = gene
  part = part[!is.na(rownames(part)),]
  part=limma::avereps(part)
  part=CreateSeuratObject(part)
  allgene= intersect(gene,allgene)
  part$sample = gsub("./", "", dir[i])
  pbmc = merge(pbmc,part)
}
pbmc = subset(pbmc,features = allgene)
save(pbmc,file="22_AAA.Rdata")

setwd("E:\\AAA_scRNA\\Control")
dir = list.dirs()
pbmc=Read10X(dir[2],gene.column = 1)
gene = rownames(pbmc)
gene = gsub("[.].*","",gene)
gene = mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL" , keytype ="ENSEMBL")
rownames(pbmc) = gene
pbmc = pbmc[!is.na(rownames(pbmc)),]
pbmc=limma::avereps(pbmc)
pbmc=CreateSeuratObject(pbmc)
pbmc$sample = gsub("./", "", dir[2])
allgene = gene
for (i in 3:12) {
  part = Read10X(dir[i],gene.column = 1)
  gene = rownames(part)
  gene = gsub("[.].*","",gene)
  gene = mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL" , keytype ="ENSEMBL")
  rownames(part) = gene
  part = part[!is.na(rownames(part)),]
  part=limma::avereps(part)
  part=CreateSeuratObject(part)
  allgene= intersect(gene,allgene)
  part$sample = gsub("./", "", dir[i])
  pbmc = merge(pbmc,part)
}
pbmc = subset(pbmc,features = allgene)
save(pbmc,file="12_Control.Rdata")