setwd("E:\\AAA_scRNA")
library(Seurat)
load("8AAA.Rdata")
aaa=pbmc
load("6Normal.Rdata")
pbmc=merge(aaa,pbmc)
rm(aaa)
markerlist = list()
for (i in unique(pbmc$cell)[1:6]) {
  part = subset(pbmc,cell== i)
  marker = FindMarkers(part,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "sample",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
  markerlist[[i]] = marker 
  }
markerlist[["all"]] = marker 
marker = FindMarkers(pbmc,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "aaa_ro_con",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
markerlist[["all"]] = marker 
save(diff,file="diff_all.Rdata")
diff=list()
for (i in names(markerlist)) {
  diff[[i]]  = cbind(markerlist[[i]],rownames(markerlist[[i]]))
    }
