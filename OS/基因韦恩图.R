library(VennDiagram)
setwd("E:\\bladder")
load("survgene.Rdata")
load("coxgene.Rdata")
load("diffgene.Rdata")
survgene = c(survgene[[1]],survgene[[2]],survgene[[3]])
coxgene = c(coxgene[[1]]$id,coxgene[[2]]$id,coxgene[[3]]$id)
diffgene = c(diffgene[[1]]$ID,diffgene[[2]]$ID)

survgene = unique(survgene)
coxgene = unique(coxgene)
diffgene = unique(diffgene)
genelist = list(
  DEG = diffgene,
  CG = coxgene,
  SG = survgene
)
col = pal_frontiers()(3)
names(col) = c("DEG","CG","SG")
venn= venn.diagram(genelist,compression = "lzw",filename = "venn.tiff",imagetype="tiff",fill=col,height = 1500,width = 1500, cat.default.pos = "outer")

center = survgene[survgene %in% coxgene]
center = center[center %in% diffgene]
