setwd("E:\\bladder\\infiltration")     #???ù???Ŀ¼                                       #?????ļ?
gmtFile="immune.gmt"                                           #GMT?ļ?
gene = "EFEMP1"
imgene = c("CTLA4","CD274","PDCD1","TIGIT","PDCD1LG2")
library(GSEABase)
library(GSVA)
library(limma)
vm = read.table("E:\\bladder\\vm.txt")
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())
load("E:\\bladder\\TCGA_tpm.Rdata")
#ssgsea????
ssgseaScore=gsva(count, geneSet, method='ssgsea', kcdf='Poisson', abs.ranking=TRUE)
vmscore=gsva(count, list(vm$V1), method='gsva', kcdf='Poisson', abs.ranking=TRUE)
ssgseaScore=rbind(ssgseaScore,vmscore)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
gsea=normalize(ssgseaScore)
imm = matrix(NA,nrow = 5,ncol = ncol(gsea))
rownames(imm) = imgene
for (i in imgene) {
  if(length(count[rownames(count) %in% i,]) != 0){
    imm[i,] = log2(count[i,]+1)
  }
}
rownames(imm) = c("CTLA-4","PD-L1","PD-1","TIGIT","PD-L2")
gsea = rbind(gsea[1:29,],imm,gsea[30,],log2(count[gene,]+1))
gsealist= list()
gsealist[["TCGA"]]=gsea

load("E:\\bladder\\GEO-BLCA.Rdata")
rname = rownames(count)
count = apply(count,2,as.numeric)
rownames(count) = rname
ssgseaScore=gsva(count, geneSet, method='ssgsea', kcdf='Poisson', abs.ranking=TRUE)
vmscore=gsva(count, list(vm$V1), method='gsva', kcdf='Poisson', abs.ranking=TRUE)
ssgseaScore=rbind(ssgseaScore,vmscore)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
gsea=normalize(ssgseaScore)
imm = matrix(NA,nrow = 5,ncol = ncol(gsea))
rownames(imm) = imgene
for (i in imgene) {
  if(length(count[rownames(count) %in% i,]) != 0){
    imm[i,] = log2(count[i,]+1)
  }
}
rownames(imm) = c("CTLA-4","PD-L1","PD-1","TIGIT","PD-L2")
gsea = rbind(gsea[1:29,],imm,gsea[30,],log2(count[gene,]+1))
gsealist[["GEO"]]=gsea

load("E:\\bladder\\im210.Rdata")
count = count210
ssgseaScore=gsva(count, geneSet, method='ssgsea', kcdf='Poisson', abs.ranking=TRUE)
vmscore=gsva(count, list(vm$V1), method='gsva', kcdf='Poisson', abs.ranking=TRUE)
ssgseaScore=rbind(ssgseaScore,vmscore)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
gsea=normalize(ssgseaScore)
imm = matrix(NA,nrow = 5,ncol = ncol(gsea))
rownames(imm) = imgene
for (i in imgene) {
  if(length(count[rownames(count) %in% i,]) != 0){
    imm[i,] = log2(count[i,]+1)
  }
}
rownames(imm) = c("CTLA-4","PD-L1","PD-1","TIGIT","PD-L2")
gsea = rbind(gsea[1:29,],imm,gsea[30,],log2(count[gene,]+1))
gsealist[["IM210"]]=gsea

clinic = list(
  TCGA = as.data.frame(clinicdata),
  GEO = as.data.frame(clinicgeo),
  IM210= clinic210
)


save(gsealist,clinic,file="ssGSEA.Rdata")

##############???ӻ?####################33
load("ssGSEA.Rdata")
TCGA = gsealist[["TCGA"]]
cli = clinic$TCGA
Normal_TCGA = TCGA[,!(cli$shortLetterCode %in% "TP" )]
Tumor_TCGA = TCGA[,cli$shortLetterCode %in% "TP"]
geo = gsealist[["GEO"]]
cli = clinic[["GEO"]]
Control_geo = geo[,cli$geotype %in% "Control"]
Surrounding_geo = geo[,cli$geotype %in% "Surrounding"]
Primary_geo = geo[,cli$geotype %in% "Primary"]
Recurrent_geo = geo[,cli$geotype %in% "Recurrent"]

im210 = gsealist[["IM210"]]
cli  = clinic[["IM210"]]
SDPD_im = im210[,cli$binaryResponse %in% "SD/PD"]
CRPR_im = im210[, cli$binaryResponse %in% "CR/PR"]
NE_im = im210[, cli$binaryResponse %in% NA]
gsea=list(
  Normal_TCGA = Normal_TCGA,
  Tumor_TCGA = Tumor_TCGA,
  Control_GEO = Control_geo,
  Surrounding_GEO=Surrounding_geo,
  Primary_GEO  = Primary_geo,
  Recurrent_geo = Recurrent_geo,
  SDPD_IM210 = SDPD_im,
  CRPR_IM210 = CRPR_im,
  NE_IM210 = NE_im
)
vmdata = matrix(NA,nrow = 34,ncol = 9)
rownames(vmdata) = rownames(TCGA[1:34,])
colnames(vmdata) = names(gsea)
genedata = matrix(NA,nrow = 34,ncol = 9)
rownames(genedata) = rownames(TCGA[1:34,])
colnames(genedata) = names(gsea)
for (i in names(gsea)) {
  part =  gsea[[i]]
  for (j in rownames(vmdata)){
    vmdata[j,i]  = cor(part[j,],part[35,])
    genedata[j,i]  = cor(part[j,],part[36,])
  }
}
vmdata=as.matrix(vmdata)
genedata = as.matrix(genedata)
data = cbind(vmdata,genedata)
library(pheatmap)
pdf(paste(gene,"checkpoint_cor_heatmap.pdf"),width =9,height = 12)
bk <- c(seq(-1,-0.01,by=0.02),seq(0,1,by=0.02))
pheatmap::pheatmap(data,na_col = "gray",cluster_rows = F,cluster_cols = F,main = "Angiogenesis & EFEMP1 Correlation with Immune checkpoints",
                   color = c(colorRampPalette(colors = c("#164194FF" ,"white"))(length(bk)/2),colorRampPalette(colors = c("white","#D51317FF" ))(length(bk)/2)),
                   ,display_numbers = TRUE,  #是否显示每个色块对应的数???(经归一化后的数???),
                   breaks = bk,
                   number_format = "%.2f",fontsize_col=12,angle_col = 315,gaps_col = 9,show_rownames = T,gaps_row = 29)
dev.off()
library(ggsci)
pal_frontiers()(9)


#################五个免疫基因#############
