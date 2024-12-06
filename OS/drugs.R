library(pRRophetic)
library(GSVA)
setwd("E:\\bladder")
gene = "EFEMP1"
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
drugname <- unique( drugData2016$Drug.name)
vm = read.table("vm.txt")
filename= "TCGA_count"
load(paste(filename,".Rdata",sep=""))
vmscore=gsva(count, list(vm$V1), method='gsva', kcdf='Poisson', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore=normalize(vmscore)
exp = log2(count[gene,]+1)
name = NULL
for ( i in drugname) {
  tryCatch({
    predictedPtype <- pRRopheticPredict(testMatrix=count, 
                                        drug=i,
                                        tissueType = "urogenital_system", 
                                        batchCorrect = "eb",
                                        selection=1,
                                        dataset = "cgp2016")
    drug = rbind(drug,predictedPtype)
    name = c(name,i)
  }, error = function(e){
    predictedPtype = rep("NaN",ncol(count))
    drug = rbind(drug,predictedPtype)
    name = c(name,i)
  })
}
rownames(drug) = name
druglist = list()
druglist[[filename]] = drug
druglist[[filename]] = rbind(druglist[[filename]],vmscore,exp)

filename= "GEO-BLCA"
load(paste(filename,".Rdata",sep=""))
names = rownames(count)
count = apply(count,2,as.numeric)
rownames(count) = names
vmscore=gsva(count, list(vm$V1), method='gsva', kcdf='Poisson', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore=normalize(vmscore)
exp = log2(count[gene,]+1)
drug = NULL
for ( i in drugname) {
  tryCatch({
    predictedPtype <- pRRopheticPredict(testMatrix=count, 
                                        drug=i,
                                        tissueType = "urogenital_system", 
                                        batchCorrect = "eb",
                                        selection=1,
                                        dataset = "cgp2016")
    drug = rbind(drug,predictedPtype)
  }, error = function(e){
    predictedPtype = rep(0,ncol(count))
    drug = rbind(drug,predictedPtype)
  })
}
druglist[[filename]] = drug
rownames(druglist[[filename]]) = name
druglist[[filename]] = rbind(druglist[[filename]],vmscore,exp)

filename= "IM210"
load(paste(filename,".Rdata",sep=""))
count  = count210
vmscore=gsva(count, list(vm$V1), method='gsva', kcdf='Poisson', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore=normalize(vmscore)
exp = log2(count[gene,]+1)
drug = NULL

name = NULL
for ( i in drugname) {
  tryCatch({
    predictedPtype <- pRRopheticPredict(testMatrix=count, 
                                        drug=i,
                                        tissueType = "urogenital_system", 
                                        batchCorrect = "eb",
                                        selection=1,
                                        dataset = "cgp2016")
    drug = rbind(drug,predictedPtype)
    name = c(name,i)
  }, error = function(e){
    predictedPtype = rep("NaN",ncol(count))
    drug = rbind(drug,predictedPtype)
    name = c(name,i)
  })
}
rownames(drug) = name
druglist[[filename]] = drug
rownames(druglist[[filename]]) = name
druglist[[filename]] = rbind(druglist[[filename]],vmscore,exp)


clinic = list(
  TCGA = as.data.frame(clinicdata),
  GEO = as.data.frame(clinicgeo),
  IM210= clinic210
)
save(clinic,druglist,file="drug.Rdata")
###################
load("drug.Rdata")

TCGA = druglist[["TCGA_count"]]
cli = clinic$TCGA
Normal_TCGA = TCGA[,!(cli$shortLetterCode %in% "TP" )]
Tumor_TCGA = TCGA[,cli$shortLetterCode %in% "TP"]
geo = druglist[["GEO-BLCA"]]
cli = clinic[["GEO"]]
Control_geo = geo[,cli$geotype %in% "Control"]
Surrounding_geo = geo[,cli$geotype %in% "Surrounding"]
Primary_geo = geo[,cli$geotype %in% "Primary"]
Recurrent_geo = geo[,cli$geotype %in% "Recurrent"]

im210 = druglist[["IM210"]]
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
vmdata = matrix(NA,nrow =237,ncol = 9)
rownames(vmdata) = rownames(TCGA[1:237,])
colnames(vmdata) = names(gsea)
genedata = matrix(NA,nrow = 237,ncol = 9)
rownames(genedata) = rownames(TCGA[1:237,])
colnames(genedata) = names(gsea)
for (i in names(gsea)) {
  part =  gsea[[i]]
  for (j in rownames(vmdata)){
    x=part[j,]
    x = x[!is.na(x)]
    y1 = part[238,colnames(part) %in% names(x)]
    y2 = part[239,colnames(part) %in% names(x)]
    vmdata[j,i]  = cor(x,y1)
    genedata[j,i]  = cor(x,y2)
  }
}
vmdata=as.matrix(vmdata)
genedata = as.matrix(genedata)
data = cbind(vmdata,genedata)
data=t(data)
library(pheatmap)
pdf(paste(gene,"drug_heatmap.pdf"),width =18,height = 6)
bk <- c(seq(-1,-0.01,by=0.02),seq(0,1,by=0.02))

pheatmap::pheatmap(data,na_col = "gray",cluster_rows = F,cluster_cols = T,main = "Correlations with Sensitivity of 237 Drugs",
                   color = c(colorRampPalette(colors = c("#0094CDFF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#D51317FF" ))(length(bk)/2)),
                   ,display_numbers =F,  #是否显示每个色块对应的数???(经归一化后的数???),
                   breaks = bk,
                   number_format = "%.2f",fontsize_col=5,fontsize_row = 12,angle_col = 90,gaps_row = 9,show_rownames = T)
dev.off()
library(ggsci)
pal_frontiers()
