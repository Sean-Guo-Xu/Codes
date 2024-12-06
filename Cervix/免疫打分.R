setwd("E:\\cervix")     #???ù???Ŀ¼                                       #?????ļ?
gmtFile="./infla/immune.gmt"                                           #GMT?ļ?
gene = "HDAC10"
imgene = c("CTLA4","CD274","PDCD1","TIGIT","PDCD1LG2")
library(GSEABase)
library(GSVA)
library(limma)
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())
load("combined.Rdata")
 #ssgsea????
count=rt
ssgseaScore=gsva(count, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
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
gsea = rbind(gsea[1:29,],imm,log2(count[gene,]+1))

gsea=list(
  Normal = gsea[, 1:91 ],
  Tumor = gsea[,92:397]
)
genedata = matrix(NA,nrow = 34,ncol = 2)
rownames(genedata) = rownames(gsea[[1]][1:34,])
colnames(genedata) = names(gsea)
for (i in names(gsea)) {
  part =  gsea[[i]]
  for (j in rownames(genedata)){
    genedata[j,i]  = cor(part[j,],part[35,])
  }
}
genedata = as.matrix(genedata)
data = genedata
library(pheatmap)
pdf(paste(gene,"checkpoint_cor_heatmap.pdf"),width =3,height = 12)
bk <- c(seq(-1,-0.01,by=0.02),seq(0,1,by=0.02))
pheatmap::pheatmap(data,na_col = "gray",cluster_rows = F,cluster_cols = F,main = gene,
                   color = c(colorRampPalette(colors = c("#164194FF" ,"white"))(length(bk)/2),colorRampPalette(colors = c("white","#D51317FF" ))(length(bk)/2)),
                   ,display_numbers = TRUE,  #是否显示每个色块对应的数???(经归一化后的数???),
                   breaks = bk,
                   number_format = "%.2f",fontsize_col=12,angle_col = 315,show_rownames = T,gaps_row = 29)
dev.off()
library(ggsci)
pal_frontiers()(9)


#################五个免疫基因#############
