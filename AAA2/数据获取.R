setwd("D:\\AAA\\multiomics")

allgene=read.table("D:\\AAA\\multiomics\\自噬genelist.txt")
allgene = as.data.frame(allgene)
allgene = c("UGCG","PLEK")
allgene = as.data.frame(allgene)
colnames(allgene) = "Name"
#########功能################


##################蛋白##################
pro=read.table("protein_omics.txt",header = T)
data=merge(allgene,pro,by="Name",all.x = T)
data$Ratio[data$Ratio < 1] = 0
data$Ratio[data$Ratio > 1] = 1
data$Ratio = as.character(data$Ratio)
data$Ratio[data$Ratio == "1"] = "UP"
data$Ratio[data$Ratio == "0"] = "DOWN"
data$Ratio[is.na(data$Ratio)] = "NON"

string=read.table("string_interactions.tsv")
string = string[string$V2 %in% data$Name,] 
p300 = string[string$V1 %in% "EP300",]
p300 = p300[,c(1,2)]
colnames(p300)[2]="Name"
data = merge(data,p300,by="Name",all.x = T)
data[is.na(data$V1),4] = "NON"
data[data$V1 %in% "EP300",4] = "EP300_Interact"
colnames(data)[4]="Interact"
data= data[,-3]
#################ChIPseq##########
setwd(".\\chipseq")
namelist=c("EP300AAA","H3K27ac_hsa_1","H3K27ac_hsa_2","H3K27ac_mmu_1","H3K27ac_mmu_2")
for (i in namelist) {
  chip=read.table(paste(i,"chipseq_omics.txt"),header=T)
  chip = chip[!duplicated(chip$Name),]
  data=merge(data,chip,by="Name",all.x = T)
  data$Score[is.na(data$Score)] = 0
  colnames(data)[colnames(data) %in% "Score"]=paste(i,"Score",sep = "_")
}

diff=read.table("diffEP300.txt",sep = "\t",header = T)
diff=diff[diff$FDR <0.05,]
diff=diff[abs(diff$logFC) >2,]
diff=diff[!duplicated(diff$GeneID),]
colnames(diff)[1]="Name"
diff=diff[,c(1,10)]
data=merge(data,diff,by="Name",all.x = T)
data$logFC[is.na(data$logFC)] ="NON"
data$logFC[grep("-",data$logFC)]="DOWN"
data$logFC[!(data$logFC %in% c("NON","DOWN"))]="UP"
colnames(data)[colnames(data) %in% "logFC"]="EP300_diff_chipseq"
colnames(data)[colnames(data) %in% "EP300AAA_Score"]="EP300_AAA_Score"
###############ATAC################
setwd("D:\\AAA\\multiomics\\ATACseq")
namelist=c("art_78320_peaks","art_78321_peaks","art_78323_peaks","art_78324_peaks","art_78322_peaks","art_78370_peaks"
           ,"pb_103834_peaks","pb_103835_peaks","pb_103836_peaks","pb_103837_peaks")
for (i in namelist) {
  chip=read.table(paste(i,"chipseq_omics.txt"),header=T)
  data=merge(data,chip,by="Name",all.x = T)
  data$Score[is.na(data$Score)] = 0
  colnames(data)[colnames(data) %in% "Score"]=paste(i,"Score",sep = "_")
}

##############RNAseq###############
setwd("D:\\AAA\\multiomics\\rnaseq")
namelist=c("GSE7084","GSE57691","GSE47472")
for (i in namelist) {
  rna=read.table(paste(i,".txt",sep = ""),header=T)
  rna=rna[rna$adj.P.Val<0.05,]
  colnames(rna)[1] = "Name"
  rna=rna[,c(1,2)]
  rna=rna[abs(rna$logFC)>0.5,]
  rna = rna[!duplicated(rna$Name),]
  data=merge(data,rna,by="Name",all.x = T)
  data$logFC[is.na(data$logFC)] ="NON"
  data$logFC[grep("-",data$logFC)]="DOWN"
  data$logFC[!(data$logFC %in% c("NON","DOWN"))]="UP"
  colnames(data)[colnames(data) %in% "logFC"]=i
}

ep300=read.table("EP300.txt",header = T)
ep300=ep300[,1:7]
ep300=ep300[!duplicated(ep300$Name),]
colnames(ep300)[1]="Name"
ep300=read.table("EP300.txt",header = T)
colnames(ep300)[1]="Name"
ep300=ep300[ep300$FDR<0.05,]
ep300=ep300[abs(ep300$logFC)>0.5,]
ep300=ep300[,c(1,8)]
data=merge(data,ep300,by="Name",all.x = T)
data$logFC[is.na(data$logFC)] ="NON"
data$logFC[grep("-",data$logFC)]="DOWN"
data$logFC[!(data$logFC %in% c("NON","DOWN"))]="UP"
colnames(data)[colnames(data) %in% "logFC"]="EP300_down_RNAseq"

###########scRNA############
setwd("D:\\AAA\\multiomics")
load("artery_diff_all.Rdata")
markerlist=diff
for (i in names(markerlist)) {

  marker = markerlist[[i]]
  marker = marker[,-6]
  marker$Name = rownames(marker)
  marker = marker[,c(2,6)]
  data=merge(data,marker,by="Name",all.x = T)
  colnames(data)[colnames(data) %in% "avg_log2FC"]=paste("artery",i,"logFC",sep = "_")
}
for (i in names(markerlist)) {
  marker = markerlist[[i]]
  marker = marker[,-1]
  marker$Name = rownames(marker)
  marker$ratio = marker$pct.1-marker$pct.2
  marker = marker[,c(6,7)]
  data=merge(data,marker,by="Name",all.x = T)
  colnames(data)[colnames(data) %in% "ratio"]=paste("artery",i,"pct_df",sep = "_")
}

###########scRNA############
setwd("D:\\AAA\\multiomics")
load("blood_diff_all.Rdata")
markerlist=diff
for (i in names(markerlist)) {
  marker = markerlist[[i]]
  marker = marker[,-6]
  marker$Name = rownames(marker)
  marker = marker[,c(2,6)]
  data=merge(data,marker,by="Name",all.x = T)
  colnames(data)[colnames(data) %in% "avg_log2FC"]=paste("blood",i,"logFC",sep = "_")
}
for (i in names(markerlist)) {
  marker = markerlist[[i]]
  marker = marker[,-1]
  marker$Name = rownames(marker)
  marker$ratio = marker$pct.1-marker$pct.2
  marker = marker[,c(6,7)]
  data=merge(data,marker,by="Name",all.x = T)
  colnames(data)[colnames(data) %in% "ratio"]=paste("blood",i,"pct_df",sep = "_")
}
write.table(data,"自噬_dashboard.txt",row.names = F,col.names = T,quote = F,sep = "\t")
save(data,file="multiomics_data.Rdata")

