library(stringr)
setwd("D:\\生信\\target数据分析\\新mi-m")
mi_m=read.table("starBase_Human_Pan-Cancer_MiRNA-Target_Interactions2023-02-12_03-38.xls",header=T)
circ_mi=read.table("starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2023-02-12_18-07.xls",header=T)
circ=read.table("circdiff.txt",header=T)
circ=rownames(circ)
cut=strsplit(circ_mi$geneName,"_")
count=1
for (i in cut) {
    circ_mi[count,3]=i[2]
  count=count+1
}
circ_mi[,3]=str_sub(circ_mi[,3],12,-1)
circ_mi[,3]=str_pad(circ_mi[,3],7,side="left","0")
circ=substr(circ,10,16)
circ_mi=circ_mi[which(circ_mi$geneName %in% circ),]
demi=read.table("demi.xls",header=T)
circ_mi=circ_mi[which(circ_mi$name %in%demi$id),]
mi=c("hsa-miR-139-5p",
"hsa-miR-7-5p",
"hsa-miR-221-3p","
hsa-miR-758-3p")
diff=read.table("all.xls",header=T)
diff=diff[abs(diff$logFC)>2,]
mi_m=mi_m[which(mi_m$geneName %in% diff$gene),]
mi_m=mi_m[which(mi_m$name %in% mi),]
write.table(mi_m[,1:2],"edges.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(c(unique(mi_m$name),unique(mi_m$geneName)),"nodes.txt",sep = "\t",quote = F,row.names = F,col.names = F)
