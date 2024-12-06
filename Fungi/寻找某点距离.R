setwd("D:\\bigslice")
load("seurat gcfmatrix.Rdata")
see=pbmc@reductions$pca@cell.embeddings
re=NULL
id="7468_f"
for (i in rownames(see)) {
part=see[rownames(see) %in% c(id,i),]
alculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
if (i !=id){
dis=alculateEuclideanDistance(c(part[rownames(part)==id,]),c(part[rownames(part)==i,]))
re=rbind(re,c(i,dis))}
}
re=as.data.frame(re)
re[,2]=as.numeric(re[,2])
re=re[order(re$V2),]
re=cbind(re,1:26824)
write.table(re,paste(id,"distance.txt"),sep = "\t",quote=F,row.names=F,col.names = F)