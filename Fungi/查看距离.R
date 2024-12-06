setwd("D:\\bigslice\\pictures")
system("R CMD SHLIB caculate.c")
dyn.load('D:/bigslice/pictures/caculate.dll')
load("D:\\bigslice\\GCF cluster\\seurat gcfmatrix.Rdata")
allout = list()
for(group in 1:36){
  
subpbmc = subset(pbmc, seurat_clusters == group)

see=subpbmc@reductions$pca@cell.embeddings
mat = matrix(data=1000,nrow = nrow(see),ncol =nrow(see))
colnames(mat) = gsub("_f","",rownames(see))
rownames(mat) = gsub("_f","",rownames(see))
see = as.vector(see)

c = see
re = .C("addme", as.double(c),as.double(c), as.integer(length(c)/20),as.double(c(rep(0,sum(c(1:(length(c)/20-1)))))))
re = re[[4]]
l=1
for(i in 1:(length(c)/20-1))
{
 
  mat[i+1,]=c(re[c(l:(l+i-1))],rep(1000,length(c)/20-i))
  l = l+i
}
out  = matrix(data=0,nrow = length(c)/20,ncol = 3)
for (i in 1:length(c)/20) {
  min = min(c(mat[,i],mat[i,]))
  node1 =rownames(mat)[i]
  if(min(mat[,i])<min(mat[i,]))
  { 
    node2 = rownames(mat)[mat[,i] %in% min]
    node2 = unique(node2)
  } else {  
  node2 = colnames(mat)[mat[i,] %in% min]
  node2 = unique(node2)
  }

out[i,]=c( node1,node2 ,min)
}
allout[[group]] = out
}
save(allout,file="GCF_O_distance.Rdata")
dyn.unload('D:/bigslice/pictures/caculate.dll')