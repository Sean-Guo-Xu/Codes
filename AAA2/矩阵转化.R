setwd("D:\\AAA\\scRNA")
load("7AAA.Rdata")

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
count = pbmc@assays$RNA@counts
memory.limit(size=10000000)
count = as_matrix(count)
colnames(count) = pbmc$cell
count=cbind(rownames(count),count)
colnames(count)[1] = "Gene symbol"
write.table(count,"scRNA counts.txt",row.names= F, col.names=T,quote = F,sep ="\t")