library(ape)
library(itol.toolkit)
library(data.table)
library(pheatmap)
setwd("D:\\bigslice\\pictures")
load("D:\\bigslice\\GCF cluster\\seurat gcfmatrix.Rdata")
see=pbmc@reductions$pca@cell.embeddings
heat = pheatmap(as.matrix(see))
tree = as.phylo(heat$tree_row)
write.tree(tree,"phylotree")
class = read.table("D:\\bigslice\\bgc_withclass.txt",sep="\t")
class[class$V10 %in% "Leiotrametes",10] =  "Cubamyces"
class[class$V10 %in% "Tinctoporellus",10] = "Porogramme"
class[class$V10 %in% "Taxomyces",10] = "Ceriporiopsis"
class[class$V10 %in% "Claroideoglomus",10] = "Entrophospora"
class[class$V10 %in% "Spiromastix",10] = "Spiromastigoides"
class[class$V10 %in% "Ochroconis",10] = "Scolecobasidium"
allclass = class
class= read.table("D:\\bigslice\\11608.txt",sep="\t")
for(i in 1:9){
  class[,i]  =gsub(" ","",class[,i])
}
class = class[!duplicated(class$V7),]
class[class$V7 %in% "Leiotrametes",7] =  "Cubamyces"
class[class$V7 %in% "Tinctoporellus",7] = "Porogramme"
class[class$V7 %in% "Taxomyces",7] = "Ceriporiopsis"
class[class$V7 %in% "Claroideoglomus",7] = "Entrophospora"
class[class$V7 %in% "Spiromastix",7] = "Spiromastigoides"
class[class$V7 %in% "Ochroconis",7] = "Scolecobasidium"
class[class$V7 %in% "Ceriporiopsis",3] = "Basidiomycota"
tree <- system.file("extdata",
                    "tree_of_itol_templates.tree",
                    package = "itol.toolkit")

tree_1 <- "https://raw.githubusercontent.com/TongZhou2017/itol.toolkit/master/inst/extdata/dataset3/abunt-tree.nwk"
hub_1 <- create_hub(tree_1)
data_file_1 <- "https://raw.githubusercontent.com/TongZhou2017/itol.toolkit/master/inst/extdata/dataset3/abunt-metadata.txt"
data_1 <- data.table::fread(data_file_1)
unit_2 <- create_unit(data = data, 
                      key = "rep_Zheng2022ep_3al_2_range", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      tree = tree_1)

tree ="D:\\bigslice\\pictures\\genus.nwk"
hub <- create_hub(tree = tree)
GCF = hub@meta.data$tip
class = read.table("D:\\bigslice\\bgc_withclass.txt",sep="\t")
class$V4 = paste(class$V4,"_f",sep = "")
class = class[class$V4 %in% GCF[,1],]
data_1 = class
unit_1 <- create_unit(data = data_1[,c(4,13)], 
                      key = "rep_Zheng2022ep_3al_1_labels", 
                      type = "LABELS",
                      tree = tree_1)
unit_1 <- create_unit(data = data_1[ ,c(4,13)], 
                      key = "rep_Zheng2022ep_3al_2_range", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      tree = tree)
hub_1 <- hub+unit_1
write_hub(hub_1,getwd())
