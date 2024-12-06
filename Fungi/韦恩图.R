library(ggvenn)
name=read.table("D:\\bigslice\\gcf-rank\\Phylum.txt",header=T)
class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
gcf=list()
for (i in name$Name){
  e=class[which(class$V6 %in% i),]
  e=e[!duplicated(e$V4),]
  gcf[[i]]=e$V4
}
mat=c(rep(0,length(name$Name)))
for (i in 1:(length(name$Name)-1)){
  mat=cbind(mat,rep(0,length(name$Name)))}

colnames(mat)=c(1:nrow(mat))
rownames(mat)=c(1:nrow(mat))
mc=5
for (i in 1:8){
  mc=c(mc,gcf[[i+2]])
}
mc=mc[!duplicated(mc)]
ve=list("Ascomycota"= gcf$Ascomycota,"Basidiomycota"=gcf$Basidiomycota ,"Others"= mc)
p=ggvenn(ve,fill_color = c("cyan3","brown3","yellow3"))
ggsave("D:\\bigslice\\phylum venn.tiff",p,height = 6,width = 5)
