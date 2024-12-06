library(circlize)
library(RColorBrewer)
library(ggsci)
library(ComplexHeatmap)
setwd("D:\\bigslice\\pictures")
data=read.table("D:\\bigslice\\pictures\\circos data.txt",sep="\t",header = T)
data=data[order(data$BGCNum,decreasing = T),]
phycol = pal_npg()(9)
phycol = c(phycol,"gray")
names(phycol)=c("Ascomycota","Basidiomycota","Mucoromycota","Chytridiomycota",
                "Zoopagomycota","Blastocladiomycota","Olpidiomycota","Microsporidia","Cryptomycota","Others")

classcol = pal_nejm()(7)
names(classcol) = colnames(data)[23:29]

knowncol = c("cyan3","brown3")
names(knowncol)=c("1","0")
data$Known.Unknown = as.character(data$Known.Unknown)

groupcol<-colorRampPalette(brewer.pal(8,'Set2'))(36)
names(groupcol)=as.character(1:36)

data$BGCNum = log10(data$BGCNum+1)
data$GenomesNum = log10(data$GenomesNum+1)
data$StrainNum = log10(data$StrainNum+1)
data$SpeciesNum = log2(data$SpeciesNum+1)
data$GenusNum = log2(data$GenusNum+1)
data$FamilyNum = log2(data$FamilyNum+1)
data$OrderNum = log2(data$OrderNum+1)
numcol= pal_lancet()(9)
names(numcol) = c("SpeciesNum","PhylumNum","FamilyNum","GenusNum","ClassNum","OrderNum","StrainNum","GenomesNum","BGCNum")  
circleorder = c("PhylumNum","ClassNum","OrderNum","FamilyNum","GenusNum","SpeciesNum","StrainNum","GenomesNum","BGCNum")
creatlabel =  function(){
   circos.track(
     track.index = get.current.track.index(), 
     panel.fun = function(x, y) {
       circos.rect(
         CELL_META$cell.xlim[2]+1 , 0,
         CELL_META$cell.xlim[2] +20, CELL_META$ylim[2],
         col = "green", border = NA
       )
       
     }, bg.border = NA
   )
}
creatwhite = function(){
  white = "white"
  names(white) = "empty"
  circos.heatmap(rep("empty",26825), col = white, track.height = 0.05,track.margin=c(0,0))
}

pdf("circos.pdf",width = 10,height = 10)
circos.heatmap.initialize(data,cluster = F)
for(i in names(phycol)){
  col = colorRamp2(c(0,1),c("white",phycol[i]))
  circos.heatmap(data[,i], col = col, track.height = 0.01,track.margin=c(0,0.005))
  }


creatwhite()
for(i in names(classcol)){
  col = colorRamp2(c(0,1),c("white",classcol[i]))
  circos.heatmap(data[,i], col = col, track.height = 0.01,track.margin=c(0,0.005))
}
creatwhite()

for (i in circleorder) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  circos.heatmap(data[,i], col = col, track.height = 0.02,track.margin=c(0,0.005))
}
dev.off()
circos.clear()
names(knowncol)=c("Known","Unknown")
phyleg=Legend(title = "Pylum Proportion", at = names(phycol), 
       legend_gp = gpar(fill = phycol))
groupleg=Legend(title = "GCF Group", at = names(groupcol), 
              legend_gp = gpar(fill = groupcol ),ncol=4)
classleg = Legend(title = "Type Proportion", at = names(classcol), 
                  legend_gp = gpar(fill = classcol))
knownleg = Legend(title="Known&Unknown",at =names(knowncol),
                  legend_gp = gpar(fill = knowncol))
listleg = list()

for (i in circleorder[1:2]) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  listleg[[i]] =Legend(title = i, col_fun = col )
}
for (i in circleorder[3:6]) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  listleg[[i]] =Legend(title = paste("log2(",i,")",sep = ""), col_fun = col )
}
for (i in circleorder[7:9]) {
  col = colorRamp2(c(0,max(data[,i])),c("white",numcol[i]))
  listleg[[i]] =Legend(title = paste("log10(",i,")",sep = ""), col_fun = col )
}
h = dev.size()[2]
lgd_list = packLegend(phyleg,groupleg,classleg,knownleg,listleg[[1]],listleg[[2]],listleg[[3]]
                      ,listleg[[4]],listleg[[5]],listleg[[6]],listleg[[7]],listleg[[8]],listleg[[9]],max_height = unit(0.9*h, "inch"))
pdf("legend.pdf",width = 8,height = 5)
draw(lgd_list)
dev.off()
