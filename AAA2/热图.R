library(circlize)
library(RColorBrewer)
library(ggsci)
library(ComplexHeatmap)
setwd("D:\\AAA\\multiomics")
load("multiomics_data.Rdata")
data=data[!duplicated(data$Name),]
data=data[order(data$EP300_AAA_Score),]
data=data[order(data$EP300_down_RNAseq),]
data=data[order(data$EP300_diff_chipseq),]
data=data[order(data$Pathway),]
creattext=function(text,col){circos.text(0:(length(data$Name)-1)+0.5, rep(1.1, length(data$Name)), text, cex = 0.5, 
            font = 0.1, 
            col = col, 
            facing = "clockwise", niceFacing = TRUE, 
            adj = c(0, 0.5))}
########蛋白组和功能准备#############
procol = c("orange","white","purple")
names(procol)= c("UP","NON","DOWN")  
phycol = pal_npg()(9)
udcol = c("#E64B35FF","white","#4DBBD5FF")
names(udcol) = c("UP","NON","DOWN")
genelink1 = data$Name
genelink1[!(data$Interact %in% "EP300_Interact")] = ""
genelink2 = data$Name
genelink2[!(data$Interact %in% "NON")] = ""
pathwaycol=c(phycol[3],phycol[4],"yellow3",phycol[6],phycol[7])
names(pathwaycol) = unique(data$Pathway)
###############chipseq##########
mat = data[,c(5:9)]
mat[mat > 0] = "1"
chipcol1 = c("white",phycol[9])
names(chipcol1) = c("0","1")
EP300chipcol = c("white","black")
names(EP300chipcol) = c("0","1")
chipcol2 = c("white","brown")
names(chipcol2) = c("0","1")
###########rnaseq###############
rnamat = data[,11:13]
rnacol = c("pink1","white","blue1")
names(rnacol) = c("UP","NON","DOWN")

rnacol2 = c("red","white","blue")
names(rnacol2) = c("UP","NON","DOWN")
##########scRNA##############
scmat1 = data[,21:27]
sccol1 = colorRamp2(c(min(scmat1[!is.na(scmat1)]),0,max(scmat1[!is.na(scmat1)])),c("blue2","white","red2"))
scmat2 = data[,28:34]
sccol2 = colorRamp2(c(min(scmat2[!is.na(scmat1)]),0,max(scmat2[!is.na(scmat1)])),c("purple3","white","yellow1"))
#############绘图###############
pdf("circos.pdf",width = 15,height = 15)

circos.heatmap(data[,2], col = pathwaycol, track.height = 0.02,track.margin=c(0.02,0),rownames.side = "outside",cluster = F)
creattext(genelink1,"red")
creattext(genelink2,"black")
circos.heatmap(data$Ratio, col = procol, track.height = 0.02,track.margin=c(0.04,0),cluster = F)
circos.heatmap(mat[,1], col = EP300chipcol, track.height = 0.01,track.margin=c(0.01,0),cluster = F)
circos.heatmap(mat[,2:3], col = chipcol1, track.height = 0.02,track.margin=c(0.01,0),cluster = F)
circos.heatmap(mat[,4:5], col = chipcol2, track.height = 0.02,track.margin=c(0.01,0),cluster = F)
circos.heatmap(data$EP300_diff_chipseq, col = udcol, track.height = 0.02,track.margin=c(0.04,0),cluster = F)
circos.heatmap(rnamat, col = rnacol, track.height = 0.03,track.margin=c(0.01,0),cluster = F)
circos.heatmap(data$EP300_down_RNAseq, col = rnacol2, track.height = 0.02,track.margin=c(0.04,0),cluster = F)
circos.heatmap(scmat1, col = sccol1, track.height = 0.07,track.margin=c(0.03,0),cluster = F)
circos.heatmap(scmat2, col = sccol2, track.height = 0.07,track.margin=c(0.01,0),cluster = F)

dev.off()
circos.clear()
##############圈图##############

nameleg=Legend(title = "GeneName", at = c("Interact with EP300"), 
              legend_gp = gpar(fill = "red"))
pathwayleg=Legend(title = "Pathways", at = names(pathwaycol), 
                legend_gp = gpar(fill = pathwaycol ))
proleg = Legend(title = "Proteome (AAA-Con)", at = names(procol[-2]), 
                  legend_gp = gpar(fill = procol[-2]))
chipleg = Legend(title="ChIPseq",at =c("EP300","H3K27ac-hsa","H3K27ac-mm"),
                  legend_gp = gpar(fill = c(phycol[9],"black","brown")))
diffchipleg= Legend(title="EP300 diffChIPseq (AAA-Con)",at =names(udcol[-2]),
                    legend_gp = gpar(fill = udcol[-2]))
gsernaleg=Legend(title="GSE RNAseq datasets (AAA-Con)",at =names(rnacol[-2]),
                 legend_gp = gpar(fill = rnacol[-2]))
ep300leg = Legend(title="RNAseq (EP300knockout-Con)",at =names(rnacol2[-2]),
                  legend_gp = gpar(fill = rnacol2[-2]))
sclogfcleg=Legend(title = "scRNAseq logFC (AAA-Con)",col_fun=sccol1)
scpctleg=Legend(title = "scRNAseq pct difference (AAA-Con)",col_fun=sccol2)

lgd_list = packLegend(nameleg,pathwayleg,proleg,chipleg,diffchipleg,gsernaleg,ep300leg,sclogfcleg,scpctleg,max_height = unit(0.6*dev.size()[2], "inch"))
pdf("legend.pdf",width = 8,height = 10)
draw(lgd_list)
dev.off()
