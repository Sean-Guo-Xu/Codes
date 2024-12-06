library(cellcall)
library(Seurat)
library(ggplot2)
setwd("D:\\AAA\\scRNA")             #设置工作目录
load("Normal.Rdata")

testAB.integrated=subset(testAB.integrated,sample=="Normal")
testAB.integrated=subset(testAB.integrated,cell!=" ")
readCounts=(testAB.integrated@assays$RNA@counts)



for (i in 1:length(colnames(readCounts))) {
  colnames(readCounts)[i] = paste(i,"_",gsub("-", " ",testAB.integrated$cell[i]),sep = "")
}

rm(testAB.integrated)
memory.limit(size = 1000000)
readCounts<-as.data.frame(readCounts)
mt <- CreateNichConObject(data=readCounts, min.feature = 0,
                          names.field = 2,
                          names.delim = "_",
                          source = "fullLength",
                          scale.factor = 10^6,
                          Org = "Homo sapiens",
                          project = "Microenvironment")
rm(readCounts)
gc()
mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')
n <- mt@data$expr_l_r_log2_scale
pathway.hyper.list <- lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
myPub.df=myPub.df[which(myPub.df$pvalue<0.05),]
myPub.df=myPub.df[which(abs(myPub.df$NES)>0.1),]
my=myPub.df
p=ggplot(my,aes(x=my$cc,y=my$pathway,color=my$NES))+geom_point(aes(size=10))+theme_bw()+theme(legend.box = "horizontal",axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid()+xlab("")+ylab("")+scale_color_gradient(low="cyan3",high="brown3")
ggsave("D:\\callNES.tiff",p,width = 18,height = 3)

pdf(file="AMLcheck.pdf",width=15,height=5)
plotBubble(myPub.df)
dev.off()
library(RColorBrewer)
library(ggsci)
cell_color <- data.frame(color=colorRampPalette((pal_npg( "nrc")(9)))(13), stringsAsFactors = FALSE)
rownames(cell_color) <- unique(mt@meta.data$cell)
cell_color <- data.frame(color=c("#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F"), stringsAsFactors = FALSE)
rownames(cell_color) <- c("Plasma","CD8","B","pDC","NK","CD4")
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 8,angle_col = "45",  
             main="score")
pdf(file="Norcheckcircle.pdf",width=10,height=10)

ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 0.01,
                
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)
dev.off()
save(mt,file="Normal.Rdata")
