setwd("E:\\bladder\\infiltration")
load("9_infiltration.Rdata")
library(GSVA)
gene = "SVIL"
vm = read.table("E:\\bladder\\vm.txt")
see=gsva(expr = count,gset.idx.list = list(vm$V1),kcdf="Gaussian", verbose=T)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

see=normalize(see)
genem = matrix(NA,nrow = 9,ncol = 8)
allvm = matrix(NA,nrow = 9,ncol = 8)
rownames(genem) = rownames(inf_list[[1]])
colnames(genem) = names(inf_list)
rownames(allvm) = rownames(inf_list[[1]])
colnames(allvm) = names(inf_list)
for (i in names(inf_list)){
  inf= inf_list[[i]]
  for (j in rownames(inf)) {
    allvm[j,i] = cor(as.numeric(inf[j,]),see[1,])
    genem[j,i] = cor(as.numeric(inf[j,]),count[gene,])
  }
}
library(ggsci)
pal_frontiers()(10)

data = cbind(allvm,genem)
data= apply(data, 2, as.numeric)
rownames(data) = rownames(allvm)
library(pheatmap)
pdf(paste(gene,"IMcell_cor_heatmap.pdf"),width =9,height = 5)
bk <- c(seq(-1,-0.01,by=0.02),seq(0,1,by=0.02))
pheatmap::pheatmap(data,na_col = "gray",cluster_rows = F,cluster_cols = F,main = paste("Angiogenesis &",gene,"Correlation with Infiltration"),
                   color = c(colorRampPalette(colors = c("#4DBBD5FF" ,"white"))(length(bk)/2),colorRampPalette(colors = c("white","#E64B35FF"))(length(bk)/2)),
                   ,display_numbers = TRUE,  #是否显示每个色块对应的数�?(经归一化后的数�?),
                   breaks = bk,
                   number_format = "%.2f",fontsize_col=12,angle_col = 315,gaps_col = 8,show_rownames = T)
dev.off()
library(ggplot2)
cellcol = pal_frontiers()(9)
names(cellcol) = rownames(inf_list[[1]])
for (i in names(inf_list)) {
  inf=inf_list[[i]]
  inf = inf[!(is.na(inf[,1])),]
 # Others = 1-apply(inf,2,sum)
 # inf = rbind(inf,Others)
#  rownames(inf)[5] = "Others"
  inf = inf[,order(count[gene,],decreasing = T)]
  part = NULL
  for (j in colnames(inf)) {
    for (cell in rownames(inf)) {
      part = rbind(part,c(j,cell,inf[cell,j],count[gene,j],see[1,j]))
    }
  }
  part = as.data.frame(part)
  colnames(part) = c("Samples","Cells","Proportion","Expression","GSVA")
  part$Samples = factor(part$Samples,levels = unique(part$Samples))
  part$Proportion=as.numeric(part$Proportion)
  part$Expression=as.numeric(part$Expression)
 ggplot( part, aes( x = Samples, weight=Proportion, fill = Cells))+scale_fill_manual(values = cellcol)+
    geom_bar( position = "stack",width = 1)+coord_flip()+theme_classic()+theme(panel.border = element_rect(colour = "black",fill=NA),axis.text.y = element_blank(),axis.ticks.y=element_blank())+scale_y_continuous(expand = c(0,0),limits = c(0,ceiling(max(apply(inf, 2, sum))*10)/10))+
    xlab("")+ylab("")+labs(title = i)
  ggsave(paste(i,"stack.tiff"),width = 4.5,height = 8)
}
part$GSVA =as.numeric(part$GSVA)
ggplot( part, aes( x = Samples, y=GSVA,fill=""))+scale_fill_manual(values =  "#F39200FF")+
  geom_bar( position="dodge", stat="identity",width = 1)+coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = 'none')+
  scale_y_continuous(expand = c(0,0),limits = c(0,ceiling(max(part$GSVA)*10)/10))+
  xlab("")+ylab("")+labs(title = paste("Angiogenesis","Score"))
ggsave("Gsva bar.tiff",width = 2,height = 8)
part$Expression =as.numeric(part$Expression)
ggplot( part, aes( x = Samples, y=Expression,fill=""))+scale_fill_manual(values =  "#D51317FF" )+
  geom_bar( position="dodge", stat="identity",width = 1)+coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = 'none')+
  scale_y_continuous(expand = c(0,0),limits = c(0,ceiling(max(part$Expression)*10)/10))+
  xlab("")+ylab("")+labs(title = paste(gene,"Expression"))
ggsave(paste(gene," bar.tiff"),width = 2,height = 8) 


