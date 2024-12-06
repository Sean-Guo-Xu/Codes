setwd("E:\\cervix")
library(tidyr)
library(immunedeconv)
load("TCGA_infiltration.Rdata")
##########????list########
cell=c("B cell","NK cell","T cell","Basophil","Eosinophil","Mast cell","Macrophage/Monocyte","Neutrophil","Myeloid dendritic cell")
inf_list[["estamite"]] = NULL
rownames(inf_list[["cibersort"]])=NULL
ciber = inf_list[["cibersort"]]
type=ciber$cell_type
ciber = ciber[ ,2:ncol(ciber)]

cibersort<- ciber[1,] +  ciber[2,]+ciber[3,]
cibersort<- rbind(cibersort,c(ciber[11,]+ciber[12,]))
cibersort<-rbind(cibersort,c(ciber[4,]+ciber[5,]+ciber[6,]+ciber[7,]+ciber[8,]+ciber[9,]+ciber[10,]))
cibersort<-rbind(cibersort,rep(NA,nrow(ciber)))
cibersort<-rbind(cibersort,ciber[21,])
cibersort<-rbind(cibersort,ciber[19,]+ciber[20,])
cibersort<-rbind(cibersort,ciber[13,]+ciber[14,]+ciber[15,]+ciber[16,])
cibersort<-rbind(cibersort,ciber[22,])
cibersort<-rbind(cibersort,ciber[17,]+ciber[18,])
rownames(cibersort) = cell
cibersort=as.data.frame(cibersort)
for (i  in names(inf_list)) {
  if(i != "cibersort"){
inf_list[[i]] <- map_result_to_celltypes(inf_list[[i]],cell,i)}
}
inf_list[["cibersort"]] = cibersort
###########################
count = read.table("TCGA_tpm.txt" , sep="\t")
library(GSVA)

count = as.matrix(count)
count = log2(count+1)

inf_list[["estimate"]]=NULL

gene = "HDAC10"
pm = matrix(NA,nrow = 9,ncol = 8)
allvm = matrix(NA,nrow = 9,ncol = 8)
rownames(pm) = rownames(inf_list[[1]])
colnames(pm) = names(inf_list)
rownames(allvm) = rownames(inf_list[[1]])
colnames(allvm) = names(inf_list)
for (i in names(inf_list)){
  inf= inf_list[[i]]
  for (j in rownames(inf)) {
   
    allvm[j,i] = cor(as.numeric(inf[j,]),count[gene,])
  if(!is.na(inf[j,1]))  {
   p= cor.test(as.numeric(inf[j,]),count[gene,])
   p = p$p.value
   pm[j,i] = p}
   }
}

feature = pm
feature[pm<0.001] = "***" 
feature[0.001<pm & pm<0.01] = "**"
feature[0.01<pm & pm<0.05] = "*"
feature[!(feature %in% c("*","***","**",NA))]= ""

feature[feature %in% NA]=""
data=allvm
data[!(feature %in% c("*","**","***"))]=NA
data = round(data,digits = 2)
feature[feature %in% c("*","**","***")]=paste(feature,"\n",data,sep = "")[feature %in% c("*","**","***")]
library(ggsci)
pal_frontiers()(10)

data = cbind(allvm)
data= apply(data, 2, as.numeric)
rownames(data) = rownames(allvm)
library(pheatmap)
pdf(paste(gene,"IMcell_cor_heatmap.pdf"),width =5.5,height = 5)
bk <- c(seq(-1,-0.02,by=0.02),seq(0,1,by=0.02))

pheatmap(data,cluster_row = F, cluster_col = F, border=NA,main = paste(gene,"Correlation with Infiltration"),
         display_numbers = feature,number_color = "black",na_col = "gray",angle_col = 45,
         color = colorRampPalette(c("#51B1B7", "#E1C855" ,"#E07B54"))(50),
       scale = "none",
)

dev.off()
library(ggplot2)
cellcol = pal_frontiers()(9)
colnames(count) = gsub("[.]","-",colnames(count))
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
      part = rbind(part,c(j,cell,inf[cell,j],count[gene,j]))
    }
  }
  part = as.data.frame(part)
  colnames(part) = c("Samples","Cells","Proportion","Expression")
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


