library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(paletteer)
library(ggplot2)
setwd("D:\\bigslice\\pictures")
class = read.table("D:\\bigslice\\bgc_withclass.txt",sep="\t")
class = class[!(class$V6 %in% "Others"),]
dt=table(class$V6,class$V14)
dt = dt[order(rownames(dt)),]
dt2 = log2(dt+1)
col_fun = colorRamp2(c(0,4,16), c( "cyan4","white", "brown4"))
mypal1 <- pal_npg("nrc", alpha =1)(9)
mypal2 = pal_nejm( alpha =1)(8)


####################表#############
library(ggplot2)
all=read.table("D:\\bigslice\\new11608.txt",sep="\t")
withclass = read.table("D:\\bigslice\\化合物\\bgc_withclass.txt",sep="\t",header=F,fill=T)
genome=unique(withclass$V1)
num=c()
for (i in genome){
  num=c(num,nrow(withclass[which(withclass$V1 %in% i),]))
}
bgc=cbind(genome,num)
allbgc=all[-which(all$V1 %in% genome),]
nobgc=cbind(allbgc$V1,rep(0,length(allbgc$V1)))
allbgc=rbind(bgc,nobgc)
class=unique(all$V3)
data=c("class",1)
for (i in class)
{
  n=all[which(all$V3 %in% i),1]
  for (j in n){
    data=rbind(data,c(i,allbgc[which(allbgc[,1] %in% j),2]))
  }
}
data=data[-1,]
data=as.data.frame(data)
data$num = as.numeric(data$num)
bgc_genome = NULL
for(i in unique(data$V1)){
  part = data[data$V1 %in% i,]
  bgc_genome = c(bgc_genome,sum(part$num)/nrow(part))
}
bgc_genome = cbind(unique(data$V1),bgc_genome)

gcf = NULL
for(i in unique(withclass$V1)){
  part = withclass[withclass$V1 %in% i,]
  part = part[!duplicated(part$V4),]
  gcf = rbind(gcf,c(part[1,1],nrow(part)))
    
}

gcf=as.data.frame(gcf)
gcf$V2 =as.numeric(gcf$V2)
gcf_genome=NULL
for(i in unique(withclass$V6)){
  part1 = all[all$V3 %in% i,]
  part2 = gcf[gcf$V1 %in% part1$V1,]
  gcf_genome = rbind(gcf_genome,c(i,sum(part2$V2)/nrow(part1)))
}

gcf_genome=as.data.frame(gcf_genome)
gcf_genome[,2] = as.numeric(gcf_genome[,2])
gcf_genome = gcf_genome[!(gcf_genome$V1 %in% "Others"),]
gcf_bgc_genome = merge(gcf_genome,bgc_genome,by="V1")


gcfnum = NULL 
for(i in unique(withclass$V6)){
  part = withclass[withclass$V6 %in% i,]
  if(any(duplicated(part$V4))){
  part = part[!duplicated(part$V4),]}
  gcfnum = rbind(gcfnum,c(i,nrow(part)))
}
gcfnum = as.data.frame(gcfnum)
gcfnum$V2 = as.numeric(gcfnum$V2)
gcfnum = gcfnum[!(gcfnum$V1 %in% "Others"),]
gcf_bgc_genome = merge(gcfnum,gcf_bgc_genome,by="V1")
gcf_bgc_genome$bgc_genome = as.numeric(gcf_bgc_genome$bgc_genome)
bar1 <- HeatmapAnnotation(
  sum1 = anno_barplot(
    colSums(dt),
    bar_width = 0.9,
    gp = gpar(col = "white", fill = mypal2),
    border = F,
    axis_param = list(at = c(0,50000,100000),
                      labels = c("","50K","100k")
    ),
    height = unit(2, "cm")), show_annotation_name = F)
bar2 = rowAnnotation(
  sum2 = anno_barplot(
    log(rowSums(dt)+1,10),
    bar_width = 0.9,
    gp = gpar(col = "white", fill = mypal1),
    border = T,
    axis_param = list(at = c(0,1,2,3,4,5),
                      labels = c("0","10","100","1k","10k","100k"),
                      labels_rot=45
    ),
    
    width = unit(3, "cm")), show_annotation_name = F)

bar3 = rowAnnotation(
  sum3 = anno_barplot(
    log(gcf_bgc_genome[,2],10),
    bar_width = 0.9,
    gp = gpar(col = "white", fill = mypal1),
    border = T,
    axis_param = list(at = c(0,1,2,3,4),
                      labels = c("0","10","100","1k","10k"),
                      labels_rot=45
    ),
    
    width = unit(3, "cm")), show_annotation_name = F)

bar4 = rowAnnotation(
  sum4 = anno_barplot(
    gcf_bgc_genome[,4],
    bar_width = 0.9,
    gp = gpar(col = "white", fill = mypal1),
    border = T,
    axis_param = list(at = c(0,10,20),
                      labels = c("0","10","20"),
                      labels_rot=45
    ),
    
    width = unit(2, "cm")), show_annotation_name = F)
bar5 = rowAnnotation(
  sum5 = anno_barplot(
    gcf_bgc_genome[,3],
    bar_width = 0.9,
    gp = gpar(col = "white", fill = mypal1),
    border = T,
    axis_param = list(at = c(0,10,20),
                      labels = c("0","10","20"),
                      labels_rot=45
    ),
    
    width = unit(2, "cm")), show_annotation_name = F)
pdf("heatmap.pdf",width = 8.5,height = 5)
Heatmap(dt2,col = col_fun,
        name = "Log2BGCs",
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        top_annotation = bar1,
        right_annotation =bar2,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 45
)+bar3+bar4+bar5
dev.off()
gcf_bgc_genome=cbind(rowSums(dt),gcf_bgc_genome)
write.table(gcf_bgc_genome,"table.xls",row.names = F,col.names = F,sep = "\t")
