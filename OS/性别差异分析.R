library(ggstatsplot)
library(ggplot2)
library(ggsci)
library(ggpubr)
setwd("E:\\bladder\\sex")
load("ssGSEA.Rdata")
gseadata=list()

count = gsealist$TCGA
count=t(count)
count=as.data.frame(count)
count$Sex =unlist( clinic$TCGA$gender)
count$Sex[count$Sex %in% "male"] = "M" 
count$Sex[count$Sex %in% "female"] = "F"
gseadata[["TCGA"]] = count

geocli = read.table("E:\\bladder\\GSE13507clinical.txt",sep = "\t",header=T)
geocli = cbind(clinic$GEO[clinic$GEO$geotype %in% "Primary",],geocli$SEX)
count = gsealist$GEO[,clinic$GEO$geotype %in% "Primary"] 
count=t(count)
count=as.data.frame(count)
count$Sex = geocli$`geocli$SEX`
gseadata[["GEO"]] = count

count = gsealist$IM210
count=t(count)
count=as.data.frame(count)
count$Sex = clinic$IM210$Sex
gseadata[["IM210"]] = count


"#D51317FF" "#F39200FF" "#EFD500FF" "#95C11FFF" "#007B3DFF" "#31B7BCFF" "#0094CDFF" "#164194FF" "#6F286AFF"

for (j in names(gseadata)) {
  count=gseadata[[j]]
data = NULL
for (i in 1:29) {
  data=rbind(data,cbind(count[,i],count[,37],rep(colnames(count)[i],nrow(count))))}
colnames(data) = c("Scores","Sex","Checkpoints")
data=as.data.frame(data)
data$Scores =as.numeric(data$Scores)
p1=ggplot(data, aes(x=Checkpoints, y=Scores,fill=Sex)) +stat_compare_means( aes(label = ..p.signif..),
                                                                             method = "t.test")+
  geom_violin(trim=FALSE) + geom_boxplot(width=0.5,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_manual(values = c("#D51317FF","#0094CDFF"))+ 
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=10), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=10,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 10,face="plain"), 
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属�?
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属�?
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  xlab("")+ggtitle(j)
ggsave(paste(j,"29box.tiff"),p1,width = 12,height = 4)
data = NULL
for (i in 30:34) {
  data=rbind(data,cbind(count[,i],count[,37],rep(colnames(count)[i],nrow(count))))}
colnames(data) = c("Scores","Sex","Checkpoints")
data=as.data.frame(data)
data$Scores =as.numeric(data$Scores)
p2=ggplot(data, aes(x=Checkpoints, y=Scores,fill=Sex)) +stat_compare_means( aes(label = ..p.signif..),
                                                                            method = "t.test")+
  geom_violin(trim=FALSE) + geom_boxplot(width=0.5,position=position_dodge(0.9))+ #绘制箱线�?
  scale_fill_manual(values = c("#D51317FF","#0094CDFF"))+ 
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=6), #设置x轴刻度标签的字体显示倾斜角度�?15度，并向下调�?1(hjust = 1)，字体簇为Times大小�?20
        axis.text.y=element_text(family="Times",size=6,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 6,face="plain"), 
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属�?
                                 size=6),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属�?
                                  size=6),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  xlab("")+ggtitle(j)
ggsave(paste(j,"5box.tiff"),p2,width = 4.5,height = 4)


}

save(gseadata,file="gseadata.Rdata")
ggbetweenstats(
  data  = count,
  x     = Sex,
  y     = V36,
  subtitle = F,
  title="TCGA",
)+geom_point(aes(color = count$Sex))+
  scale_color_manual(values = col)

load("E:\\bladder\\infiltration\\9_infiltration.Rdata")

sex=cbind(colnames(inf_list$abis),unlist(clinic$TCGA$gender[clinic$TCGA$shortLetterCode %in% "TP"]))
sex=as.data.frame(sex)
meanm = matrix(NA,nrow = 9,ncol = 8)
rownames(meanm) = rownames(inf_list[[1]])
colnames(meanm) = names(inf_list)
pm = matrix(NA,nrow = 9,ncol = 8)
rownames(pm) = rownames(inf_list[[1]])
colnames(pm) = names(inf_list)


for (i in names(inf_list)) {
  m = inf_list[[i]]
  for (j in rownames(m)) {
    female = m[j,sex$V2 %in% "female"]
    male = m[j,sex$V2 %in% "male"]
    if(is.na(female)){
      pm[j,i] =NA
      meanm[j,i] = NA
    }
    else{
   re =t.test(female,male)
    pm[j,i] = re$p.value
    meanm[j,i] = mean(as.numeric(female))-mean(as.numeric(male))}
   }
}

library(pheatmap)
pmm=round(pm,digits = 2)
feature = pm
feature[pm<0.001] = "***" 
feature[0.001<pm & pm<0.01] = "**"
feature[0.01<pm & pm<0.05] = "*"
feature[!(feature %in% c("*","***","**",NA))]= pmm[!(feature %in% c("*","***","**",NA))]
bk <- c(seq(-1,-0.01,by=0.02),seq(0,1,by=0.02))
pdf("sex infiltration.pdf",width = 6,height = 5)
pheatmap(meanm,scale = "column",cluster_row = F, cluster_col = F, border=NA,
         display_numbers = feature,number_color = "black",na_col = "gray",angle_col = 45,
          color = c(colorRampPalette(colors = c("#0094CDFF" ,"white"))(length(bk)/2),colorRampPalette(colors = c("white","#D51317FF" ))(length(bk)/2))
         )
dev.off()

ggbetweenstats(
  data  = count,
  x     = Sex,
  y     = V36,
  subtitle = F,
  title="TCGA",
)+geom_point(aes(color = count$Sex))+
  scale_color_manual(values = col)
