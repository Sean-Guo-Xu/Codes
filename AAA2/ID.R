setwd("D:/学习/大创/韩彦槊老师/转录组学/R/AAA/GSE7084")
rm(list=ls())

# Export the raw data
load(file='GSE7084_raw.Rdata')
write.csv(raw_Set, file = 'rawSet.csv')

# IDs are changed by using Excel
# 由于平台GPL2507没有对应的注释包，在上一步得到.csv未注释原始文件后，
# 使用Excel的Vlookup函数连接了annotation和raw完成注释，
# 再将表格降序排序，对第一列（基因名）删除重复数据，得到了最终的ID转换后的表格ID.txt(删除重复基因名，保留平均表达值最大的那一元组)。
# 载入ID.txt，重新定义行名称。
raw_Set <- read.csv("ID.txt", sep = "\t")
row.names(raw_Set) <- raw_Set[,1]
raw_Set <- raw_Set[,-1]
raw_Set[1:4,1:4]
group_list <- c(rep("AAA", 7), rep("Control", 8))
save(raw_Set,group_list,file='GSE7084_ID.Rdata')
load('GSE7084_ID.Rdata')

# 表达量取log2，画统计图
raw_Set <- read.csv("ID.txt", sep = "\t")
gname <- raw_Set$X
raw_Set <- raw_Set[,-1]
newSet <- log2(raw_Set)
newSet <- cbind(gname, newSet)
library(reshape2)
newSet_L <- melt(newSet)

colnames(newSet_L)=c('probe','sample','value')
newSet_L$group=rep(group_list,each=nrow(newSet))
head(newSet_L)

library(ggplot2)
p=ggplot(newSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)
p=ggplot(newSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
print(p)
p=ggplot(newSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)
p=ggplot(newSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
print(p)
p=ggplot(newSet_L,aes(value,col=group))+geom_density()
print(p)
p=ggplot(newSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="grey")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
print(p)

# load
load('GSE7084_ID.Rdata')
newSet <- log2(raw_Set)
## hclust
colnames(newSet)=paste(group_list,1:ncol(newSet),sep='_')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")
hc=hclust(dist(t(newSet)))
par(mar=c(5,5,5,10))
png('hclust.png',res=90)
plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
dev.off()

## PCA
library(ggfortify)
df=as.data.frame(t(newSet))
df$group=group_list
png('pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()

setwd("C:/Users/HP15/Documents/WeChat Files/wxid_n1j6mponrrgp22/FileStorage/File/2022-03")
kegg <- read.csv(file = "KEGGId.txt",sep = "\t")

