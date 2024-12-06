setwd("D:\\bigslice")
library(ggplot2)
len = read.table("11608bgc.tsv",header=T)
class = read.table("bgc_withclass.txt",sep='\t')
colnames(class)[2]="id"
data =merge(len,class,by="id")
data=cbind(data$length_nt,data$V14)
data=as.data.frame(data)
data$V1 = as.numeric(data$V1)
colnames(data)=c("Length","Type")
col=read.table("color_class.txt",sep = "\t", comment.char = "",row.names = 1)
data$Type = factor(data$Type,levels = rownames(col))
library(scales)
ggplot(data, aes(x = Length, fill = Type)) +theme_classic()+  geom_histogram(position = "stack", alpha = 0.8, bins = 250)+
  scale_y_sqrt(limits = c(0,60000),expand=c(0,2),labels=c("0","20,000","40,000","60,000"))+scale_x_continuous(expand = c(0,0),breaks = c(0,50000,100000,150000,200000,250000),limits = c(0,250000),labels = c("0","50","100","150","200","250"))+
 scale_fill_manual(values =col$V2 )+
  xlab("BGC Length (kb)")+ylab("Count")
ggsave("histogram.pdf",width = 10,height = 4)
