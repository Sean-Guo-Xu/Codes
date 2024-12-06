setwd("D:\\bigslice")
library(ggplot2)
class = read.table("bgc_withclass.txt",sep='\t')

num=NULL
for(i in 1:26825){
  count = nrow(class[class$V4 %in% i, ])
  num = c(num,count)
}
num = as.data.frame(num)
library(scales)
ggplot(num, aes(x = num)) +theme_classic()+  geom_histogram(position = "stack", alpha = 0.8, bins = 5)+
  scale_y_sqrt(expand = c(0,0))+scale_x_log10(expand=c(0,0))+
  scale_fill_manual(values =col$V2 )+                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
  xlab("BGC Length (kb)")+ylab("Count")
ggsave("histogram.pdf",width = 10,height = 4)
