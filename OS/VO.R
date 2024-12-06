data<-diff
library(ggplot2)
setwd("D:\\生信\\target数据分析\\vol")
data$Significant<- factor(ifelse(data$FDR < 0.05 & abs(data$logFC) > 2, 
                              ifelse(data$logFC > 2,'Up','Down'),'Not'))
p = ggplot(data, aes(logFC, -log10(FDR)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("cyan", "black", "pink"))+
  labs(title = " ")+
  geom_hline(yintercept = -log10(0.05),linetype=2,cex=1)+  #添加辅助线
  geom_vline(xintercept = c(-2,2),linetype=2,cex=1)+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
p
pdf("circRNA.vol.pdf",width=5.5,height=4)
print(p)
dev.off()



