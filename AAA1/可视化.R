library(Seurat)
library(RColorBrewer)
load("F:\\AAA_scRNA\\7AAA.Rdata")
setwd("F:\\AAA_scRNA\\scDRS")
score = read.table("Aorta.score",sep = "\t",header = T)
meta = data.frame(pbmc@meta.data, pbmc@reductions$umap@cell.embeddings)
meta$X = rownames(meta)
meta = merge(meta,score,by="X")
#meta$norm_score[meta$pval >0.05] = NA("#62B197", "#62B197")
lowcolor ="blue3"
highcolor =  "red3"
library(ggplot2)

  conver = function(x){
    x[x %in% 0] = NA
    return(x)
  }
  data = data.frame(meta$norm_score, meta$UMAP_1,meta$UMAP_2)
  data=apply(data, 2, conver)
  data = as.data.frame(data)
  colnames(data)=c("drsscore","UMAP_1","UMAP_2")
  ggplot() +  geom_point(data=data[is.na(data[,1]),],aes(UMAP_1, UMAP_2),color = "white", size=1.5) +geom_point(data=data[!is.na(data[,1]),],aes(UMAP_1, UMAP_2, color=drsscore), size=1.5)+scale_colour_gradient2(low = lowcolor, mid = "white",high = highcolor)  + theme_light(base_size = 15)+labs(title = "")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+labs(color = "Expression")+
    theme(plot.title = element_text(hjust = 0.5))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank(),axis.text.x = element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) + xlab("")+ylab("")+ggtitle("Scores in Aorta") +labs(color = "scDRS Score")
ggsave("Aorta_DRS_feature.tiff",width = 6,height =5)

data = cbind(meta$norm_score,meta$cell)
data = as.data.frame(data)
data$V1 = as.numeric(data$V1)
colnames(data) =  c("drsscore","cell")
mean_scores <- data %>%
  group_by(cell) %>%
  summarise(mean_score = mean(drsscore)) %>%
  arrange(desc(mean_score))
data$cell <- factor(data$cell, levels = mean_scores$cell)
color_palette <-color_palette <- rev(brewer.pal(11, "Reds"))
ggplot(data, aes(x = cell, y = drsscore, fill = cell)) +
  geom_violin(trim = FALSE, alpha = 0.8) +  # 绘制小提琴图，设置透明度
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +  # 在小提琴图上叠加箱线图，不显示离群值
  stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", ref.group = ".all.") +  # # 添加标题和轴标签
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中并加粗
    axis.title = element_text(face = "bold"),  # 坐标轴标签加粗
    axis.text = element_text(size = 12),  # 坐标轴刻度字体大小
    legend.position = "none",  # 移除图例
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.border = element_blank(),  # 移除面板边框
    axis.text.x = element_text(angle = 45, hjust = 1)
    )+
  labs(title = "", x = "", y = "DRS Scores") +  # 添加标题和轴标签
  theme(legend.position = "none") + scale_fill_manual(values = color_palette) 
ggsave("Aorta_DRS_vlnplot.tiff",width = 8,height =5)
