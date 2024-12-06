library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(GSVA)

rt=read.table("matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

#seq生成步长�?2的等差序列，将正常组按列提取到前面，肿瘤组提取到后面
rt=rt[,
      c(
        seq(1,115,2) , seq(2,116,2)
        )
      ]

#判断原始数据是否去了log
max(rt)
if(max(rt)>30) rt=log2(rt+1)     #rt最大值大�?30则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))

#未标准化
cols=rainbow(ncol(rt)) ###针对24个样本，设置颜色，整体呈现彩虹色
par(cex = 0.7)
if(ncol(rt)>40) par(cex = 0.5)   ###设置字体大小
#pdf(file = "raw.pdf",width=5,height = 4)
boxplot(rt,las=2,col =cols ) ###绘图
#dev.off()

#标准�?
cols=rainbow(ncol(rt1)) ###针对24个样本，设置颜色，整体呈现彩虹色
par(cex = 0.5)
if(ncol(rt1)>40) par(cex = 0.5)   ###设置字体大小
pdf(file = "nor.pdf",width=5,height = 4.5)
boxplot(rt1,las=2,col =cols ) ###绘图
dev.off()

#保存标准化后结果
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file="norexp.txt",sep="\t",quote=F,col.names = F)

####差异分析
data=rt1
#data=rt

#控制组的数量，因为前面已经把正常组提到了最前面几列（一定要根据自己的数据集情况操作），所以直接提就行�?
afcon=63
conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#limma差异标准流程
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))
#保存所有基因的差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file="DIFF_all.xls",sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (abs(logFC)>1 & adj.P.Val < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="DIFF_af.xls",sep="\t",quote=F,col.names=F)
rt<-exp
Diff<-diff
#热图展示差异最大的�?50个基�?
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(100)){
  afGene=diffGene[c(1:50,(diffLength-50+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=exp[afGene,]
#分组标签
Type=c(rep("L",48),rep("H",47))
names(Type)=colnames(exp)
Type=as.data.frame(Type)
#分组标签的注释颜�?
anncolor=list(Type=c(H=pal_npg()(1),L=pal_npg()(2)[2]))
afExp<-apply(afExp,2,as.numeric)
pdf(file="DIFF_heatmap.pdf",height=7,width=10)
pheatmap(afExp,                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c(pal_npg()(2)[2],"white", pal_npg()(1)))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类�?
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 12,
         fontsize_row=12,
         fontsize_col=12,
         annotation_colors=anncolor
)
dev.off()

#火山图差异标准设�?
adjP=0.05
aflogFC=0
Significant=ifelse((Diff$P.Value<adjP & abs(Diff$logFC)>aflogFC), ifelse(Diff$logFC>aflogFC,"Up","Down"), "Not")
#开始绘�?
p = ggplot(Diff, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant),size=3)+
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1)))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(adjP)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash",size=1)
#添加标记，按�?
point.Pvalue=0.01
point.logFc=3
#继续绘制
Diff$symbol=rownames(Diff)
pdf("DIFF_vol.pdf",width=6.5,height=6)
p=p+theme_bw()
for_label <- Diff %>% 
  filter(abs(logFC) >point.logFc & P.Value< point.Pvalue )
p+geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",
    label.size =0.1
  )
dev.off()



####GSEA分析
deg=Diff
logFC_t=0.5
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始富集分�?
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )

class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)
#GO����
GO <- gseGO(
  geneList, #gene_fc
  ont = "ALL",
  OrgDb = org.Hs.eg.db,#����ע�ͻ���
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",#pֵУ������
)
class(GO)
colnames(GO@result)
kegg_result <- as.data.frame(GO)
rownames(GO@result)[head(order(GO@result$enrichmentScore))]
af=as.data.frame(GO@result)
write.table(af,file=paste0("2.","all_GSEAGO.xls"),sep="\t",quote=F,col.names=T)


#排序后分别取GSEA结果的前5个和�?5�?
num=5
pdf(paste0("2.","down_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
dev.off()
pdf(paste0("2.","up_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
dev.off()
#排序后取�?5个和�?5个一起展�?
num=3
pdf(paste0("2.","all_GSEA.pdf"),width = 8,height =9.5)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))],color =
c("#00FFFF","#FF4040","#8EE5EE","#00CDCD","#CD3333","#8B2323"),ES_geom = "dot",rel_heights = c(1.5, 0.3, 0.3),
base_size = 14)

            dev.off()
  num=3
            pdf(paste0("2.","all_GSEA.pdf"),width = 8,height =9.5)
            gseaplot2(GO, geneSetID = rownames(GO@result)[c(head(order(GO@result$enrichmentScore),num),tail(order(GO@result$enrichmentScore),num))],color =
                        c("#00FFFF","#8EE5EE","#FF4040","#CD3333","#00CDCD","#8B2323"),ES_geom = "dot",rel_heights = c(1.5, 0.3, 0.3),
                      base_size = 14)
            
            dev.off()
#单独展示,自行保存
gseaplot2(kk2,
          title = "name",  #设置标题
          "hsa04936", #绘制hsa04658通路的结果，通路名称与编号对�?
          color="red", #线条颜色
          base_size = 20, #基础字体的大�?
          subplots = 1:3, 
          pvalue_table = T) # 显示p�?

#山脊图，展示10个，自行保存
library(stringr)
kk2@result$Description=gsub("HALLMARK_","",kk2@result$Description)
ridgeplot(kk2,showCategory = 20)
