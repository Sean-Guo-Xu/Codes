setwd("D:\\AAA\\qtl_scRNA")
ac=read.table("Artery_coloc_result.txt",header=T)
bc=read.table("Blood_coloc.txt",header=T)
as=read.table("eqtl_Artery_smr.txt",header=T)
bs=read.table("eqtl_Blood_smr.txt",header=T)
mq=read.table("mqtl_Blood_smr.txt",header = T)

ac=ac[!is.na(ac$Gene),]

bc=bc[!is.na(bc$Gene),]
ac$Symbol = ac$Gene
bc$Symbol = bc$Gene
ac$gene = ac$Symbol
bc$gene = bc$Symbol
as$gene = as$Symbol
bs$gene = bs$Symbol
mq$gene = mq$Gene
as=as[as$p_SMR<0.05,]
as=as[as$p_GWAS<0.0001,]

bs=bs[bs$p_SMR<0.05,]
bs = bs[bs$p_GWAS<0.0001,]

allgene = c(ac$Symbol,bc$Symbol,as$Symbol,bs$Symbol,mq$Gene)
allgene = unique(allgene)
allgene = allgene[!is.na(allgene)]
qtlm = matrix(0,ncol = length(allgene),nrow =5)
qtl = list(
  eqtl_Artery_coloc = ac,
  eqtl_Blood_coloc  = bc,
  eqtl_Artery_smr   = as,
  eqtl_Blood_smr    = bs,
  mqtl_Blood_smr    = mq
)
colnames(qtlm) = allgene
rownames(qtlm) = names(qtl)

for (i in names(qtl)[1:5]) {
  
  data = qtl[[i]]
  for (j in data$gene) {
    qtlm[i,j]=1
  }
  
}
qtlm = qtlm[1:5,order(apply(qtlm, 2, sum),decreasing = T)]
library(pheatmap)
pdf("eqtl_gwas_heat.pdf",width = 6,height =6)
pheatmap::pheatmap(qtlm,na_col = "gray",cluster_rows = F,cluster_cols = F,legend = F, border_color = NA,
                   color=c("#6BB7CA","#E18E6D"),
                   #鏄惁鏄剧ず姣忎釜鑹插潡瀵瑰簲鐨勬暟???(缁忓綊涓�鍖栧悗鐨勬暟???),
                   fontsize_col=10,angle_col = 315,show_rownames = T)
dev.off()
write.table(allgene,"eqtl_gene.txt",row.names = F,col.names = F,quote=F)

library(VennDiagram)
category = list(
  eqtl_Artery_Coloc = ac$gene,
  eqtl_Artery_SMR = as$gene,
  eqtl_Blood_Color = bc$gene,
  eqtl_Blood_SMR = bs$gene,
  mqtl_Blood_SMR = mq$gene
  )
pdf("venn.pdf",width = 6,height = 5)
venn.plot <- venn.diagram(
  x = category,
  alpha = 0.8,
  col= c("#6BB7CA","#E18E6D","#6BB7CA","#E18E6D","#A8D08D"),
  fill = c("#6BB7CA","#E18E6D","#6BB7CA","#E18E6D","#A8D08D"),
  filename = NULL, # 使用NULL在RStudio的查看器中直接显示图形http://127.0.0.1:24419/graphics/plot_zoom_png?width=424&height=371
  output = TRUE, # 使函数返回一个图形对象
  euler.d = TRUE, # 使用欧拉图来更准确地表示重叠区域的大小
  scaled = TRUE, # 根据集合大小调整区域大小
  cat.pos = c(0, 0, 270, 135, 0)
)
grid.draw(venn.plot)
dev.off()

