library(pscl)
library(msigdbr)
library(igraph)
library(Seurat)
library(STRINGdb)
library(doParallel)
library(foreach)

setwd("E://AAA_rupture")
load("Aorta.Rdata")
m_df = msigdbr(species = "Homo sapiens")
m_df = m_df[m_df$gs_subcat %in% c("CP:KEGG" ),]
m_df = m_df[m_df$gene_symbol %in% rownames(pbmc),]
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=400, input_directory="")
#计算通路相关基因的概率富集得分
prob_caculate<-function(zinb_data){

  # 为每个基因估计ZINB参数
  estimate_zinb_params <- function(counts) {
    # 使用zeroinfl函数估计ZINB模型参数
    fit <- zeroinfl(counts ~ 1 | 1, dist = "negbin")
    return(fit)
  }
  zero_inflation_prob <- function(zero_coeff) {
    exp(zero_coeff) / (1 + exp(zero_coeff))
  }
  # 定义EM算法以估计ZINB 参数
  
  adjusted_pnbinom <- function(counts, mu, theta, pi) {
    # 计算负二项分布的CDF
    cdf_negbin <- pnbinom(counts, size = theta, mu = mu)
    
    # 去除零膨胀的影响: 
    # 如果count = 0，CDF 应该是零膨胀概率 pi
    # 如果count > 0, CDF 应该是负二项部分的累积分布值, 排除掉零膨胀概率
    adjusted_cdf <- ifelse(counts == 0, pi, (cdf_negbin - pi) / (1 - pi))
    
    # 确保CDF在0到1之间
    adjusted_cdf <- pmin(pmax(adjusted_cdf, 0), 1)
    
    return(adjusted_cdf)
  }
  zinb_params <- foreach(i = 1:ncol(zinb_data), .packages = 'pscl') %dopar% {
    estimate_zinb_params(zinb_data[, i])
  }
  
  
  pis <- sapply(zinb_params, function(fit) zero_inflation_prob(fit$coefficients$zero))
  mus <- sapply(zinb_params, function(fit) exp(fit$coefficients$count[1]))
  thetas <- sapply(zinb_params, function(fit) fit$theta)
  
  # 计算每个基因-细胞组合的累积分布函数值
  adjusted_cdf_matrix <- matrix(0, nrow = nrow(zinb_data), ncol = ncol(zinb_data))
  for (i in 1:ncol(zinb_data)) {
    adjusted_cdf_matrix[, i] <- adjusted_pnbinom(zinb_data[, i], mus[i], thetas[i], pis[i])
  }
  colnames(adjusted_cdf_matrix) = colnames(zinb_data)
  rownames(adjusted_cdf_matrix) = rownames(zinb_data)
  return(adjusted_cdf_matrix)
}
num_cores <- 4  # 使用可用核心数减1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

numpath = length(unique(m_df$gs_name))
cellnum = length(colnames(pbmc))
score = matrix(NA,ncol = numpath,nrow = cellnum)
colnames(score) = unique(m_df$gs_name)
rownames(score) = colnames(pbmc)

for (i in 1:5) {
pathway = colnames(score)[i]
print(paste("Caculating Scores of",pathway,"(",i,"/",length(colnames(score)),")"))

gene_set =  m_df[m_df$gs_name %in% pathway,]
gene_set = gene_set$gene_symbol
mapped_genes <- string_db$map(data.frame(gene=gene_set), "gene", removeUnmappedRows = TRUE)
ppi_data <- string_db$get_interactions(mapped_genes$STRING_id)
colnames(mapped_genes)[2] = "from"
ppi_data = merge(ppi_data,mapped_genes,by="from")
colnames(mapped_genes)[2] = "to"
ppi_data = merge(ppi_data,mapped_genes,by="to")
ppi_network <- graph_from_data_frame(ppi_data[, c("gene.x", "gene.y")], directed = FALSE)
genes_in_network <- unique(c(ppi_data$gene.x, ppi_data$gene.y))
isolated_genes <- setdiff(gene_set, genes_in_network)
for (gene in isolated_genes) {
  ppi_network <- add_vertices(ppi_network, nv = 1, name = gene)
}
pagerank_scores <- page_rank(ppi_network)$vector

adjusted_cdf_matrix = prob_caculate(pbmc@assays$RNA$counts)
res = adjusted_cdf_matrix %*% pagerank_scores
score[,i] = res
}

stopCluster(cl)
save(score,file="Aorta_score.Rdata")# 关闭并行环境



source("E://scProEnrich//R//ProScore.R")
score = ProScore(pbmc@assays$RNA$counts,4,"HallMark")
pbmc$score = score[,1]
FeaturePlot(pbmc,features = "score")


setwd("E://AAA_rupture")
load("Aorta_score.Rdata")
save(score,file="Blood.score")
load("Blood.Rdata")

library(pheatmap)
anno = cbind(pbmc$sample,pbmc$cell)
anno = anno[order(factor(pbmc$sample,levels = c("rAAA","AAA","Normal"))),]
cellrank = anno[,2]
anno = anno[order(anno[,2]),]
score = score[order(factor(pbmc$sample,levels = c("rAAA","AAA","Normal"))),]
score = score[order(cellrank),]
anno = as.data.frame(anno)
pdf("heatmap.pdf",width = 12,height = 10)
pheatmap(t(score),annotation_col = anno,cluster_rows = T,cluster_cols = F,show_colnames = F)
dev.off()

load("Aorta_score.Rdata")
load("Aorta_KEGG_score.Rdata")
library(pheatmap)
anno = cbind(pbmc$sample,pbmc$cell)
anno = anno[order(factor(pbmc$sample,levels = c("rAAA","AAA","Normal"))),]
cellrank = anno[,2]
anno = anno[order(anno[,2]),]
score = score[order(factor(pbmc$sample,levels = c("rAAA","AAA","Normal"))),]
score = score[order(cellrank),]
anno = as.data.frame(anno)
pdf("heatmap.pdf",width = 12,height = 16)
pheatmap(t(score),annotation_col = anno,cluster_rows = T,cluster_cols = F,show_colnames = F,fontsize_row=4)
dev.off()

library(Seurat)
library(ggpubr)
library(ggplot2)
#VlnPlot(subset(pbmc,cell == "Myeloid"),features = "score",group.by = "sample")+  
res = NULL
for(i in unique(pbmc$cell)){
  names = colnames(pbmc)[pbmc$cell %in% i]
  partscore = score[rownames(score) %in% names,]
  mean = apply(partscore,2,mean)
  res = cbind(res,mean)
  }
colnames(res) = unique(pbmc$cell)

for (i in unique(pbmc$cell)) {
  res = res[order(res[,i],decreasing = T),] 
  toppath = rownames(res)[1:6]
  padata =NULL
  for (j in toppath) {
    partdata =cbind(pbmc$cell[pbmc$cell %in% i],pbmc$sample[pbmc$cell %in% i],score[pbmc$cell %in% i,j],rep(j,length(pbmc$cell[pbmc$cell %in% i])))
    padata = rbind(padata,partdata)
    }
  padata=as.data.frame(padata)
  padata$V3=as.numeric(padata$V3)
  colnames(padata) = c("cell","sample","score","pathway")
  padata$pathway = gsub("HALLMARK_","",padata$pathway)
  padata$sample = factor(padata$sample,levels = c("rAAA","AAA","Normal"))
ggviolin(padata,x="sample",y="score",fill = "sample",facet.by = "pathway") +
  stat_compare_means(comparisons = list(c("rAAA", "AAA"),
                                 c("rAAA", "Normal"),
                                 c("AAA", "Normal")), y_position = c(8.0, 8.5, 7.5), 
              tip_length = c(0.1, 0.05, 0.1,0.05,0.1,0.05), label = "p.signif")+
  scale_fill_manual(values = c("#D95D39", "#FF8C42", "#4A90E2"))+theme(legend.position = "None")+xlab("")+ylab("")
ggsave(paste(i,"pathway_vln.png"),height = 9,width = 8)
}

