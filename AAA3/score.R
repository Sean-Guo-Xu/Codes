library(msigdbr)
library(igraph)
library(Seurat)
library(STRINGdb)
setwd("E://AAA_rupture")
load("Aorta.Rdata")
m_df = msigdbr(species = "Homo sapiens")
m_df = m_df[m_df$gs_subcat %in% c("CP:KEGG" ),]
m_df = m_df[m_df$gene_symbol %in% rownames(pbmc),]
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=400, input_directory="")

pro_caculate<-function(gene_expression){
  # 假设的先验参数
  alpha_prior <- 2  # Gamma 分布的形状参数先验
  beta_prior <- 1   # Gamma 分布的比例参数先验
  alpha_prior_p <- 2  # Beta 分布的形状参数先验
  beta_prior_p <- 2   # Beta 分布的形状参数先验
  
  # 更新后得到的后验参数
  r_post <- alpha_prior + sum(gene_expression)  # r 的后验形状参数
  beta_post <- beta_prior + length(gene_expression)  # r 的后验比例参数
  alpha_post_p <- alpha_prior_p + sum(gene_expression)  # p 的后验形状参数
  beta_post_p <- beta_prior_p + length(gene_expression)  # p 的后验形状参数
  
  # 计算每个细胞的后验累积概率
  posterior_cdf <- sapply(gene_expression, function(x) {
    # 使用更新后的参数估计负二项分布的累积分布函数
    p_post <- alpha_post_p / (alpha_post_p + beta_post_p)  # 计算后验 p
    r_post <- r_post  / beta_post  # 计算后验 r
    pnbinom(x, size = r_post, prob = p_post)  # 负二项分布的累积分布函数
  })
  
  # 查看每个细胞的累积概率值
  return = posterior_cdf
}
pathway = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION" 
score<- function(pathway){ 
gene_set =  m_df[m_df$gs_name %in% pathway,]
gene_set = gene_set$gene_symbol
gene_expression <- FetchData(pbmc,vars =gene_set ,slot = "counts")  
pro = apply(gene_expression, 2, pro_caculate)
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
res = pro %*% pagerank_scores
res = apply(res, 1, sum)
pbmc$score = res
}
score(pathway )
FeaturePlot(pbmc,features = "score")
table(pbmc$score)
x = names(table(pbmc$score))
x= as.numeric(x)
plot(x,table(pbmc$score))
exp = FetchData(pbmc,vars = "CD14")
hist(exp$CD14)
hist(pbmc$score[pbmc$cell %in% "M2-M1"])
hist(pbmc$score[pbmc$cell %in% "M2-M2"])


rzinb <- function(n, lambda, k, pstr0) {
  # 生成零膨胀部分
  zero_inflated <- rbinom(n, size = 1, prob = pstr0)
  
  # 生成负二项部分
  negbin_part <- rnbinom(n, size = k, mu = lambda)
  
  # 如果是零膨胀部分则取0，否则取负二项生成的数值
  result <- ifelse(zero_inflated == 1, 0, negbin_part)
  
  return(result)
}

zinb_data <- matrix(0, nrow = N, ncol = G)
for (i in 1:G) {
  # 随机生成参数
  mu <- runif(1, 5, 20)
  size <- runif(1, 2, 10)
  pi <- runif(1, 0.2, 0.6)
  
  zinb_data[, i] <- rzinb(N, lambda = mu, k = size, pstr0 = pi)
}