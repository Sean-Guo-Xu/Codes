gene_expression <- c(1, 0, 3, 2, 0, 5, 0, 4, 6, 2)
library(rstan)

stan_model_code <- "
data {
  int<lower=0> N;  // 细胞数量
  int<lower=0> y[N];  // 基因表达数据
}
parameters {
  real<lower=0> r;  // 负二项分布的形状参数
  real<lower=0, upper=1> p;  // 成功概率
}
model {
  r ~ gamma(2, 0.5);  // r 的先验分布 (Gamma 分布)
  p ~ beta(2, 2);  // p 的先验分布 (Beta 分布)
  for (n in 1:N)
    y[n] ~ neg_binomial_2(r, p);  // 似然函数
}
"
data_list <- list(N = length(gene_expression), y = gene_expression)

# 编译并运行模型
fit <- stan(model_code = stan_model_code, data = data_list, iter = 2000, chains = 4, seed = 123)

# 查看结果
print(fit)

posterior_samples <- extract(fit)
r_samples <- posterior_samples$r
p_samples <- posterior_samples$p

# 计算每个细胞基因表达量的后验概率
posterior_samples <- extract(fit)
r_samples <- posterior_samples$r
p_samples <- posterior_samples$p

# 计算每个细胞基因表达量的后验概率
posterior_probs <- sapply(1:length(gene_expression), function(i) {
  mean(dnbinom(4, size = r_samples, prob = p_samples))
})

# 查看每个细胞的概率值
posterior_probs

# 查看每个细胞的概率值
posterior_probs
