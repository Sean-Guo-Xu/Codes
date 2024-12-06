setwd("D:\\AAA\\MR")
library(VariantAnnotation)
library(gwasglue)
library(dplyr)
library(tidyr)
library(CMplot)
library(TwoSampleMR)
library(MendelianRandomization) 
library(TwoSampleMR)
see=available_outcomes()
inputFile="prot-a-2300.vcf.gz"     #输入文件(根据下载暴露数据的文件名称进行修改)
     #设置工作目录

#读取输入文件, 并对输入文件进行格式转换
vcfRT <- readVcf(inputFile)
alldata=gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
setwd("D:\\AAA\\MR\\prot")
#根据pvalue<5e-08对结果进行过滤
outTab<-subset(alldata, pval.exposure<5e-06)
write.csv(outTab, file="exposure.pvalue.csv", row.names=F)

#准备绘制曼哈顿图的数据
data=alldata[,c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
colnames(data)=c("SNP","CHR","BP","pvalue")

#绘制线性的曼哈顿图
CMplot(data,  plot.type="m",
       LOG10=TRUE, threshold=5e-08, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,
       chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,50),
       file="pdf", file.output=TRUE, width=15, height=9, verbose=TRUE)

exposure_dat<-read_exposure_data(filename="exposure.pvalue.csv",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 eaf_col = "eaf.exposure",
                                 samplesize_col = "samplesize.exposure",
                                 clump = F)

#去除连锁不平衡的SNP

exposure_dat_clumped <- clump_data(outcome_dat, clump_kb=1000, clump_r2=0.2)
write.csv(exposure_dat_clumped, file="exposure.LD.csv", row.names=F)
Ffilter=10        #F鍊艰繃婊ゆ潯浠?
dat<-read.csv("exposure.LD.csv", header=T, sep=",", check.names=F)

N=dat[1,"samplesize.exposure"]    
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))     #璁＄畻R2
dat=transform(dat,F=beta.exposure^2/se.exposure^2)      #璁＄畻F妫�楠屽�?
outTab=dat[dat$F>Ffilter,]
write.csv(dat, "exposure.F.csv", row.names=F)

    #引用包        #输入文件
dat=read.csv("exposure.F.csv" , header=T, sep=",", check.names=F)

#对SNP分组
snpId=dat$SNP
y=seq_along(snpId)
chunks <- split(snpId, ceiling(y/100))

#对分组进行循环,每次得到一个分组
outTab=data.frame()
for(i in names(chunks)){
  #混杂因素分析
  confounder=phenoscanner(
    snpquery = chunks[[i]],
    catalogue = "GWAS",
    pvalue = 1e-05,
    proxies = "None",
    r2 = 0.8,
    build = 37)
  outTab=rbind(outTab, confounder$results)
}
#输出SNP相关性状的表格
write.csv(outTab, "confounder.result.csv", row.names=F)

delSnp=unique(outTab$snp)    #混杂SNP的名称(需修改)
dat=dat[!dat$SNP %in% delSnp,]
write.csv(dat, "exposure.confounder.csv", row.names=F)

exposureName="PLEK"                         #图形中展示暴露数据的名称
outcomeName="AAA"       #图形中展示结局数据的名称

exposureFile="exposure.confounder.csv"     #暴露数据输入文件

outcomeFile="ieu-a-7.vcf.gz"               #结局数据输入文件(改成我们下载的结局数据文件)

#读取暴露数据
exposure_dat<-read_exposure_data(filename=exposureFile,
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 eaf_col = "eaf.exposure",
                                 clump = F)
out = read.table("D:\\AAA\\MR\\Finngen_R9_I9_ABAORTANEUR",sep="\t")
#读取结局数据的vcf文件,并对数据进行格式转换
outcome_dat = format_data(
  out,
  type='outcome',
  snp_col = "V5",
  beta_col = "V9",
  se_col = "V10",
  effect_allele_col ="V4",
  other_allele_col = "V3",
  eaf_col = "V11",
  pval_col = "V7"
)

#将暴露数据和结局数据合并
exposure_dat$exposure=exposureName
outcome_dat$outcome=outcomeName
dat<-harmonise_data(exposure_dat=exposure_dat,
                    outcome_dat=outcome_dat)

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#MR-PRESSO异常值检测(偏倚的SNP)
presso=run_mr_presso(dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file="table.MR-PRESSO.csv")

#孟德尔随机化分析
mrResult=mr(dat)
#选择孟德尔随机化的方法
#mr_method_list()$obj
#mrResult=mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)
save(dat,mrResult,file="MR ATHSCLE.Rdata")
library(ggplot2)
#绘制散点图
scatter=function (mr_results, dat) 
{
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - 
                                                                                                                           se.outcome, ymax = beta.outcome + se.outcome), 
                                                                                                            colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - 
                                                                                                                                                                                 se.exposure, xmax = beta.exposure + se.exposure), 
                                                                                                                                                                  colour = "grey", height = 0) + ggplot2::geom_point() + 
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                           slope = b, colour = method), show.legend = TRUE) + theme_classic()+
                           ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                                                   "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                                                   "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test", 
                                                                                                                     x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", 
                                                                                                                                                                          d$outcome[1])) + ggplot2::theme(legend.position = "top", 
                                                                                                                                                                                                          legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                       })
  mrres
}
scatter(mrResult[mrResult$method%in%c("MR Egger","Inverse variance weighted"),],dat)
ggsave(paste(exposureName,"scatter.tiff"),width = 5,height = 5.5)

#森林图
library(ggplot2)
forest=function (singlesnp_results, exponentiate = FALSE) 
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), 
                     function(d) {
                       d <- plyr::mutate(d)
                       if (sum(!grepl("All", d$SNP)) < 2) {
                         return(blank_plot("Insufficient number of SNPs"))
                       }
                       levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
                       levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
                       am <- grep("All", d$SNP, value = TRUE)
                       d$up <- d$b + 1.96 * d$se
                       d$lo <- d$b - 1.96 * d$se
                       d$tot <- 0.01
                       d$tot[d$SNP %in% am] <- 1
                       d$SNP <- as.character(d$SNP)
                       nom <- d$SNP[!d$SNP %in% am]
                       nom <- nom[order(d$b)]
                       d <- rbind(d, d[nrow(d), ])
                       d$SNP[nrow(d) - 1] <- ""
                       d$b[nrow(d) - 1] <- NA
                       d$up[nrow(d) - 1] <- NA
                       d$lo[nrow(d) - 1] <- NA
                       d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                       xint <- 0
                       if (exponentiate) {
                         d$b <- exp(d$b)
                         d$up <- exp(d$up)
                         d$lo <- exp(d$lo)
                         xint <- 1
                       }
                       ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + 
                         ggplot2::geom_vline(xintercept = xint, linetype = "dotted") + 
                         ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                              xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                 height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                         ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                               "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                "red")) + ggplot2::scale_size_manual(values = c(0.3, 
                                                                                                                                                                                                1)) + ggplot2::theme(legend.position = "none", 
                                                                                                                                                                                                                     axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                     axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                     axis.title.x = ggplot2::element_text(size = 8)) + 
                         ggplot2::labs(y = "", x = paste0("MR effect size for\n'", 
                                                          d$exposure[1], "' on '", d$outcome[1], "'"))+theme_classic()
                     })
  res
}
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
forest(res_single)
ggsave(paste(exposureName,"forest.tiff"),width = 6,height = 6)
#漏斗图
funnel = function (singlesnp_results) 
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), 
                     function(d) {
                       d <- plyr::mutate(d)
                       if (sum(!grepl("All", d$SNP)) < 2) {
                         return(blank_plot("Insufficient number of SNPs"))
                       }
                       am <- grep("All", d$SNP, value = TRUE)
                       d$SNP <- gsub("All - ", "", d$SNP)
                       am <- gsub("All - ", "", am)
                       ggplot2::ggplot(subset(d, !SNP %in% am), ggplot2::aes(y = 1/se, 
                                                                             x = b)) + ggplot2::geom_point() + ggplot2::geom_vline(data = subset(d, 
                                                                                                                                                 SNP %in% am), ggplot2::aes(xintercept = b, colour = SNP)) + 
                         ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                                                 "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                                                 "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                                                 "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(y = expression(1/SE[IV]), 
                                                                                                                   x = expression(beta[IV]), colour = "MR Method") + 
                         ggplot2::theme(legend.position = "top", legend.direction = "vertical")+theme_bw()
                     })
  res
}
funnel(singlesnp_results = res_single)
ggsave(paste(exposureName,"funnel.tiff"),width = 6,height = 5)

#留一法敏感性分析
leaveoneout = function (leaveoneout_results) 
{
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.3, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), axis.title.x = ggplot2::element_text(size = 8)) + 
                                                ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", 
                                                                                 d$exposure[1], "' on '", d$outcome[1], "'"))+theme_classic()
                                            })
  res
}


leaveoneout(leaveoneout_results = mr_leaveoneout(dat))
ggsave(paste(exposureName,"leaveoneout.tiff"),width = 6,height = 6)

