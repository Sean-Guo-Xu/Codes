setwd("D:/学习/大创/R/batch_effect")
rm(list=ls())

library(mgcv)
library(nlme)
library(genefilter)
library(BiocParallel)
library(sva)
library(dplyr)


# 输入的文件为已经ID转换为symbol后且标准化后的表达矩阵，case在前，control在后，不可乱序
# grouplist为分组信息
grouplist1 <- c(rep("AAA", 7),rep("Control", 8))
grouplist2 <- c(rep("AAA", 14), rep("Control", 8))
grouplist3 <- c(rep("AAA", 49), rep("Control", 10))
s1 <- read.csv(file = "7084_exp.txt", sep = "\t", row.names = 1, header = T)
s2 <- read.csv(file = "47472_exp.txt", sep = "\t", row.names = 1, header = T)
s3 <- read.csv(file = "57691_exp.txt", sep = "\t", row.names = 1, header = T)

rmbatch <- function(s1, s2, s3, grouplist1, grouplist2, grouplist3){
  colnames(s1) <- paste(grouplist1,1:ncol(s1),sep='_')
  colnames(s2) <- paste(grouplist2,1:ncol(s2),sep='_')
  colnames(s3) <- paste(grouplist3,1:ncol(s3),sep='_')
  f1 <- as.vector(summary(as.factor(grouplist1)))
  f2 <- as.vector(summary(as.factor(grouplist2)))
  f3 <- as.vector(summary(as.factor(grouplist3)))
  ## intersect函数提取共有的symbol
  rowSym1 <- intersect(rownames(s1),rownames(s2))
  rowSym <- intersect(rowSym1, rownames(s3))
  ## 利用cbind合并
  expr <- cbind(s1[rowSym,],s2[rowSym,],s3[rowSym,])
  # 数据分别编号，第一套所有的case和control为1，第二套为2...
  batch <- paste0("batch", rep(c(1,2,3),c(ncol(s1),ncol(s2),ncol(s3))))
  # 数据分别的case,control，第一套有x1个case，y1个control，第二套有x2个case，y2个control...
  tissue <- rep(rep(c("cse","ctl"),3), c(f1, f2, f3))
  mod <- model.matrix(~tissue)
  expr_batch <- ComBat(dat = expr, batch = batch, mod = mod)
  # 整理数据，把两组的case与control分离，重新排列
  expr_batch <- data.frame(expr_batch[,1 : f1[1]],
                           expr_batch[,(f1[1]+f1[2]+1) : (f1[1]+f1[2]+f2[1])],
                           expr_batch[,(f1[1]+f1[2]+f2[1]+f2[2]+1) : (f1[1]+f1[2]+f2[1]+f2[2]+f3[1])],
                           expr_batch[,(f1[1]+1) : (f1[1]+f1[2])],
                           expr_batch[,(f1[1]+f1[2]+f2[1]+1) : (f1[1]+f1[2]+f2[1]+f2[2])],
                           expr_batch[,(f1[1]+f1[2]+f2[1]+f2[2]+f3[1]+1) : (f1[1]+f1[2]+f2[1]+f2[2]+f3[1]+f3[2])])
  colnames(expr_batch) <- c(paste(grouplist1[f1[1]], 1 : f1[1], sep="_"),
                            paste(grouplist2[f2[1]], (f1[1]+1) : (f1[1]+f2[1]), sep="_"),
                            paste(grouplist3[f3[1]], (f1[1]+f2[1]+1) : (f1[1]+f2[1]+f3[1]), sep="_"),
                            paste(grouplist1[(f1[1]+1):length(grouplist1)], 1 : f1[2], sep="_"),
                            paste(grouplist2[(f2[1]+1):length(grouplist2)], (f1[2]+1) : (f1[2]+f2[2]), sep="_"),
                            paste(grouplist3[(f3[1]+1):length(grouplist3)], (f1[2]+f2[2]+1) : (f1[2]+f2[2]+f3[2]), sep="_"))
  return(expr_batch)
}

expr_batch <- rmbatch(s1, s2, s3, grouplist1, grouplist2, grouplist3)
write.table(expr_batch, file = "3GSE_batchrm_expr.txt", sep = "\t", quote = F)
# 完成后记得在结果文件的首部加个"\t"!!!!!!!!!!!!!!!!!!!!
