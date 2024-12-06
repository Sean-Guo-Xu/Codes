library(amanida)
setwd("D:/学习/大创/AAA_AD/代谢组学/AAA R")

#import file
coln = c("Name", "Pvalue", "Foldchange", "Ntotal", "References")
datafile <- amanida_read("newdata.txt", mode = "quan", coln, separator = ";")

#get result
amanida_result <- compute_amanida(datafile)

#see result
amanida_result@stat

#volcano plot
volcano_plot(amanida_result, cutoff = c(0.05, 1.75))

#amanida结果输出到文件
write.csv(amanida_result@stat, file = "result_stat.txt")
write.csv(amanida_result@vote, file = "result_vote.txt")

#得到names
result_name <- read.csv("newresult_name.txt", header = TRUE)
write.csv(amanida_result@stat[1], file = "result_name.txt", row.names = FALSE)



new <- amanida_result@stat[amanida_result@stat$pval<0.05&(amanida_result@stat$fc>=1.75 | amanida_result@stat$fc<=(1/1.75)),]



write.csv(new, file = "阈值调整后.txt")
write.csv(new[,1], file='new_AAA.txt', row.names = F, quote = F)

explore_plot(amanida_result@stat, type = "mix", 1)


