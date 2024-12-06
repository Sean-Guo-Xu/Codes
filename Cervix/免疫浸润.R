library(ggplot2)
library(tidyr)
library(immunedeconv)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(CIBERSORT)
library(limma)
setwd("E:\\cervix")
load("TCGA-CESC.Rdata")
count <- assay(data,"tpm_unstrand")
see  =rowData(data)
rownames(count) = see$gene_name
see = colData(data)
see = see$shortLetterCode
count = avereps(count)
clinic=colData(data)
clinic = clinic@listData
clinicdata = NULL
for (i in names(clinic)) {
  clinicdata=cbind(clinicdata,clinic[[i]])
}
colnames(clinicdata) = names(clinic)
clinic = unlist(clinic)

count = count[,see %in% "TP"]
quant = deconvolute(count, "quantiseq", tumor = TRUE)
consensus = deconvolute(count, indications = c(rep("BLCA",ncol(count))),"consensus_tme")
timer = deconvolute(count, "timer",indications=c(rep("BLCA",ncol(count))))
mcp = deconvolute(count, "mcp_counter")
lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
ciber = cibersort(lm22f, 
                        "TCGA_tpm.txt" , 
                        perm = 1000, 
                        QN = T)
write.table(count,"TCGA_tpm.txt",sep="\t",quote = F)
xcell = deconvolute(count,"xcell")
epic = deconvolute(count,"epic")
abis = deconvolute(count,"abis")
consensus = deconvolute(count,"consensus_tme")
estimate = deconvolute(count,"estimate")
ciber = t(ciber)
ciber = as.data.frame(ciber)
ciber = cbind(rownames(ciber),ciber)
colnames(ciber)[1] = "cell_type"
estscore = deconvolute_estimate(count)
inf_list = list(
  quantiseq = quant,
  cibersort = ciber,
  consensus_tme = consensus,
  abis = abis,
  timer = timer,
  xcell = xcell,
  epic = epic,
  mcp_counter = mcp,
estimate  =estimate
  )
save(inf_list,file="TCGA_infiltration.Rdata")
