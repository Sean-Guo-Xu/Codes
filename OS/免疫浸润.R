library(ggplot2)
library(tidyr)
library(immunedeconv)
library(TCGAbiolinks)
library(SummarizedExperiment)
setwd("E:\\bladder")
load("TCGA-BLCA.Rdata")
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
save(count,clinicdata,file="E:\\bladder\\TCGA_tpm.Rdata")

count = count[,see %in% "TP"]
quant = deconvolute(count, "quantiseq", tumor = TRUE)
consensus = deconvolute(count, indications = c(rep("BLCA",ncol(count))),"consensus_tme")
timer = deconvolute(count, "timer",indications=c(rep("BLCA",ncol(count))))
mcp = deconvolute(count, "mcp_counter")
write.table(count,"TCGA_tpm.txt",sep="\t",quote = F)

source("E:\\bladder\\cibersort\\CIBERSORT.R")
cibersort = CIBERSORT("TCGA_tpm.txt", "E:\\bladder\\cibersort\\ref.txt", perm=100, QN=TRUE)

xcell = deconvolute(count,"xcell")
epic = deconvolute(count,"epic")
abis = deconvolute(count,"abis")
consensus = deconvolute(count,"consensus_tme")
estimate = deconvolute(count,"estimate")

estscore = deconvolute_estimate(count)
inf_list = list(
  quantiseq = quant,
  consensus_tme = consensus,
  abis = abis,
  cibersort =cibersort
)
save(count,inf_list,estimate,file="TCGA_infiltration.Rdata")
