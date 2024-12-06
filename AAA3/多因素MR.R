library(TwoSampleMR)
library(tidyverse)
setwd("E:\\AAA_gwas")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
see=available_outcomes()
gwasname = see[grep("HLA",see$trait),]
unique(gwasname$author)
gwasname = gwasname[gwasname$author %in% "Sun BB",]
load("Europe_gwas.Rdata")
allgwas = gwas 
gwas = data.frame(
  SNP = allgwas$V5,
  chr  = allgwas$V1,
  pos = allgwas$V2,
  beta.outcome= allgwas$V9,
  se.outcome = allgwas$V10,
  samplesize.outcome= rep(429209,nrow(allgwas)) ,
  pval.outcome = allgwas$V7,
  eaf.outcome = allgwas$V11,
  effect_allele.outcome =allgwas$V4,
  other_allele.outcome = allgwas$V3,
  outcome= rep("AAA",nrow(allgwas)),
  id.outcome= rep("AAA",nrow(allgwas))
)
exposure_dat <-extract_instruments(gwasname$id)
outcome_dat  <- gwas[gwas$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat,  outcome_dat)
mr_results <- mr(dat)


exposure_dat <-mv_extract_exposures(gwasname$id[gwasname$id %in% unique(mr_results$id.exposure)])
outcome_dat  <- gwas[gwas$SNP %in% exposure_dat$SNP,]
mvdat <- mv_harmonise_data(exposure_dat,  outcome_dat)
res <- mv_multiple(mvdat)
write.csv(res,"mv_mr_res.csv")

write.csv(mr_results,"mr_res.csv")
