
setwd("E:\\AAA_gwas")
or = read.table("Asian_OR.txt",sep = "\t",header = T)
load("Europe_gwas.Rdata")
library(MendelianRandomization)     #???ð?
library(TwoSampleMR)
allgwas = gwas
allgwas = data.frame(
  SNP = allgwas$V5,
  chr  = allgwas$V1,
  pos = allgwas$V2,
  beta.outcome= allgwas$V9,
  se.outcome = allgwas$V10,
  samplesize.outcome= rep(31684,nrow(allgwas)) ,
  pval.outcome = allgwas$V7,
  eaf.outcome = allgwas$V11,
  effect_allele.outcome =allgwas$V4,
  other_allele.outcome = allgwas$V3,
  outcome= rep("AAA",nrow(allgwas)),
  id.outcome= rep("AAA",nrow(allgwas))
)
#https://gwas.mrcieu.ac.uk/datasets/
#gene=c("EIF2B2","CETP","LDAH","GDF7","NEK9","HMGCR","CERT1","MLH3","MRC2","SH3PXD2A","LRP1")
#gene = c("LDAH","NEK9","MLH3","GDF7")
gene =unique(or$exposure)
gene = gene[!is.na(gene)]
library(org.Hs.eg.db)
id=mapIds(x = org.Hs.eg.db,#ע?Ͱ?
          keys = gene, #??Ҫת???Ļ???ID
          keytype = "SYMBOL", #??Ҫת????????
          column = "ENSEMBL")  #??Ҫת??Ϊ??


gene = cbind(gene,id)
exposure_id=paste0('eqtl-a-',id)

exposure_dat <- extract_instruments(exposure_id, p1=5e-08, clump=TRUE)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
exposure_dat <- extract_instruments(exposure_id, p1=5e-08, clump=TRUE)
library(phenoscanner)
#??SNP????
snpId=exposure_dat$SNP
y=seq_along(snpId)
chunks <- split(snpId, ceiling(y/100))
outTab=data.frame()
for(i in names(chunks)){
  #???????ط???
  confounder=phenoscanner::phenoscanner(
    snpquery = chunks[[i]],
    catalogue = "GWAS",
    pvalue = 1e-05,
    proxies = "None",
    r2 = 0.8,
    build = 37)
  outTab=rbind(outTab, confounder$results)
}

write.csv(outTab, "confounder.result.csv", row.names=F)
delSnp=c("rs1532624", "rs12720820")     #????SNP??????(???޸?)
exposure_dat=exposure_dat[!exposure_dat$SNP %in% delSnp,]

##################################################
R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2=R2a/(R2a+R2b)
exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)
exposure_dat=exposure_dat[exposure_dat$F_statistics>10,]
exposure_dat$exposure=stringr::str_sub(exposure_dat$exposure,1,15)

#########SNP????????######
outcome_dat = allgwas[allgwas$SNP %in% exposure_dat$SNP,]
direction =directionality_test(harmonise_data(exposure_dat, outcome_dat))
steiger<-steiger_filtering(harmonise_data(exposure_dat, outcome_dat))
library(tidyverse)
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]
exposure_dat$ENSEMBL =gsub("eqtl-a-","",exposure_dat$id.exposure)

id=mapIds(x = org.Hs.eg.db,#ע?Ͱ?
          keys = exposure_dat$ENSEMBL, #??Ҫת???Ļ???ID
          keytype = "ENSEMBL", #??Ҫת????????
          column = "SYMBOL")

exposure_dat$exposure  =id
#mv_harmonise_data(exposure_dat,  outcome_dat)
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)
harmonised_dat <- harmonised_dat[harmonised_dat$SNP %in% steiger$SNP[steiger$steiger_pval<0.05],]
harmonised_dat <- harmonised_dat[harmonised_dat$id.exposure %in% direction$id.exposure[direction$steiger_pval<0.05],]

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}
mr_res <- mr_modified(harmonised_dat, prop_var_explained = T)
mr_res <- mr_res[mr_res$method %in% c("Inverse variance weighted","Wald ratio"),]
heter_dat=mr_heterogeneity(harmonised_dat)
pleio_dat=mr_pleiotropy_test(harmonised_dat)


table1 <- mr_res %>% 
  filter(pval < 1,
         method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  left_join(exposure_dat, by = "exposure")

table1 <- table1 %>% 
  generate_odds_ratios()%>% 
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
         `P value` = scales::scientific(pval),
         `PVE` = paste0(sprintf("%.2f",100*pve),"%"),
         `F statistics` = sprintf("%.2f",F_statistics)) %>% 
  dplyr::select(Gene = exposure, `ENSEMBL ID` = ENSEMBL,
                SNP, `Effect allele` = effect_allele.exposure, 
                `OR (95% CI)`, `P value`, 
                PVE, `F statistics`)


data = mr_res %>% 
  filter(method %in% c("Wald ratio","Inverse variance weighted")) 
p_thershold <- 0.05
data = cbind(data$b,-log10(data$pval),ifelse(data$pval < p_thershold, data$exposure, NA), data$pve)
data = as.data.frame(data)
colnames(data) = c("x","y","label","size")
data[,1]  =as.numeric(data[,1])
data$y = as.numeric(data$y)
data$size = as.numeric(data$size)
library(ggplot2)
ggplot(data,aes( x = x, y = y)) +
  geom_point(data = data[data$x>0,],aes(size = size), alpha = 1, color = "#E07B54") +
  geom_point(data = data[data$x<0,],aes(size = size), alpha = 1, color = "#6BB7CA") +
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6.5),
        legend.position = legend.position)+
  labs(x = "ln(OR)", 
       y = parse(text = "-log[10]*(italic(P)-value)"),
       title = title) +
  scale_size(name = "PVE",
             breaks = c(0.2*1:3)) +
  ggrepel::geom_label_repel(aes(label = label),size = 3)
ggsave("volcano.tiff",width = 5,height = 5)

mr_res2=mr_res[mr_res$exposure %in% table1$Gene,]
result_or=generate_odds_ratios(mr_res2)
write.table(result_or[,4:ncol(result_or)],"Europe_test_OR.txt",row.names = F,sep = "\t",quote = F)









library(grid)
library(forestploter)

mydata=read.table("Europe_test_OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
pdf("Europe_test_forest.pdf",width = 12,height = 12)
colnames(mydata)[c(1:3,6)]= c("Exposure","Method","nSNP","Pval")
forest(data=mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 2),
)
dev.off()
write.csv(pleio_dat,"Europe_test_pleiotropy_test.csv",row.names = F)
write.csv(heter_dat,"Europe_test_heterogeneity_test.csv",row.names = F)
write.csv(direction,"Europe_test_direction_test.csv",row.names = F)
write.csv(steiger,"Europe_test_steiger_test.csv",row.names = F)
#########????????
library(TwoSampleMR)
bimr_SLE = data.frame(
  pos = allgwas$pos,
  chr = allgwas$chr,
  SNP = allgwas$SNP,
  beta.exposure= allgwas$beta.outcome,
  se.exposure = allgwas$se.outcome,
  effect_allele.exposure =allgwas$effect_allele.outcome ,
  other_allele.exposure = allgwas$other_allele.outcome,
  eaf.exposure = allgwas$eaf.outcome,
  pval.exposure = allgwas$pval.outcome,
  samplesize.exposure= rep("31684",nrow(allgwas)) 
)
#bimr_SLE<- extract_instruments('ebi-a-GCST90011866', p1=5e-08, clump=TRUE)

bimr_SLE = bimr_SLE[bimr_SLE$pval.exposure< 5e-8,]
bimr_SLE$id.exposure = rep("AAA",nrow(bimr_SLE))
bimr_SLE$exposure = rep("AAA",nrow(bimr_SLE))

allresult = NULL
id = mr_res[mr_res$pval<0.05,]$id.exposure

for (i in unique(id)) {
  
  outcome_gene<- extract_outcome_data(snps=bimr_SLE$SNP, outcomes=i)
  if(is.null(outcome_gene)){
    next
  }
  bimr_SLEE =bimr_SLE [bimr_SLE $SNP %in% outcome_gene$SNP,]
  harmonised_SLE_gene <- harmonise_data(bimr_SLEE, outcome_gene)
  bimr_mr_SLE_gene <- mr(harmonised_SLE_gene)
  result_or=generate_odds_ratios(bimr_mr_SLE_gene)
  allresult = rbind(allresult,result_or)
}
allresult$id.outcome = gsub("eqtl-a-","",allresult$id.outcome)
library(org.Hs.eg.db)
id=mapIds(x = org.Hs.eg.db,#ע?Ͱ?
          keys =allresult$id.outcome, #??Ҫת???Ļ???ID
          keytype = "ENSEMBL", #??Ҫת????????
          column = "SYMBOL")  #??Ҫת??Ϊ??
allresult$id.outcome = id
allresult$Outcome = id
write.table(allresult[,4:ncol(allresult)],"Europe_test_reverse_bi_OR.txt",row.names = F,sep = "\t",quote = F)

mydata=read.table("Europe_test_reverse_bi_OR.txt",header = T,sep = "\t")

mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95,  mydata$or_uci95))


library(grid) 
library(forestploter)
pdf("Europe_test_reverse_mr.pdf",width = 10,height = 12)
mydata$exposure = "AAA"
mydata=mydata[mydata$method %in% c("Inverse variance weighted"),]
colnames(mydata)[c(1:2,6)] = c("Exposure","Method","Pval")
forest(mydata[,c(1:2,6,12,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5,
       ref_line = 1,
       xlim = c(0.05, 3),
)
dev.off()

topgene = mydata[mydata$Pval>0.05,colnames(mydata) %in% "Outcome"]
mydata=read.table("Europe_test_OR.txt",header = T,sep = "\t")
mydata = mydata[mydata$exposure %in% topgene,]
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
pdf("Europe_test_final forest.pdf",width = 12,height = 12)
colnames(mydata)[c(1:3,6)]= c("Exposure","Method","nSNP","Pval")
forest(data=mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 2),
)
dev.off()                      
