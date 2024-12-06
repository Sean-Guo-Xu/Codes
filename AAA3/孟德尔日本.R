
setwd("E:\\AAA_gwas")
load("asian_gwas.Rdata")
ac=read.table("Asian_Artery_coloc.txt",header=T)
bc=read.table("Asian_Blood_coloc.txt",header=T)
as=read.table("Asian_Artery_smr.txt",header=T)
bs=read.table("Asian_Blood_smr.txt",header=T)
mq=read.table("Asian_mqtl_smr.txt",header = T)
bloodgene = unique(c(bc$Gene,bs$Gene,mq$Gene))
arterygene = unique(c(ac$Gene,as$Gene))
colocgene = unique(c(ac$Gene,bc$Gene))
allgene =arterygene[arterygene %in% bloodgene]
setwd("E:\\AAA_gwas")
library(MendelianRandomization)     #???ð?
library(TwoSampleMR)
allgwas = data.frame(
  SNP = gwas$SNPID,
  chr  = gwas$CHR,
  pos = gwas$POS,
  beta.outcome= gwas$BETA,
  se.outcome = gwas$SE,
  samplesize.outcome= gwas$N,
  pval.outcome = gwas$p.value,
  eaf.outcome = gwas$AF_Allele2,
  effect_allele.outcome = gwas$Allele2,
  other_allele.outcome = gwas$Allele1,
  outcome= rep("AAA",nrow(gwas)),
  id.outcome= rep("AAA",nrow(gwas))
)
#https://gwas.mrcieu.ac.uk/datasets/
#gene=c("EIF2B2","CETP","LDAH","GDF7","NEK9","HMGCR","CERT1","MLH3","MRC2","SH3PXD2A","LRP1")
#gene = c("LDAH","NEK9","MLH3","GDF7")
gene =colocgene
gene = gene[!is.na(gene)]
library(org.Hs.eg.db)
id=mapIds(x = org.Hs.eg.db,#ע?Ͱ?
          keys = gene, #??Ҫת???Ļ???ID
          keytype = "SYMBOL", #??Ҫת????????
          column = "ENSEMBL")  #??Ҫת??Ϊ??


gene = cbind(gene,id)
exposure_id=paste0('eqtl-a-',id)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
exposure_dat <- extract_instruments(exposure_id, p1=5e-08, clump=TRUE)

pheno = read.table("gwas_catalog_v1.0.2-associations_e112_r2024-05-20.tsv",header = T,sep = "\t",fill = T,na.strings = "")
pheno = pheno[pheno$SNPS %in% exposure_dat$SNP,]
write.csv(pheno, "Asian_confounder.result.csv", row.names=F)
#????ȥ?????????صĽ???
delSnp=c("rs1532624", "rs12720820")     #????SNP??????(???޸?)
exposure_dat=exposure_dat[!exposure_dat$SNP %in% delSnp,]
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

#绘制散点图
lineplot=function (mr_results, dat) 
{
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
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
                                                                           slope = b, colour = method), show.legend = TRUE) + ggplot2::theme_classic()+
                           ggplot2::scale_colour_manual(values = c("#E07B54", 
                                                                   "#6BB7CA", "#b2df8a", "#33a02c", "#fb9a99", 
                                                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                                                   "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test", 
                                                                                                                     x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", 
                                                                                                                                                                          d$outcome[1])) + ggplot2::theme(legend.position = "top", 
                                                                                                                                                                                                          legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                       })
  mrres
}
for(i in unique(mr_res[mr_res$pval<0.052,]$exposure)){
  mr_res1=mr_res[mr_res$exposure==i,]
  harmonised_dat1=harmonised_dat[harmonised_dat$exposure==i,]
  lineplot(mr_res1, harmonised_dat1)
  ggsave(paste(i,"scatter.tiff"),height = 5.5,width = 5)
}
harmonised_dat = harmonised_dat[!(harmonised_dat$gene %in% heter_dat$exposure[heter_dat$Q_pval<0.05]),]


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
write.table(result_or[,4:ncol(result_or)],"Asian_OR.txt",row.names = F,sep = "\t",quote = F)









library(grid)
library(forestploter)

mydata=read.table("Asian_OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
pdf("Asian_forest.pdf",width = 12,height = 12)
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
write.csv(pleio_dat,"Asian_pleiotropy_test.csv",row.names = F)
write.csv(heter_dat,"Asian_heterogeneity_test.csv",row.names = F)
write.csv(direction,"Asian_direction_test.csv",row.names = F)
write.csv(steiger,"Asian_steiger_test.csv",row.names = F)
#########????????
library(TwoSampleMR)
bimr_SLE = data.frame(
  pos = gwas$POS,
  chr = gwas$CHR,
  SNP = gwas$SNPID,
  beta.exposure= gwas$BETA,
  se.exposure = gwas$SE,
  effect_allele.exposure =gwas$Allele2 ,
  other_allele.exposure = gwas$Allele1,
  eaf.exposure = gwas$AF_Allele2,
  pval.exposure = gwas$p.value,
  samplesize.exposure= gwas$N) 
)
#bimr_SLE<- extract_instruments('ebi-a-GCST90011866', p1=5e-08, clump=TRUE)

bimr_SLE = bimr_SLE[bimr_SLE$pval.exposure< 5e-6,]
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
write.table(allresult[,4:ncol(allresult)],"Asian_reverse_bi_OR.txt",row.names = F,sep = "\t",quote = F)

mydata=read.table("Asian_reverse_bi_OR.txt",header = T,sep = "\t")

mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95,  mydata$or_uci95))


library(grid) 
library(forestploter)
pdf("Asian_reverse_mr.pdf",width = 10,height = 12)
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
mydata=read.table("Asian_OR.txt",header = T,sep = "\t")
mydata = mydata[mydata$exposure %in% topgene,]
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
pdf("Asian_final forest.pdf",width = 12,height = 12)
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
