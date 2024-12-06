
setwd("E:\\AAA_gwas")

load("Europe_gwas.Rdata")
ac=read.table("Europe_Artery_coloc.txt",header=T)
bc=read.table("Europe_Blood_coloc.txt",header=T)
as=read.table("Europe_Artery_smr.txt",header=T)
bs=read.table("Europe_Blood_smr.txt",header=T)
mq=read.table("Europe_mqtl_smr.txt",header = T)
bloodgene = unique(c(bc$Gene,bs$Gene,mq$Gene))
arterygene = unique(c(ac$Gene,as$Gene))
colocgene = unique(c(ac$Gene,bc$Gene))
allgene =arterygene[arterygene %in% bloodgene]
library(MendelianRandomization)     #???ð?
library(TwoSampleMR)
allgwas = gwas
allgwas = data.frame(
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
#https://gwas.mrcieu.ac.uk/datasets/
#gene=c("EIF2B2","CETP","LDAH","GDF7","NEK9","HMGCR","CERT1","MLH3","MRC2","SH3PXD2A","LRP1")
#gene = c("LDAH","NEK9","MLH3","GDF7")
smrgene = unique(c(as$Gene,bs$Gene,mq$Gene))
gene =colocgene
gene = gene[!is.na(gene)]
library(org.Hs.eg.db)
id=mapIds(x = org.Hs.eg.db,#ע?Ͱ?
          keys = gene, #??Ҫת???Ļ???ID
          keytype = "SYMBOL", #??Ҫת????????
          column = "ENSEMBL")  #??Ҫת??Ϊ??


myeqtl = read.table("Artery_Aorta.v8.EUR.egenes.txt",header=T)
myeqtl$phenotype_id = gsub("[.].*","",myeqtl$phenotype_id)
myeqtl = myeqtl[myeqtl$phenotype_id %in% id,]
myeqtl$eaf = myeqtl$maf
myeqtl$eaf[myeqtl$ref_factor == -1] = 1 - myeqtl$eaf[myeqtl$ref_factor == -1]
myeqtl = myeqtl[!duplicated(myeqtl$variant_id),]
alleqtl = read.table("E:\\AAA_gwas\\Artery_eqtl.txt",header=T)
alleqtl = alleqtl[alleqtl$variant_id %in% myeqtl$variant_id,]
alleqtl$phenotype_id = gsub("[.].*","",alleqtl$phenotype_id)
alleqtl = alleqtl[alleqtl$phenotype_id %in% myeqtl$phenotype_id,]
myeqtl = merge(alleqtl,myeqtl,by="variant_id")
myeqtl = myeqtl[myeqtl$pval_nominal.x<5e-8,]
exposure_dat = data.frame(
  pos.exposure = myeqtl$V2.x,
  chr.exposure =myeqtl$V1,
  SNP = myeqtl$V2.y,
  beta.exposure= myeqtl$slope.x,
  se.exposure = myeqtl$slope_se.x,
  effect_allele.exposure =myeqtl$V4 ,
  other_allele.exposure = myeqtl$V3,
  eaf.exposure = myeqtl$eaf,
  pval.exposure = myeqtl$pval_nominal.x,
  samplesize.exposure= rep(387,nrow(myeqtl)),
  id.exposure = myeqtl$phenotype_id.x,
  exposure = myeqtl$phenotype_id.x
)


pheno = read.table("gwas_catalog_v1.0.2-associations_e112_r2024-05-20.tsv",header = T,sep = "\t",fill = T,na.strings = "")
pheno = pheno[pheno$SNPS %in% exposure_dat$SNP,]
outsnp = "rs11172113"
exposure_dat = exposure_dat[!(exposure_dat$SNP %in% outsnp), ]

write.csv(pheno, "Europe__aorta_confounder.result.csv", row.names=F)

#????ȥ?????????صĽ???
#delSnp=c("rs1532624", "rs12720820")     #????SNP??????(???޸?)
#exposure_dat=exposure_dat[!exposure_dat$SNP %in% delSnp,]
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
#harmonised_dat <- harmonised_dat[harmonised_dat$SNP %in% steiger$SNP[steiger$steiger_pval<0.05],]
#harmonised_dat <- harmonised_dat[harmonised_dat$id.exposure %in% direction$id.exposure[direction$steiger_pval<0.05],]

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
write.table(mr_res,"Aorta_MR_res.txt",row.names = F,col.names = T,sep = "\t",quote = F)
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
for(i in unique(mr_res[mr_res$pval<0.05,]$exposure)){
  mr_res1=mr_res[mr_res$exposure==i,]
  harmonised_dat1=harmonised_dat[harmonised_dat$exposure==i,]
  lineplot(mr_res1, harmonised_dat1)
  ggsave(paste(i,"scatter.tiff"),height = 5.5,width = 5)
}
#harmonised_dat = harmonised_dat[!(harmonised_dat$gene %in% heter_dat$exposure[heter_dat$Q_pval<0.05]),]


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
        legend.text = element_text(size = 6.5))+
  labs(x = "ln(OR)", 
       y = parse(text = "-log[10]*(italic(P)-value)")) +
  scale_size(name = "PVE",
             breaks = c(0.2*1:3)) +
  ggrepel::geom_label_repel(aes(label = label),size = 2)
ggsave("aorta_volcano.tiff",width = 11,height = 10)

mr_res2=mr_res[mr_res$exposure %in% table1$Gene,]
result_or=generate_odds_ratios(mr_res2)
write.table(result_or[,4:ncol(result_or)],"aorta_OR.txt",row.names = F,sep = "\t",quote = F)
write.table(table1,"aorta_MR_SNP_res.txt",row.names = F,col.names = T,sep = "\t")








library(grid)
library(forestploter)

mydata=read.table("OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
pdf("forest.pdf",width = 12,height = 12)
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
write.csv(pleio_dat,"aorta_pleiotropy_test.csv",row.names = F)
write.csv(heter_dat,"aorta_heterogeneity_test.csv",row.names = F)
write.csv(direction,"aorta_direction_test.csv",row.names = F)
write.csv(steiger,"aorta_steiger_test.csv",row.names = F)
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
bimr_SLE$id.exposure = rep("aaa",nrow(bimr_SLE))
bimr_SLE$exposure = rep("aaa",nrow(bimr_SLE))
id = mr_res[mr_res$pval<0.05,]$id.exposure
myeqtl = read.table("Artery_Aorta.v8.EUR.egenes.txt",header=T)
myeqtl$phenotype_id = gsub("[.].*","",myeqtl$phenotype_id)
myeqtl = myeqtl[myeqtl$phenotype_id %in% id,]
myeqtl$eaf = myeqtl$maf
myeqtl$eaf[myeqtl$ref_factor == -1] = 1 - myeqtl$eaf[myeqtl$ref_factor == -1]
myeqtl = myeqtl[!duplicated(myeqtl$variant_id),]
alleqtl = read.table("E:\\AAA_gwas\\Artery_eqtl.txt",header=T)
alleqtl = alleqtl[alleqtl$variant_id %in% myeqtl$variant_id,]
alleqtl$phenotype_id = gsub("[.].*","",alleqtl$phenotype_id)
alleqtl = alleqtl[alleqtl$phenotype_id %in% myeqtl$phenotype_id,]
myeqtl = merge(alleqtl,myeqtl,by="variant_id")
myeqtl = myeqtl[myeqtl$V2.y %in% bimr_SLE$SNP,]
eqtl =myeqtl
#setwd("E:\\smr")
#commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery/Artery_Aorta --query 1 --out eqtl_print")
#system(commond,intern = T)
#eqtl = read.table("eqtl_print.txt",header=T)

allresult = NULL

for (i in unique(id)) {

  myeqtl = eqtl[eqtl$phenotype_id.x %in% i,]
  if(nrow(myeqtl) == 0){
    next
  }

  outcome_gene = data.frame(
    pos.outcome = myeqtl$V2.x,
    chr.outcome =myeqtl$V1,
    SNP = myeqtl$V2.y,
    beta.outcome= myeqtl$slope.x,
    se.outcome = myeqtl$slope_se.x,
    effect_allele.outcome =  myeqtl$V4,
    other_allele.outcome = myeqtl$V3,
    pval.outcome = myeqtl$pval_nominal.x,
    samplesize.outcome= rep(387,nrow(myeqtl)),
    id.outcome = myeqtl$phenotype_id.y,
    outcome = myeqtl$phenotype_id.x,
    eaf.outcome = myeqtl$eaf
  )

  bimr_SLEE =bimr_SLE [bimr_SLE$SNP %in% outcome_gene$SNP,]
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
write.table(allresult[,4:ncol(allresult)],"reverse_bi_OR.txt",row.names = F,sep = "\t",quote = F)

mydata=read.table("reverse_bi_OR.txt",header = T,sep = "\t")

mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95,  mydata$or_uci95))


library(grid) 
library(forestploter)
pdf("reverse_mr.pdf",width = 10,height = 12)
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
mydata=read.table("OR.txt",header = T,sep = "\t")
mydata = mydata[mydata$exposure %in% topgene,]
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))

pdf("final forest.pdf",width = 12,height = 12)
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