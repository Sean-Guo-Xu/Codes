setwd("E:\\cervix")
library(survival)
library(survminer)
library(ggprism)
library(glmnet)
library(TCGAbiolinks)
library(SummarizedExperiment)
coxPfilter=0.05          #??????cox?????????Թ??˱?׼     #???ù???Ŀ¼

#??????cox????
load("TCGA-CESC.Rdata")
count <- assay(data,"unstranded")
see  =rowData(data)
gene = c("EFEMP1","MICB","HDAC10")
rownames(count) = see$gene_name
see = colData(data)
over1 = see$days_to_death
over2 = see$days_to_last_follow_up
over1[is.na(over1)] = over2[is.na(over1)]
time  = cbind(see$barcode,see$vital_status,over1,see$shortLetterCode)
time = time[!(time[,3] %in% c(NA,"[Discrepancy]","[Not Available]")),]
time = time[time[,3] >= 0,]
count = count[,colnames(count) %in% time[,1]]
time = time[order(time[,1]),]
count = count[,order(colnames(count))]
vm = read.table("E:\\cervix\\vm.txt")
count = count[rownames(count) %in% vm$V1,]
count = t(count)
time[time[,2] %in% "Alive",2] = 0
time[time[,2] %in% "Dead",2] = 1
time =as.data.frame(time)
time$V2 = as.numeric(time$V2)
time$over1 = as.numeric(time$over1)
count = cbind(time[,2:3],count)
colnames(count)[1:2] = c("fustat","futime")
rt = count
rt = as.data.frame(rt)

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,i])<0.001){next}
  #??????cox????
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  sigGenes=NULL
    #?????????Ի???
  if(coxP<1){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
outTab$pvalue = as.numeric(outTab$pvalue)
coxgene = list(TCGA  = outTab[outTab$pvalue<0.05 ,])
write.table(outTab,file="TCGA_uni.Cox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="TCGA.SigExp.txt",sep="\t",row.names=F,quote=F)
outTab1=outTab
#############GEO##############3
TCGA = NULL

for (i in colnames(rt[,3:ncol(rt)])) {
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
  diff=survdiff(Surv(futime, fustat) ~group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if((pValue <0.05) | (i %in% gene)){
    ggsurvplot(fit, 
               data=rt,
               pval=pValue,
               pval.size=5,
               surv.median.line = "hv",
               legend.labs=c("High", "Low"),
               legend.title=i,
               xlab="Time(days)",
               palette=c("#D51317FF", "#0094CDFF"),
               risk.table.height=.15)
    ggsave(paste("E:\\cervix\\",i,".tiff",sep=""),width=6,height=5)
    
    TCGA = c(TCGA, i)}
}
