library(survival)
library(ggplot2)
library(survminer)
library(ExperimentHub)
library(SummarizedExperiment)
options(digits = 4)
setwd("E:\\bladder")
load("TCGA-BLCA.Rdata")
count <- assay(data,"unstranded")
see  =rowData(data)
gene = c("EFEMP1","SVIL")
rownames(count) = see$gene_name
see = colData(data)
time  = cbind(see$barcode,see$vital_status,see$`paper_Combined days to last followup or death`,see$shortLetterCode)
time = time[!(time[,3] %in% c(NA,"[Discrepancy]","[Not Available]")),]
time = time[time[,3] >= 0,]
count = count[,colnames(count) %in% time[,1]]
time = time[order(time[,1]),]
count = count[,order(colnames(count))]
vm = read.table("vm.txt")
count = count[rownames(count) %in% vm$V1,]
count = t(count)
time[time[,2] %in% "Alive",2] = 0
time[time[,2] %in% "Dead",2] = 1
time =as.data.frame(time)
time$V2 = as.numeric(time$V2)
time$V3 = as.numeric(time$V3)
count = cbind(time[,2:3],count)
colnames(count)[1:2] = c("fustat","futime")
rt = count
rt = as.data.frame(rt)
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
ggsave(paste("E:\\bladder\\survival\\TCGA\\",i,".tiff",sep=""),width=6,height=5)

TCGA = c(TCGA, i)}
}
group=cbind(see$barcode,see$gender)
group=group[order(see$barcode),]
group = as.data.frame(group)
group = group[group$V1 %in% rownames(count), ]
group=group$V2
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
ggsurvplot(fit, 
           data=rt,
           pval = pValue,
           surv.median.line = "hv",
           pval.size=5,
           legend.labs=c("F", "M"),
           legend.title="",
           xlab="Time(days)",
           palette=c("#D51317FF", "#0094CDFF"),
           risk.table.height=.15)
ggsave(paste("E:\\bladder\\survival\\TCGA\\","sex",".tiff",sep=""),width=6,height=5)

exp=ifelse(rt[,"EFEMP1"]>median(rt[,"EFEMP1"]),"EFEMP1+","EFEMP1-")
group = paste(group,exp)
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
ggsurvplot(fit, 
           data=rt,
           pval = pValue,
           surv.median.line = "hv",
           pval.size=5,
           legend.labs = c("F EFEMP1-","F EFEMP1+","M EFEMP1-","M EFEMP1+"),
           legend.title="",
           xlab="Time(days)",
           palette=c("#D51317FF", "#0094CDFF", "#EFD500FF" ,"#95C11FFF"),
           risk.table.height=.15)
ggsave(paste("E:\\bladder\\survival\\TCGA\\","sex exp",".tiff",sep=""),width=6,height=5)

##############GEO############
load("GEO-BLCA.Rdata")
geo=count
geo=geo[,c(69:(256-23))]
clinic = read.table("GSE13507clinical.txt",header = T,sep = "\t")
geo = t(geo)
geo = cbind(clinic$overall.survival,clinic$survivalMonth,geo)
geo = as.data.frame(geo)
colnames(geo)[1:2]=c("fustat","futime")
rt = geo
vm = read.table("vm.txt")
rt = rt[,colnames(rt) %in% vm$V1]
rt = cbind(geo[,1:2],rt)
rt = apply(rt,2,as.numeric)
rt=as.data.frame(rt)
rt$futime = rt$futime*365/12
GEO = NULL
for (i in colnames(rt[,3:ncol(rt)])) {
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
  diff=survdiff(Surv(futime, fustat) ~group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if((pValue <0.05) | (i %in% gene)){
    ggsurvplot(fit, 
               data=rt,
               pval = pValue,
               surv.median.line = "hv",
               pval.size=5,
               legend.labs=c("High", "Low"),
               legend.title=i,
               xlab="Time(days)",
               palette=c("#D51317FF", "#0094CDFF"),
               risk.table.height=.15)
    ggsave(paste("E:\\bladder\\survival\\GEO\\",i,".tiff",sep=""),width=6,height=5)
    
      GEO = c(GEO,i)
    }
}
group=clinic$SEX
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
  ggsurvplot(fit, 
             data=rt,
             pval = pValue,
             surv.median.line = "hv",
             pval.size=5,
             legend.labs=c("F", "M"),
             legend.title="",
             xlab="Time(days)",
             palette=c("#D51317FF", "#0094CDFF"),
             risk.table.height=.15)
  ggsave(paste("E:\\bladder\\survival\\GEO\\","sex",".tiff",sep=""),width=6,height=5)

group=ifelse(rt[,"EFEMP1"]>median(rt[,"EFEMP1"]),"EFEMP1+","EFEMP1-")
group = paste(clinic$SEX,group)
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
  ggsurvplot(fit, 
             data=rt,
             pval = pValue,
             surv.median.line = "hv",
             pval.size=5,
              legend.labs = c("F EFEMP1-","F EFEMP1+","M EFEMP1-","M EFEMP1+"),
             legend.title="",
             xlab="Time(days)",
             palette=c("#D51317FF", "#0094CDFF", "#EFD500FF" ,"#95C11FFF"),
             risk.table.height=.15)
  ggsave(paste("E:\\bladder\\survival\\GEO\\","sex exp",".tiff",sep=""),width=6,height=5)

################IM210#################

load("im210.Rdata")
count210 = t(count210)
vm=read.table("vm.txt")
count210 = count210[,colnames(count210) %in% vm$V1]
rt = cbind(clinic210$censOS,clinic210$os,count210)
colnames(rt)[1:2] = c("fustat","futime")
rt = as.data.frame(rt)
rt$futime = rt$futime*365
IM210 = NULL
for (i in colnames(rt[,3:ncol(rt)])) {
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
  diff=survdiff(Surv(futime, fustat) ~group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if((pValue <0.05) | (i %in% gene)){
  ggsurvplot(fit, 
             data=rt,
             surv.median.line = "hv",
             pval.size=5,
             pval =pValue,
             legend.labs=c("High", "Low"),
             legend.title=i,
             xlab="Time(days)",
             palette=c("#D51317FF", "#0094CDFF"),
             risk.table.height=.15)
  ggsave(paste("E:\\bladder\\survival\\IM210\\",i,".tiff",sep=""),width=6,height=5)
  
    IM210 = c(IM210,i)
  }
}
write.table(rt,"IM210exp.txt",col.names = T,row.names = T,sep = "\t",quote = F)
survgene = list(
  TCGA = TCGA,
  GEO=GEO,
  IM210=IM210
)

group=clinic210$Sex
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
ggsurvplot(fit, 
           data=rt,
           pval = pValue,
           surv.median.line = "hv",
           pval.size=5,
           legend.labs=c("F", "M"),
          
           xlab="Time(days)",
           palette=c("#D51317FF", "#0094CDFF"),
           risk.table.height=.15)
ggsave(paste("E:\\bladder\\survival\\IM210\\","sex",".tiff",sep=""),width=6,height=5)

group=ifelse(rt[,"EFEMP1"]>median(rt[,"EFEMP1"]),"EFEMP1+","EFEMP1-")
group = paste(clinic210$Sex,group)
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
ggsurvplot(fit, 
           data=rt,
           pval = pValue,
           surv.median.line = "hv",
           pval.size=5,
           legend.labs = c("F EFEMP1-","F EFEMP1+","M EFEMP1-","M EFEMP1+"),
          legend.title="",
           xlab="Time(days)",
           palette=c("#D51317FF", "#0094CDFF", "#EFD500FF" ,"#95C11FFF"),
           risk.table.height=.15)
ggsave(paste("E:\\bladder\\survival\\IM210\\","sex exp",".tiff",sep=""),width=6,height=5)

rt = rt[!(clinic210$binaryResponse %in% NA),]
group = cbind(as.character(clinic210$Sex),as.character(clinic210$binaryResponse))
group = group[!is.na(group[,2]),]
group = paste(group[,1],group[,2])
fit =  survfit(Surv(futime, fustat) ~ group,data=rt)
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
ggsurvplot(fit, 
           data=rt,
           pval = pValue,
           surv.median.line = "hv",
           pval.size=5,
           legend.labs = c("F CR/PR","F SD/PD","M CR/PR","M SD/PD"),
           legend.title="",
           xlab="Time(days)",
           palette=c("#D51317FF", "#0094CDFF", "#EFD500FF" ,"#95C11FFF"),
           risk.table.height=.15)
ggsave(paste("E:\\bladder\\survival\\IM210\\","sex antiPDL1",".tiff",sep=""),width=6,height=5)



save(survgene,file="survgene.Rdata")

