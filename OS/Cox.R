setwd("E:\\bladder\\survival")
library(survival)
library(survminer)
library(ggprism)
library(glmnet)
coxPfilter=0.05          #单因素cox方法显著性过滤标准     #设置工作目录
rt=read.table("TCGAexp.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt[,3:ncol(rt)] = log2(rt[,3:ncol(rt)]+1)
#单因素cox分析
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,i])<0.001){next}
  #单因素cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  #保留显著性基因
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
coxPfilter=0.05          #单因素cox方法显著性过滤标准     #设置工作目录
rt=read.table("GEOexp.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt[,3:ncol(rt)] = log2(rt[,3:ncol(rt)]+1)
rt$futime = rt$futime*365/12
#单因素cox分析
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,i])<0.001){next}
  #单因素cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  #保留显著性基因
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
coxgene$GEO  = outTab[outTab$pvalue<0.05 ,]
write.table(outTab,file="GEO_uni.Cox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="GEO.SigExp.txt",sep="\t",row.names=F,quote=F)
outTab2=outTab

##################IM210##############
coxPfilter=0.05          #单因素cox方法显著性过滤标准     #设置工作目录
rt=read.table("IM210exp.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt[,3:ncol(rt)] = log2(rt[,3:ncol(rt)]+1)
#单因素cox分析
rt$futime = rt$futime*365
sigGenes=c("futime","fustat")
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,i])<0.001){next}
  #单因素cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  #保留显著性基因
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
coxgene$IM210  = outTab[outTab$pvalue<0.05 ,]
save(coxgene,file = "coxgene.Rdata")
write.table(outTab,file="IM210_uni.Cox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="IM210.SigExp.txt",sep="\t",row.names=F,quote=F)
outTab3=outTab
outTab1=as.data.frame(outTab1)
outTab2=as.data.frame(outTab2)
outTab3=as.data.frame(outTab3)
outTab = merge(outTab1,outTab2,by="id")
outTab = merge(outTab,outTab3,by="id")
library(forestploter)
dt=outTab
outTab[,2:13] = apply(outTab[,2:13],2,as.numeric)
dt= outTab[,c(5,9,13)]
dt$TCGA = ""
dt$GEO = ""
dt$IM210 = ""
dt = dt[,c(4,1,5,2,6,3)]
dt$Gene = outTab$id
colnames(dt)[c(2,4,6)] = c("TCGA.PVal","GEO.PVal","IM210.PVal")
colnames(dt)[1] = "            TCGA            "
colnames(dt)[3] = "            GEO            "
colnames(dt)[5] = "            IM210            "
pdf("forest.pdf",width = 12,height = 12)
tm=forest_theme(base_size = 10)
plot=forest(dt[,c(7,1:6)],
            est = list(outTab$HR.x,outTab$HR.y,outTab$HR),
            lower = list(outTab$HR.95L.x,outTab$HR.95L.y,outTab$HR.95L),
            upper = list(outTab$HR.95H.x,outTab$HR.95H.y,outTab$HR.95H),
            ci_column = c(2,4,6),
            arrow_lab = c("Placebo Better", "Treatment Better"),
            ref_line = 1,
            xlim = c(0.6, 1.4),
            ticks_at = c(0.6, 0.8, 1,1.2, 1.4),
            theme = tm)
add_border(plot, row = NULL, col = NULL, part = "header")
dev.off()
library(grid)

plot=edit_plot(plot,row = c(1,7,9),col = 2,which = c("ci"),gp=gpar(col="#D51317FF"))
edit_plot(plot,row = c(2,3,4,5,6,8,10,11),col = 2,which = c("ci"),gp=gpar(col="#0094CDFF"))
