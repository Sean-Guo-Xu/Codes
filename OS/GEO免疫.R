setwd("E:\\bladder\\infiltration")
library(ggstatsplot)
library(ggplot2)
load("E:\\bladder\\GEO-BLCA.Rdata")
vm = read.table("E:\\bladder\\vm.txt")  
library(immunedeconv)
library(ggsci)
gene = "EFEMP1"
name = rownames(count)
count = apply(count,2,as.numeric)
rownames(count)=name
estimate = deconvolute_estimate(count)
geneexp = log2(count["EFEMP1",]+1)
vmexp = count[rownames(count) %in% vm$V1,]
estimate = rbind(estimate,geneexp)
estimate = t(estimate)
estimate = as.data.frame(estimate)
colnames(estimate)[5] = "GeneExp"
vmscore=gsva(count, list(vm$V1), method='ssgsea', kcdf='Poisson', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore = normalize(vmscore)
clinicgeo=as.data.frame(clinicgeo)
vmscore = t(vmscore)
estimate = cbind(estimate,vmscore,clinicgeo$geotype)
col = c("#0094CDFF","#95C11FFF", "#D51317FF", "#F39200FF" )
names(col) = unique(estimate$`clinicgeo$geotype`)
estimate$`clinicgeo$geotype` = factor(estimate$`clinicgeo$geotype`,levels = unique(estimate$`clinicgeo$geotype`))
p1=ggscatterstats(
  data  = estimate,
  y     = StromalScore,
  x     = GeneExp,
  title = "GEO EFEMP1",
)+geom_point(aes(color = estimate$`clinicgeo$geotype`))+guides(color = guide_legend(title = 'Type'))+
  scale_color_manual(values = col)
p2=ggscatterstats(
  data  = estimate,
  y     = ImmuneScore,
  x     = GeneExp,
)+geom_point(aes(color = estimate$`clinicgeo$geotype`))+theme(legend.position = 'none')+
  scale_color_manual(values = col)
p3=ggscatterstats(
  data  = estimate,
  y     = ESTIMATEScore,
  x     = GeneExp,
)+geom_point(aes(color = estimate$`clinicgeo$geotype`))+theme(legend.position = 'none')+
  scale_color_manual(values = col)
p4=ggscatterstats(
  data  = estimate,
  y     = TumorPurity,
  x     = GeneExp,
)+geom_point(aes(color = estimate$`clinicgeo$geotype`))+theme(legend.position = 'none')+
  scale_color_manual(values = col)
p1+p2+p3+p4
ggsave("GEO_gene_point.tiff",width=14,height = 11)
 p5=ggscatterstats(
   data  = estimate,
   y     = StromalScore,
   x     = GeneExp,
   title = "GEO Angiogenesis",
 )+geom_point(aes(color = estimate$`clinicgeo$geotype`))+guides(color = guide_legend(title = 'Type'))+
   scale_color_manual(values = col)
 p6=ggscatterstats(
   data  = estimate,
   y     = ImmuneScore,
   x     = GeneExp,
 )+geom_point(aes(color = estimate$`clinicgeo$geotype`))+theme(legend.position = 'none')+
   scale_color_manual(values = col)
 p7=ggscatterstats(
   data  = estimate,
   y     = ESTIMATEScore,
   x     = GeneExp,
 )+geom_point(aes(color = estimate$`clinicgeo$geotype`))+theme(legend.position = 'none')+
   scale_color_manual(values = col)
 p8=ggscatterstats(
   data  = estimate,
   y     = TumorPurity,
   x     = GeneExp,
 )+geom_point(aes(color = estimate$`clinicgeo$geotype`))+theme(legend.position = 'none')+
   scale_color_manual(values = col)
 p5+p6+p7+p8
 ggsave("GEO_vm_point.tiff",width=14,height = 11)

 
 
 colnames(estimate)[7] ="Type"
 p1=ggbetweenstats(
   data  = estimate,
   x     = Type,
   y     = StromalScore,
   subtitle = F,
   title="GEO",
 )+geom_point(aes(color = estimate$Type))+
   scale_color_manual(values = col)
 
 p2=ggbetweenstats(
   data  = estimate,
   x     = Type,
   y     = ImmuneScore,
   subtitle = F,
 )+geom_point(aes(color = estimate$Type))+
   scale_color_manual(values = col)
 
 p3=ggbetweenstats(
   data  = estimate,
   x     = Type,
   y     = TumorPurity,
   subtitle = F,
 )+geom_point(aes(color = estimate$Type))+
   scale_color_manual(values = col)
 
 p4=ggbetweenstats(
   data  = estimate,
   x     = Type,
   y     = ESTIMATEScore,
   subtitle = F,
 )+geom_point(aes(color = estimate$Type))+
   scale_color_manual(values = col)
 
 p5=ggbetweenstats(
   data  = estimate,
   x     = Type,
   y     = vmscore,
   subtitle = F,
 )+geom_point(aes(color = estimate$Type))+
   scale_color_manual(values = col)
 
 p6=ggbetweenstats(
   data  = estimate,
   x     = Type,
   y     = GeneExp,
   subtitle = F,
 )+geom_point(aes(color = estimate$Type))+
   scale_color_manual(values = col)
p1+p2+p3+p4+p5+p6
ggsave("GEO_boxplot.tiff",width = 18,height = 10)
library(limma)
exp = read.table("GEO_DDR.txt",sep = "\t",header = T)
name = exp$SYMBOL
exp=exp[,-1]
exp=as.matrix(exp)
rownames(exp) = name
exp = avereps(exp)
vmscore=gsva(exp, list(vm$V1), method='ssgsea', kcdf='Poisson', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore = normalize(vmscore)
exp = rbind(log2(exp[gene,]+1),vmscore)
data = NULL
for (i in colnames(exp)) {
  data = rbind(data,c(i,exp[2,i],exp[1,i])) 
  
}
data = data[order(data[,1]),]
data = as.data.frame(data)

col = c("#0094CDFF","#D51317FF" )
names(col) = c("DDR","Untreated")
data$Type  = c("DDR","DDR","Untreated","DDR","Untreated","DDR","Untreated","DDR","Untreated")
colnames(data)[1:3] = c("Sample","VMscore",paste(gene,"Exp",sep = ""))
data$VMscore = as.numeric(data$VMscore)
data$EFEMP1Exp = as.numeric(data$EFEMP1Exp)
p1=ggbetweenstats(
  data  = data,
  x     = Type,
  y     = VMscore,
  title = "VMscore under DDR treatment"
)+geom_point(aes(color = data$Type))+
  scale_color_manual(values = col)

p2=ggbetweenstats(
  data  = data,
  x     = Type,
  y     = EFEMP1Exp,
  title = "EFEMP1Exp under DDR treatment"
)+geom_point(aes(color = data$Type))+
  scale_color_manual(values = col)
p1+p2
ggsave("DDR_boxplot.tiff",width = 8,height = 5)

