library(immunedeconv)
library(GSVA)
library(ggstatsplot)
library(ggplot2)
library(ggsci)
setwd("E:\\bladder")
gene = "EFEMP1"
load("TCGA_tpm.Rdata")
clinicdata=as.data.frame(clinicdata)
clinic = cbind(unlist(clinicdata$sample_type),unlist(clinicdata$ajcc_pathologic_stage),unlist(clinicdata$ajcc_pathologic_t),unlist(clinicdata$ajcc_pathologic_n),unlist(clinicdata$ajcc_pathologic_m))
colnames(clinic) = c("Type","Stage","T","N","M")
vm = read.table("vm.txt")
estimate = deconvolute_estimate(count)
count = log2(count+1)
vmscore=gsva(count, list(vm$V1), method='ssgsea',  kcdf="Gaussian" ,abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore = normalize(vmscore)
geneexp = count[gene,]
estimate = rbind(estimate,vmscore,geneexp)
rownames(clinic) = colnames(estimate)
estimate = rbind(estimate,t(clinic))
estimate = t(estimate)
estimate = as.data.frame(estimate)
colnames(estimate)[5:6]=c("VMscore","GeneExp")
estimate[,1:6] = apply(estimate[,1:6],2,as.numeric)
col = c("#D51317FF","#0094CDFF" )
names(col) = c("Primary Tumor" ,"Solid Tissue Normal")
p1=ggscatterstats(
  data  = estimate,
  y     = StromalScore,
  x     = GeneExp,
  title = "TCGA EFEMP1",
)+geom_point(aes(color = estimate$Type))+guides(color = guide_legend(title = 'Type'))+
  scale_color_manual(values = col)

p2=ggscatterstats(
  data  = estimate,
  y     = ImmuneScore,
  x     = GeneExp,
)+geom_point(aes(color = estimate$Type))+theme(legend.position = 'none')+
  scale_color_manual(values = col)

p3=ggscatterstats(
  data  = estimate,
  y     = ESTIMATEScore,
  x     = GeneExp,
)+geom_point(aes(color = estimate$Type))+theme(legend.position = 'none')+
  scale_color_manual(values = col)

p4=ggscatterstats(
  data  = estimate,
  y     = TumorPurity,
  x     = GeneExp,
)+geom_point(aes(color = estimate$Type))+theme(legend.position = 'none')+
  scale_color_manual(values = col)
p1+p2+p3+p4
ggsave("TCGA_gene_point.tiff",width=14,height = 11)
p1=ggscatterstats(
  data  = estimate,
  y     = StromalScore,
  x     = VMscore,
  title = "TCGA Angiogenesis",
)+geom_point(aes(color = estimate$Type))+guides(color = guide_legend(title = 'Type'))+
  scale_color_manual(values = col)

p2=ggscatterstats(
  data  = estimate,
  y     = ImmuneScore,
  x     = VMscore,
)+geom_point(aes(color = estimate$Type))+theme(legend.position = 'none')+
  scale_color_manual(values = col)

p3=ggscatterstats(
  data  = estimate,
  y     = ESTIMATEScore,
  x     = VMscore,
)+geom_point(aes(color = estimate$Type))+theme(legend.position = 'none')+
  scale_color_manual(values = col)

p4=ggscatterstats(
  data  = estimate,
  y     = TumorPurity,
  x     = VMscore,
)+geom_point(aes(color = estimate$Type))+theme(legend.position = 'none')+
  scale_color_manual(values = col)
p1+p2+p3+p4
ggsave("TCGA_vm_point.tiff",width=14,height = 11)

#######œ‰œﬂÕº######
col = c("#0094CDFF" ,"#D51317FF")
names(col) = c("Solid Tissue Normal","Primary Tumor" )
p1=ggbetweenstats(
  data  = estimate,
  x     = Type,
  y     = GeneExp,
  title = "TCGA EFEMP1",
)+geom_point(aes(color = estimate$Type))+
  scale_color_manual(values = col)

col = pal_frontiers()(4)
names(col) = c("Stage I","Stage II"  ,"Stage III" ,"Stage IV"  )
data = estimate[!(estimate$Stage %in% NA),]
p2=ggbetweenstats(
  data  = data,
  x     = Stage,
  y     = GeneExp,
  subtitle = F,
)+geom_point(aes(color = data$Stage))+
  scale_color_manual(values = col)

col = c("#0094CDFF" ,"#D51317FF")
names(col) = c("Solid Tissue Normal","Primary Tumor" )
p3=ggbetweenstats(
  data  = estimate,
  x     = Type,
  y     = VMscore,
  title = "TCGA Angiogenesis",
)+geom_point(aes(color = estimate$Type))+
  scale_color_manual(values = col)

col = pal_frontiers()(4)
names(col) = c("Stage I","Stage II"  ,"Stage III" ,"Stage IV"  )
data = estimate[!(estimate$Stage %in% NA),]
p4=ggbetweenstats(
  data  = data,
  x     = Stage,
  y     = VMscore,
  subtitle = F,
)+geom_point(aes(color = data$Stage))+
  scale_color_manual(values = col)

p1+p2+p3+p4
ggsave("TCGA_vm_box.tiff",width=14,height = 10)
