library(immunedeconv)
library(GSVA)
library(ggstatsplot)
library(ggplot2)
library(ggsci)
setwd("E:\\bladder")
gene = "EFEMP1"
load("IM210.Rdata")
vm = read.table("vm.txt")
count=count210
estimate = deconvolute_estimate(count)
vmscore=gsva(count, list(vm$V1), method='ssgsea',  kcdf="Poisson" ,abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
vmscore = normalize(vmscore)
count = log2(count+1)
geneexp = count[gene,]
estimate = rbind(estimate,vmscore,geneexp)
estimate = t(estimate)
estimate = as.data.frame(estimate)
colnames(estimate)[5:6]=c("VMscore","GeneExp")
estimate[,1:6] = apply(estimate[,1:6],2,as.numeric)
clinic = clinic210[,c("binaryResponse","Immune phenotype","Neoantigen burden per MB","TCGA Subtype")]
estimate = cbind(estimate,clinic)

colnames(estimate)[7] ="Type"

col = c("#D51317FF","#0094CDFF" )
names(col) = c("CR/PR" ,"SD/PD")
p1=ggscatterstats(
  data  = estimate,
  y     = StromalScore,
  x     = GeneExp,
  title = "IM210 EFEMP1",
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
ggsave("IM210_gene_point.tiff",width=14,height = 11)

p1=ggscatterstats(
  data  = estimate,
  y     = StromalScore,
  x     = VMscore,
  title = "IM210 Angiogenesis",
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
ggsave("IM210_vm_point.tiff",width=14,height = 11)
colnames(estimate)[9] = "Neoantigen_Burden"
data = estimate[!(estimate$Neoantigen_Burden %in% NA),]
p1=ggscatterstats(
  data  = data,
  y     = Neoantigen_Burden,
  x     = GeneExp,
  title = "IM210 EFEMP1",
)+geom_point(aes(color = data$Type))+guides(color = guide_legend(title = 'Type'))+
  scale_color_manual(values = col)
p2=ggscatterstats(
  data  = data,
  y     = Neoantigen_Burden,
  x     = VMscore,
  title = "IM210 Angiogenesis",
)+geom_point(aes(color = data$Type))+guides(color = guide_legend(title = 'Type'))+
  scale_color_manual(values = col)
p1+p2

ggsave("IM210_Neoantigen_Burden_point.tiff",width = 13,height = 6)
############œ‰Õº#################
col = c("#D51317FF","#0094CDFF" )
names(col) = c("CR/PR" ,"SD/PD")
estimate$Type = as.character(estimate$Type)
data = estimate[!(estimate$Type %in% NA),]
p1=ggbetweenstats(
  data  = data,
  x     = Type,
  y     = GeneExp,
  title = "IM210 EFEMP1",
)+geom_point(aes(color = data$Type))+
  scale_color_manual(values = col)

colnames(estimate)[8] ="Immune_Type"
estimate$Immune_Type=as.character(estimate$Immune_Type)
data2 = estimate[!(estimate$Immune_Type %in% NA),]
col = pal_frontiers()(3)
names(col) = c("desert","excluded","inflamed")
p2=ggbetweenstats(
  data  = data2,
  x     = Immune_Type,
  y     = GeneExp,
)+geom_point(aes(color = data2$Immune_Type))+
  scale_color_manual(values = col)

colnames(estimate)[10] = "TCGA_Stage"
estimate$Immune_Type=as.character(estimate$Immune_Type)
col = pal_frontiers()(4)
names(col) = c("I","II","III","IV")
p3=ggbetweenstats(
  data  = estimate,
  x     = TCGA_Stage,
  y     = GeneExp,
)+geom_point(aes(color = estimate$TCGA_Stage))+
  scale_color_manual(values = col)
p1+p2+p3

col = c("#D51317FF","#0094CDFF" )
names(col) = c("CR/PR" ,"SD/PD")
estimate$Type = as.character(estimate$Type)
data3 = estimate[!(estimate$Type %in% NA),]
p4=ggbetweenstats(
  data  = data3,
  x     = Type,
  y     = VMscore,
  title = "IM210 Angiogenesis",
)+geom_point(aes(color = data3$Type))+
  scale_color_manual(values = col)

colnames(estimate)[8] ="Immune_Type"
estimate$Immune_Type=as.character(estimate$Immune_Type)
data4 = estimate[!(estimate$Immune_Type %in% NA),]
col = pal_frontiers()(3)
names(col) = c("desert","excluded","inflamed")
p5=ggbetweenstats(
  data  = data4,
  x     = Immune_Type,
  y     = GeneExp,
)+geom_point(aes(color = data4$Immune_Type))+
  scale_color_manual(values = col)

colnames(estimate)[10] = "TCGA_Stage"
estimate$Immune_Type=as.character(estimate$Immune_Type)
col = pal_frontiers()(4)
names(col) = c("I","II","III","IV")
p6=ggbetweenstats(
  data  = estimate,
  x     = TCGA_Stage,
  y     = GeneExp,
)+geom_point(aes(color = estimate$TCGA_Stage))+
  scale_color_manual(values = col)
p1+p2+p3+p4+p5+p6
ggsave("IM210_boxplot.tiff",width=18,height = 12)
