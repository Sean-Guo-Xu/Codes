library("tricycle")
library("Seurat")
library(ggplot2)
library(scattermore)
library(scater)
library("ggsci")
library(scales)

setwd("D:\\AAA\\scRNA")

load("5Normal.Rdata")
pbmc = subset(pbmc,cell =="Endothelial cell")
pbmc$EC = rep("Normal EC",length(pbmc$orig.ident))

neurosphere_example<- as.SingleCellExperiment(pbmc,assay = "RNA")
gocc_sce.o <- run_pca_cc_genes(neurosphere_example,
                               exprs_values = "logcounts", species = "human",gname.type = "SYMBOL")

new.ref <- attr(reducedDim(gocc_sce.o, 'PCA'), 'rotation')[, seq_len(2)]
neurosphere_example <- project_cycle_space(neurosphere_example,ref.m = new.ref,gname.type = "SYMBOL",species = "human")
neurosphere_example <- estimate_cycle_position(neurosphere_example, ref.m  = new.ref,
                                               dimred = 'tricycleEmbedding2')

neurosphere_example <- estimate_cycle_position(neurosphere_example,g)
neurosphere_example <- estimate_Schwabe_stage(neurosphere_example,
                                              gname.type = 'SYMBOL',
                                              species = 'human')

nor_neu = neurosphere_example
nor = pbmc

load("ECseurat.Rdata")
neurosphere_example<- as.SingleCellExperiment(pbmc,assay = "RNA")
gocc_sce.o <- run_pca_cc_genes(neurosphere_example,
                               exprs_values = "logcounts", species = "human",gname.type = "SYMBOL")

new.ref <- attr(reducedDim(gocc_sce.o, 'PCA'), 'rotation')[, seq_len(2)]
neurosphere_example <- project_cycle_space(neurosphere_example,ref.m = new.ref,gname.type = "SYMBOL",species = "human")
neurosphere_example <- estimate_cycle_position(neurosphere_example, ref.m  = new.ref,
                                               dimred = 'tricycleEmbedding2')

neurosphere_example <- estimate_cycle_position(neurosphere_example,g)
neurosphere_example <- estimate_Schwabe_stage(neurosphere_example,
                                              gname.type = 'SYMBOL',
                                              species = 'human')
pbmc = merge(nor,pbmc)
neurosphere_example = cbind(nor_neu,neurosphere_example)





plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$EC, 'EC type',
                    bw = 10, fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)+scale_color_npg()
neurosphere_example$sample =c(rep("Normal",6000),rep("AAA",919))

plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$sample, 'EC type',
                    bw = 10, fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)+scale_color_npg()
plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$sample, 'EC type', type = "circular",
                    bw = 10,  fig.title = "Kernel density of \u03b8",line.size=1) +
  theme_bw(base_size = 14)
fit.l <- fit_periodic_loess(neurosphere_example$tricyclePosition,
                            assay(neurosphere_example, 'logcounts')["UGCG",],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(neurosphere_example), ")"))
fit.l$fig + theme_bw(base_size = 14)


see = pbmc$EC
see[see %in% c("EC1","EC2","EC5","EC6","EC7")] = "Others AAA EC"
neurosphere_example$draw = see
plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$draw, 'EC type',
                    bw = 10, fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)+scale_color_npg()


plot_ccposition_den(neurosphere_example$tricyclePosition,
                    neurosphere_example$draw, 'EC type', type = "circular",
                    bw = 10,  fig.title = "Kernel density of \u03b8",line.size=1) +
  theme_bw(base_size = 14)
ggsave("cycle.tiff",width = 6,height = 5)
library(cowplot)
p <- plot_emb_circle_scale(neurosphere_example, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave("EC cycle.tiff",width = 7,height = 5)

EC1 = neurosphere_example[,neurosphere_example$EC =="EC4"]
p <- plot_emb_circle_scale(EC1, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave("EC4 cycle.tiff",width = 7,height = 5)

EC1 = neurosphere_example[,neurosphere_example$EC =="EC3"]
p <- plot_emb_circle_scale(EC1, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
ggsave("EC3 cycle.tiff",width = 7,height = 5)

library(ggpubr)

data = cbind(pbmc$sphingolipid,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Sphingolipid","Cycle","Sample")
data$Sphingolipid = as.numeric(data$Sphingolipid)
data$Cycle = as.numeric(data$Cycle)
ggplot(data, aes(y=Sphingolipid, x=Cycle)) +
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")
ggsave(paste("Sphingolipid","cycle.tiff"),width=5,height = 5)

data = cbind(pbmc$apoptosis,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Apoptosis","Cycle","Sample")
data$Apoptosis = as.numeric(data$Apoptosis)
data$Cycle = as.numeric(data$Cycle)
ggplot(data, aes(y=Apoptosis, x=Cycle)) +
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")+ scale_x_continuous(breaks =seq(0, 2*pi, pi),labels=c('0', expression(pi), expression(2*pi)))
ggsave(paste("Apoptosis","cycle.tiff"),width=5,height = 5)

data = cbind(cds$Pseudotime,neurosphere_example$tricyclePosition,pbmc$sample)
data= as.data.frame(data)
colnames(data) = c("Apoptosis","Cycle","Sample")
data$Apoptosis = as.numeric(data$Apoptosis)
data$Cycle = as.numeric(data$Cycle)
ggplot(data, aes(y=Apoptosis, x=Cycle)) +
  geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")+ scale_x_continuous(breaks =seq(0, 2*pi, pi),labels=c('0', expression(pi), expression(2*pi)))



load(file="ECcds.Rdata")
library(monocle)
for(gene in c("PSAP","SGSM2","UGCG")){
  pbmc$SGSM2 = log2( exprs(pbmc)[gene,]+1)
  data = cbind(pbmc$SGSM2,neurosphere_example$tricyclePosition,pbmc$sample)
  data= as.data.frame(data)
  colnames(data) = c("SGSM2","Position","Sample")
  data$SGSM2 = as.numeric(data$SGSM2)
  data$Position = as.numeric(data$Position)
  ggplot(data, aes(y=SGSM2, x=Position)) +ylab(gene)+xlab("Cycle")+
    geom_point()+theme_bw()+ stat_smooth(formula=y~x,method="loess")+ scale_x_continuous(breaks =seq(0, 2*pi, pi),labels=c('0', expression(pi), expression(2*pi)))
  ggsave(paste(gene,"cycle.tiff"),width=5,height = 5)}