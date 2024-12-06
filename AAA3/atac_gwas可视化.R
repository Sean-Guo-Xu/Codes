setwd("E:\\AAA_gwas")
library(rtracklayer)
library(Gviz)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ggplot2)
library(ComplexHeatmap)
library(biomaRt)
library(dplyr)
library(cowplot)
library(ggbio)
library(showtext)
library(EnsDb.Hsapiens.v86)
library(ggsci)
library(ggrepel)
ensdb <- EnsDb.Hsapiens.v86
######
font_add('Arial','/Library/Fonts/Arial.ttf') #加载字体，MAC 中字体库在 /Library/Fonts
showtext_auto() 
windowsFonts(Times_New_Roman=windowsFont("Times New Roman"))
options(ucscChromosomeNames=FALSE) 
ensembl <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp",GRCh = "38")
ensemblgene <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")##eqtl
######处理gene####
mygenes <- c("CSK","PLEKHJ1","AMH","LDAH","NEK9","PSRC1","SCAPER","HLA-DRB1","FAM66C")  # 示例基因
gene_info <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = mygenes,
  mart = ensemblgene
)
allgene_info = gene_info[gene_info$chromosome_name %in% c(1:22,"x","y"),]
#########ATAC#######
aaa1<- import("ATAC\\AAA1.bw", format="bigWig")
aaa2<- import("ATAC\\AAA2.bw", format="bigWig")
normal<- import("ATAC\\N.bw", format="bigWig")

aaa1_c = read.csv("AAA1__N_diff_peaks_anno.xls",sep = "\t")
aaa2_c = read.csv("AAA2__N_diff_peaks_anno.xls",sep = "\t")



atac_disease1_df <- as.data.frame(aaa1)
atac_disease2_df <- as.data.frame(aaa2)
atac_normal_df <- as.data.frame(normal)
#######eqtl#######
setwd("E:\\SMR")
commond = paste("smr-1.3.1-win.exe --beqtl-summary Artery/Artery_Aorta --query 5e-8 --out eqtl_print")
system(commond,intern = TRUE)
eqtl = read.table("eqtl_print.txt",header=T)
gene_ids <- mapIds(org.Hs.eg.db, keys = eqtl$Probe, column = "SYMBOL", keytype = "ENSEMBL")
eqtl$Gene =gene_ids
eqtl  =eqtl[,c(1,14,10)]

colnames(eqtl) =c("SNP","pvalue","Gene")
snp_info <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start"),
  filters = "snp_filter",
  values = eqtl$SNP,
  mart = ensembl
)
eqtl_data <- merge(eqtl, snp_info, by.x = "SNP", by.y = "refsnp_id")
eqtl_data <- eqtl_data %>% rename(chr = chr_name, pos = chrom_start)
all_eqtl_data  = eqtl_data
save(all_eqtl_data,file="eqtl_5e-8.Rdata")



setwd("E:\\SMR")
commond = paste("smr-1.3.1-win.exe --beqtl-summary Blood/Blood --query 5e-8 --out eqtl_print")
system(commond,intern = T)
eqtl = read.table("eqtl_print.txt",header=TRUE)
eqtl  =eqtl[,c(1,14,10)]
colnames(eqtl) =c("SNP","pvalue","Gene")
snp_info <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start"),
  filters = "snp_filter",
  values = eqtl$SNP,
  mart = ensembl
)

eqtl_blood_data <- merge(eqtl, snp_info, by.x = "SNP", by.y = "refsnp_id")
eqtl_blood_data <- eqtl_blood_data %>% rename(chr = chr_name, pos = chrom_start)
all_eqtl_blood_data  = eqtl_blood_data
gene = mapIds(org.Hs.eg.db,keys= all_eqtl_blood_data$Gene,column = "SYMBOL",keytype = "ENSEMBL")
all_eqtl_blood_data$Gene = gene
all_eqtl_blood_data = all_eqtl_blood_data[!is.na(all_eqtl_blood_data$Gene),]
setwd("E:\\AAA_gwas")
save(all_eqtl_blood_data,file="eqtl_blood_5e-8.Rdata")

load("eqtl_blood_5e-8.Rdata")
load("eqtl_5e-8.Rdata")
##############GWAS##########
setwd("E:\\AAA_gwas")

load("Europe_gwas.Rdata")
genomeAxis <- GenomeAxisTrack()
gwas = gwas[gwas$V7 < 0.0001,c(5,7)]
colnames(gwas) =c("SNP","pvalue")
snp_info <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start"),
  filters = "snp_filter",
  values = gwas$SNP,
  mart = ensembl
)
gwas_data <- merge(gwas, snp_info, by.x = "SNP", by.y = "refsnp_id")
gwas_data <- gwas_data %>% rename(chr = chr_name, pos = chrom_start)
save(gwas_data,file="gwas_1e-4.Rdata")

load("gwas_1e-4.Rdata")
#####################
setwd("E:\\AAA_gwas")
gwas_qtl_col = c("red3","#E18E6D","#62B197","gray")
names(gwas_qtl_col) = c("Both","Only Artery","Only PBMC", "None")
atac_col = c("#E18E6D","#62B197")
peak_col = c("red3","blue3")
biotype_col = ggsci::pal_npg()(9)

for (i in mygenes) {
gene_info = allgene_info[allgene_info$hgnc_symbol %in% i,]
eqtl_data = all_eqtl_data[all_eqtl_data$Gene %in% i, ] 
eqtl_blood_data = all_eqtl_blood_data[all_eqtl_blood_data$Gene %in% i,]
chromosome <- gene_info$chromosome_name
start_position <- gene_info$start_position - 100000
end_position <- max(gene_info$start_position + 100000,gene_info$end_position)

gwas_subset <- gwas_data %>% dplyr::filter(chr == chromosome & pos >= start_position & pos <= end_position)
eqtl_subset <- eqtl_data %>% dplyr::filter(chr == chromosome & pos >= start_position & pos <= end_position)
eqtl_blood_subset <- eqtl_blood_data %>% dplyr::filter(chr == chromosome & pos >= start_position & pos <= end_position)

atac_disease1_subset <- atac_disease1_df %>% dplyr::filter(seqnames == chromosome & start >= start_position & end <= end_position)
atac_disease2_subset <- atac_disease2_df %>% dplyr::filter(seqnames == chromosome & start >= start_position & end <= end_position)
atac_normal_subset <- atac_normal_df %>% dplyr::filter(seqnames == chromosome & start >= start_position & end <= end_position)
aaa1_high_peak = aaa1_c[aaa1_c$GeneName %in% i & aaa1_c$logFC>1,]
aaa1_low_peak = aaa1_c[aaa1_c$GeneName %in% i & aaa1_c$logFC<(-1),]

aaa2_high_peak = aaa2_c[aaa2_c$GeneName %in% i & aaa2_c$logFC>1,]
aaa2_low_peak =  aaa2_c[aaa2_c$GeneName %in% i & aaa2_c$logFC<(-1),]
c_high_peak = rbind(aaa1_high_peak,aaa2_high_peak)
c_low_peak = rbind(aaa1_low_peak,aaa2_low_peak)
artery_common_snps <- intersect(gwas_subset$SNP, eqtl_subset$SNP)
blood_common_snps<-intersect(gwas_subset$SNP,eqtl_blood_subset$SNP)
both_common_snps<-intersect(artery_common_snps,blood_common_snps)
judgesnp<- function(gwas_subset){
  gwas_subset$common = "None"
  gwas_subset$common[gwas_subset$SNP %in% artery_common_snps]="Only Artery"
  gwas_subset$common[gwas_subset$SNP %in% blood_common_snps]="Only PBMC"
  gwas_subset$common[gwas_subset$SNP %in% both_common_snps]="Both"
  return(gwas_subset)
}
gwas_subset = judgesnp(gwas_subset)
eqtl_blood_subset=judgesnp(eqtl_blood_subset)
eqtl_subset=judgesnp(eqtl_subset)
eqtl_blood_subset = eqtl_blood_subset[!duplicated(eqtl_blood_subset$SNP),]
eqtl_subset = eqtl_subset[!duplicated(eqtl_subset$SNP),]

eqtl_blood_subset$pos = as.numeric(eqtl_blood_subset$pos)
eqtl_blood_subset$pvalue[eqtl_blood_subset$pvalue %in% 0] = 1e-300

min_pvalue_snp = gwas_subset[gwas_subset$common %in% "Both",]
min_pvalue_snp = min_pvalue_snp[which.min(min_pvalue_snp$pvalue),]
min_eqtl_blood = eqtl_blood_subset[eqtl_blood_subset$SNP %in% min_pvalue_snp$SNP,]
min_eqtl= eqtl_subset[eqtl_subset$SNP %in% min_pvalue_snp$SNP,]

gwas_plot <- ggplot(gwas_subset, aes(x = pos, y = -log10(pvalue), color = common)) +
  geom_point(size = 1.5) +
  geom_text_repel(data = min_pvalue_snp, aes(label = SNP), size = 3, color = "black")+
  scale_color_manual(values = gwas_qtl_col) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    axis.text.y.left = element_blank(),  # 隐藏左侧坐标文本
    axis.ticks.y.left = element_blank(), # 隐藏左侧刻度
    axis.text.y.right = element_text(size = 12),  # 显示右侧坐标文本
    legend.position = "none",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "AAA GWAS", x = "", y = "-log10(P-value)") +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  scale_y_continuous(name = "-log10(P-value)", sec.axis = sec_axis(~ ., name = NULL))

eqtl_plot <- ggplot(eqtl_subset, aes(x = pos, y = -log10(pvalue), color = common)) +
  geom_point(size = 1.5) +
  geom_text_repel(data = min_eqtl, aes(label = SNP), size = 3, color = "black")+
  scale_color_manual(values = gwas_qtl_col) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    axis.text.y.left = element_blank(),  # 隐藏左侧坐标文本
    axis.ticks.y.left = element_blank(), # 隐藏左侧刻度
    axis.text.y.right = element_text(size = 12),  # 显示右侧坐标文本
    legend.position = "none",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "eQTL - Artery", x = "", y = "-log10(P-value)") +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  scale_y_continuous(name = "-log10(P-value)", sec.axis = sec_axis(~ ., name = NULL))

eqtl_blood_plot <- ggplot(eqtl_blood_subset, aes(x = pos, y = -log10(pvalue), color = common)) +
  geom_point(size = 1.5) +
  geom_text_repel(data = min_eqtl_blood, aes(label = SNP), size = 3, color = "black")+
  scale_color_manual(values = gwas_qtl_col) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    axis.text.y.left = element_blank(),  # 隐藏左侧坐标文本
    axis.ticks.y.left = element_blank(), # 隐藏左侧刻度
    axis.text.y.right = element_text(size = 12),  # 显示右侧坐标文本
    legend.position = "none",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "eQTL - PBMC", x = "", y = "-log10(P-value)") +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  scale_y_continuous(name = "-log10(P-value)", sec.axis = sec_axis(~ ., name = NULL))

maxsignal = max(c(atac_disease1_subset$score,atac_disease2_subset$score,atac_normal_subset$score))
minsignal =  min(c(atac_disease1_subset$score,atac_disease2_subset$score,atac_normal_subset$score))

atac_plot_disease1 <-  ggplot() +
  geom_line(data = atac_disease1_subset, aes(x = start, y = score), color = atac_col[1], size = 0.1) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  geom_rect(data = as.data.frame(aaa1_high_peak), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = peak_col[1], alpha = 1) +
  geom_rect(data = as.data.frame(aaa1_low_peak), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = peak_col[2], alpha = 1) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    axis.text.y.left = element_blank(),  # 隐藏左侧坐标文本
    axis.ticks.y.left = element_blank(), # 隐藏左侧坐标刻度
    axis.text.y.right = element_text(size = 12),  # 显示右侧坐标文本
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "AAA 1", x = "", y = "AAA1") +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  scale_y_continuous(name = "AAA1", sec.axis = sec_axis(~ ., name = NULL),limits = c(minsignal,maxsignal),expand = c(0,0))

atac_plot_disease2 <- ggplot() +
  geom_line(data = atac_disease2_subset, aes(x = start, y = score), color = atac_col[1], size = 0.1) +
  geom_rect(data = as.data.frame(aaa2_high_peak), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = peak_col[1], alpha = 1) +
  geom_rect(data = as.data.frame(aaa2_low_peak), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = peak_col[2], alpha = 1) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    axis.text.y.left = element_blank(),  # 隐藏左侧坐标文本
    axis.ticks.y.left = element_blank(), # 隐藏左侧刻度
    axis.text.y.right = element_text(size = 12),  # 显示右侧坐标文本
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "AAA 2", x = "", y = "AAA2") +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  scale_y_continuous(name = "AAA1", sec.axis = sec_axis(~ ., name = NULL),limits = c(minsignal,maxsignal),expand = c(0,0))

atac_plot_normal <- ggplot() +
  geom_line(data = atac_normal_subset, aes(x = start, y = score), color = atac_col[2], size = 0.1) +
  geom_rect(data = as.data.frame(c_high_peak), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = peak_col[1], alpha = 1) +
  geom_rect(data = as.data.frame(c_low_peak), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = peak_col[2], alpha = 1) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    axis.text.y.left = element_blank(),  # 隐藏左侧坐标文本
    axis.ticks.y.left = element_blank(), # 隐藏左侧刻度
    axis.text.y.right = element_text(size = 12),  # 显示右侧坐标文本
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "Normal", x = "", y = "Control") +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  scale_y_continuous(name = "AAA1", sec.axis = sec_axis(~ ., name = NULL),limits = c(minsignal,maxsignal),expand = c(0,0))

gr <- GRanges(seqnames=chromosome, IRanges(start_position, end_position))
genes <- genes(ensdb, filter = GRangesFilter(gr, "any"))
genes_df <- as.data.frame(genes)
genes_df <- genes_df[genes_df$gene_biotype %in% c("protein_coding","antisense"),]
genes_df$start[genes_df$start<start_position] = start_position
genes_df$end[genes_df$end>end_position] = end_position
genes_df = genes_df[!(genes_df$symbol %in% c("LLfos-48D6.2","AC004490.1","LINGO3","RP11-950C14.7","RP11-950C14.3")),]
gene_axis_plot <- ggplot(genes_df) +
  geom_segment(aes(x = start, xend = end, y = gene_id, yend = gene_id), color = "#E18E6D", size = 2) +
  geom_text(aes(x = (start + end) / 2, y = gene_id, label = gene_name),position = "dodge", color ="#62B197", size = 3, vjust = -0.5) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line.x = element_line(colour = "black", size = 1, lineend = "square")) +
  scale_x_continuous(limits = c(start_position, end_position), expand = c(0, 0)) +
  ylab("Gene Axis")

atac_title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "ATACseq Signal", 
    fontface = 'bold', 
    x = 0.5, 
    hjust = 0.5, 
    size = 16
  )
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    i, 
    fontface = 'bold', 
    x = 0.5, 
    hjust = 0.5, 
    size = 18
  )
final_plot <- cowplot::plot_grid(title,gwas_plot, eqtl_plot,eqtl_blood_plot,atac_title ,atac_plot_disease1,atac_plot_disease2,atac_plot_normal,gene_axis_plot, ncol = 1, rel_heights = c(0.18, 2,2,2,0.1,1,1,1,1.5),align = 'v')

setwd("E:\\AAA_gwas")
pdf(paste(i,"combined_plot.pdf"), width = 9, height = 14)
print(final_plot)
dev.off()
}


gwas_eqtl_leg=Legend(title = "Common SNP", at = names(gwas_qtl_col), 
                     legend_gp = gpar(col =c(gwas_qtl_col)),type = "points",background = "white")
atac_leg=Legend(title = "ATACseq Singal", at = c("AAA","Control"), 
                legend_gp = gpar(fill =atac_col))
peak_leg=Legend(title = "Diff-Peaks", at = c("High in AAA","Low in AAA"), 
                legend_gp = gpar(col =peak_col),type = "lines",background = "white")

lgd_list = packLegend(gwas_eqtl_leg,atac_leg,peak_leg)
draw(lgd_list)

pdf(paste("legend.pdf"),width =5,height = 6)
draw(lgd_list)
dev.off()
