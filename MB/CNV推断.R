library(biomaRt)
setwd("E:\\brain\\CNV")
Sys.setenv(JAGS_HOME="C:\\Program Files\\JAGS\\JAGS-4.3.1")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
library("rjags")
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 
infercnv::run(infercnv_obj,
              cutoff=1, 
              out_dir="try",
              cluster_by_groups=TRUE, 
              denoise=TRUE,
              HMM=TRUE,
              num_threads=4
)

load("E:\\brain\\4type.Rdata")
gene_list <- pbmc@assays$RNA@counts@Dimnames[[1]] # 示例基因列表
# 获取基因位置数据
gene_positions <- getBM(
  attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
  filters = "external_gene_name",
  values = gene_list,
  mart = ensembl
)
gene_positions = gene_positions[gene_positions$chromosome_name %in% c(1:12,"Y","X"),]
write.table(gene_positions, "gene_order_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cell_names <- colnames(pbmc)
cell_types <- pbmc$cell
annotations <- data.frame(cell_names, cell_types)
write.table(annotations, "annotations_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
expr_matrix <- GetAssayData(pbmc, assay = "RNA", slot = "counts")
write.table(expr_matrix, "expression_matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
