library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
install.packages("Signac")
library(Signac)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
annotations <- annotations[seqnames(annotations) %in% standardChromosomes(annotations), ]

setwd("E:\\brain\\scATAC")
barcode_file <- "GSM5952324_MB3076_ATAC_barcodes.tsv"
matrix_file <- "GSM5952324_MB3076_ATAC_matrix.mtx"
peaks_file <- "GSM5952324_MB3076_ATAC_peaks.bed"
atac_counts <- Matrix::readMM(file = matrix_file)
barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)
peaks <- read.delim(peaks_file, header = FALSE, stringsAsFactors = FALSE)
rownames(atac_counts) <- paste(peaks$V1, peaks$V2, peaks$V3, sep = "-")
colnames(atac_counts) <- barcodes$V1

colnames(peaks) <- c("chr", "start", "end")
peaks <- makeGRangesFromDataFrame(peaks)
atac<- CreateChromatinAssay(counts = atac_counts, sep = c("-", "-"), genome = 'hg38', ranges = peaks)
atac <- subset(x = atac, subset = nFeature_ATAC > 1000 & nFeature_ATAC < 25000)
genome(annotations) <- "hg38"
Annotation(atac) <- annotations
atac=CreateSeuratObject(atac,assay = "ATAC")
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() 
gene.activities <- GeneActivity(atac, features = VariableFeatures(pbmc))
transfer.anchors <- FindTransferAnchors(reference = pbmc, query = atac, features = VariableFeatures(object = pbmc),
                                        reference.assay = "RNA", query.assay = "ATAC", reduction = "cca")
atac <- FindClusters(atac, resolution = 1, graph.name = 'ATAC_snn')
DimPlot(atac)

pbmc  =Read10X("MB3076")
pbmc = CreateSeuratObject(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc$all = "all"
VlnPlot(pbmc,features = c("percent.mt","percent_hb","nFeature_RNA","percent_hb"),group.by = "all")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <7000 & percent.mt < 10 & percent_hb<5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst")
pbmc=ScaleData(pbmc)
pbmc=RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1, cluster.name = "unintegrated_clusters")
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc), ref = ref, labels = ref$label.main)
pbmc$SingleR = cellpred$labels
DimPlot(pbmc,label = T,group.by = c("seurat_clusters","SingleR"))
DimPlot(pbmc)
marker = FindAllMarkers(pbmc,logfc.threshold = 0.5,min.pct = 0.2)
pbmc$cell = "all"
pbmc$cell[pbmc$seurat_clusters %in% c(0,1,2,3,4,5,7)] = "GP4"
pbmc$cell[pbmc$seurat_clusters %in% c(6)] = "Neuron"
library(Matrix)
atac_counts <- as.matrix(atac_counts)
atac_counts <- as(atac_counts,"dgCMatrix")
atac@assays$ACTIVITY =CreateAssayObject(atac_counts)
DefaultAssay(atac) = "ACTIVITY"
atac  = NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))
transfer.anchors <- FindTransferAnchors(reference = pbmc, query = atac, features = VariableFeatures(object = pbmc), 
                                        reference.assay = "RNA", query.assay = "ATAC", reduction = "cca")
