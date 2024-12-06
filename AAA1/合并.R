library(Seurat)
setwd("E:\\AAA_scRNA")
part=read.csv("H2_expression_counts.csv",row.names = 1)
part2=read.csv("H3_expression_counts.csv",row.names = 1)
part3=read.csv("H4_expression_counts.csv",row.names = 1)
part4=read.csv("H5_expression_counts.csv",row.names = 1)
part5=read.csv("H6_expression_counts.csv",row.names = 1)
part6=read.csv("H7_expression_counts.csv",row.names = 1)
part  =cbind(part,part2,part3,part4,part5,part6)
pbmc <- CreateSeuratObject(counts = part,min.cells = 3)
pbmc$sample = c(rep("Blood1",ncol(part)),rep("Blood2",ncol(part2)),rep("Blood3",ncol(part3)),rep("Blood4",ncol(part4)),rep("Blood5",ncol(part5)),rep("Blood6",ncol(part6)))
save(pbmc,file="N_blood.Rdata")
load("N_blood.Rdata")
pbmc$cell = allcell$cell_type
pbmc$index = rownames(pbmc@meta.data)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200& nFeature_RNA <5500 & percent.mt < 15 & percent_hb<5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.8)                  #对细胞分组,优化标准模块化
pbmc <- RunUMAP(pbmc,dims = 1:20)
DimPlot(pbmc,label=T,group.by = c("sample","seurat_clusters"))
library(DoubletFinder)
pbmc_list <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
pbmc_stats <- summarizeSweep(pbmc_list, GT = FALSE)
bcmvn <- find.pK(pbmc_stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] 
pK_bcmvn = as.character(pK_bcmvn) 
pK_bcmvn = as.numeric(pK_bcmvn)
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)    
nExp_poi <- round(0.075*nrow(pbmc@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
memory.limit(size=1000000)
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.3_3272")
pbmc = subset(pbmc,DF.classifications_0.25_0.3_3272 == "Singlet" )
DimPlot(pbmc,label = T,group.by = "DF.classifications_0.25_0.3_3272")
save(pbmc,file = "N_blood.Rdata")
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE)
############手动注释#########
library(Seurat)
library(SingleR)
load("N_marker.Rdata")
load("N_blood.Rdata")
DimPlot(pbmc,group.by = c("seurat_clusters","cell","singleR"),label=T)
data = pbmc@assays$RNA@data
data= as.matrix(data)
for( i in 1:8){
  write.csv(data[,(5000*(i-1)+1):(5000*i)],paste("part",i,".csv",sep = ""),col.names = T,row.names = T)
}
write.csv(data[,((5000*i)+1):41678],paste("part",9,".csv",sep = ""),col.names = T,row.names = T)
cell = NULL
for(i in 1:9){
  part= read.csv(paste("human_Blood_part",i,".csv",sep=""))
  cell = rbind(cell,part)
}
pbmc$cell = cell$cell_type
library(SingleR)
library(celldex)
ref=celldex::BlueprintEncodeData()
cellpred <- SingleR(test =GetAssayData(pbmc, "data"), ref = ref, labels = ref$label.fine)
pbmc$singleR = cellpred$labels

label=cbind(pbmc$cell,pbmc$seurat_clusters)
label[label[,2] %in% c(48),1] = "HSC"
label[label[,2] %in% c(42,7,9),1] = "B cell" 
label[label[,2] %in% c(39),1]="T cell" 
label[label[,2] %in% c(41,21,16,34,29,28,49,31),1]="Monocyte" 
label[label[,2] %in% c(26,47),1]="Neutrophil"
label[label[,2] %in% c(43),1]= "Plasma cell"
label[label[,2] %in% c(46),1]="Unsure" 
label[label[,2] %in% c(24,32),1]="Megakaryocyte"
pbmc$cell = label[,1]


library(Seurat)
setwd("E:\\AAA_scRNA\\blood")
load("N_blood.Rdata")
con = pbmc
load("A_blood.Rdata")
pbmc$aaa_ro_con= rep("AAA",length(pbmc$orig.ident))
con$aaa_ro_con = rep("Normal",length(con$orig.ident))
pbmc = merge(pbmc,con)
unique(pbmc$cell)
pbmc  = subset(pbmc, cell %in% c("T cell", "B cell", "Monocyte","Neutrophil" ,"Megakaryocyte","Plasma cell"))
markerlist = list()
for (i in unique(pbmc$cell)) {
  part = subset(pbmc,cell== i)
  marker = FindMarkers(part,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "aaa_ro_con",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
  markerlist[[i]] = marker 
}
marker = FindMarkers(pbmc,assay="RNA",ident.1="AAA",ident.2="Normal",group.by = "aaa_ro_con",  only.pos = FALSE,logfc.threshold=0,min.pct=0)
markerlist[["all"]] = marker 
diff = list()
for (i in names(markerlist)){
  diff[[i]] = cbind(markerlist[[i]],rownames(markerlist[[i]]))
}
save(diff,file="blood_diff_all.Rdata")

gene = pbmc@assays$RNA@counts@Dimnames[[1]]
gene = gene[gene %in% con@assays$RNA@counts@Dimnames[[1]]]
for (i in names(bdiff)) {
   part = bdiff[[i]] 
   part = part[part$`rownames(markerlist[[i]])` %in% gene,]
   bdiff[[i]] = part
}
save(bdiff,file="E:\\AAA_scRNA\\Blood\\blood_diff_all.Rdata")
