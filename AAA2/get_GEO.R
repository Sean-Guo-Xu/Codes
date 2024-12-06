setwd("D:/学习/大创/AAA_AD/转录组学/AAA/GSE7084")
rm(list=ls())

# Package Installation
# BiocManager::install("KEGG.DB",ask = F,update = F)
# BiocManager::install(c("GEOquery","limma","impute" ),
#                      ask = F,update = F)
# BiocManager::install(c("org.Hs.eg.db","hgu133plus2.db"),
#                      ask = F,update = F)

# Data Downloading
library(GEOquery)
eSet <- getGEO('GSE7084', destdir = ".")
save(eSet,file = 'GSE7084_Set.Rdata')

# Data Loading
load('GSE7084_Set.Rdata')
b <- eSet[[1]]
raw_Set <- exprs(b)
phe = pData(b)

# title
library(stringr)
group_list <- str_split(as.character(phe$title),' ',simplify = T)[,1]
group_list[which(group_list=="Abdominal", arr.in=T)] <- "AAA"
save(raw_Set,group_list,
     file = 'GSE7084_raw_Set.Rdata')

write.table(raw_Set, file = 'rawSet.txt', quote = F, sep = "\t")
load(file='GSE7084_raw_Set.Rdata')




