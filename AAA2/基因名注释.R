  smr=read.table("D:\\AAA\\eqtl\\Artery.smr",header = T)
library(clusterProfiler)
library(org.Hs.eg.db)
name <- bitr(smr$probeID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
colnames(name) = c("probeID","Symbol")
smr = merge(smr,name,by="probeID")
write.table(smr,"D:\\AAA\\eqtl\\Artery.smr",col.names = T,row.names = F,sep = "\t",quote = F)
