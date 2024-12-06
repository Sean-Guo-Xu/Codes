setwd("E:\\output")
library(RSQLite)
id = c("1","23")
con <- dbConnect(RSQLite::SQLite(), "E:\\output\\result\\data.db")
dbListTables(con)
cds=dbReadTable(con,"cds")
hmm=dbReadTable(con,"hmm")
bgc=dbReadTable(con,"bgc_features")
hsp=dbReadTable(con,"hsp")
hsp_alignment = dbReadTable(con,"hsp_alignment")

cds = cds[cds$bgc_id %in% id,]
bgc=bgc[bgc$bgc_id %in% id,]
bgc= merge(bgc,cds,by="bgc_id")
hsp = hsp[hsp$cds_id %in% bgc$id,]
hsp_alignment = hsp_alignment[hsp_alignment$hsp_id %in% hsp$id,]
hsp_alignment = merge(hsp_alignment,hsp,by.x = "hsp_id",by.y = "id")
hsp_alignment = hsp_alignment[!duplicated(hsp_alignment$cds_id),]
hsp_alignment = hsp_alignment[,c(5,6,8)]
hsp$link = paste(hsp$cds_id,hsp$hmm_id,sep = "_")
bgc$link = paste(bgc$id,bgc$hmm_id,sep="_")
bgc = bgc[bgc$link %in% hsp$link,]
hmm=hmm[,c(1,2,3,5)]
bgc = merge(bgc,hmm,by.x = "hmm_id",by.y = "id")
bgc= bgc[order(bgc$bgc_id),]
bgc = merge(bgc,hsp_alignment,by.x = "id",by.y ="cds_id",all.x = T )
bgc = bgc[,-c(1,4,12)]
write.csv(bgc,"out.csv",col.names = T,row.names = F)
