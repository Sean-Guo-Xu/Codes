setwd("D:\\AAA\\omics\\beta")
id="84703"
up=read.table(paste(id,"_uptarget_associate_peaks.bed",sep=""),header = T)
down=read.table(paste(id,"_downtarget_associate_peaks.bed",sep = ""),header=T)
up=up[order(up$Score,decreasing = T),]
up=up[!duplicated(up$Symbol),]
down=down[order(down$Score,decreasing = T),]
down=down[!duplicated(down$Symbol),]
down$Score=-down$Score
up1 = up
down1 = down


id="88612"
up=read.table(paste(id,"_uptarget_associate_peaks.bed",sep=""),header = T)
down=read.table(paste(id,"_downtarget_associate_peaks.bed",sep = ""),header=T)
up=up[order(up$Score,decreasing = T),]
up=up[!duplicated(up$Symbol),]
down=down[order(down$Score,decreasing = T),]
down=down[!duplicated(down$Symbol),]
down$Score=-down$Score
up2=up
down2=down
up=merge(up1,up2,by="Symbol")
down=merge(down1,down2,by="Symbol")
all=rbind(up,down)
write.table(cbind(all$Symbol,all$Score.x,all$Score.y),"mmH3K_chipseq_omics.txt",quote = F,sep="\t",row.names = F,col.names  =F)


id="EP300AAA"
up=read.table(paste(id,"_uptarget_associate_peaks.bed",sep=""),header = T)
down=read.table(paste(id,"_downtarget_associate_peaks.bed",sep = ""),header=T)
up=up[order(up$Score,decreasing = T),]
up=up[!duplicated(up$Symbol),]
down=down[order(down$Score,decreasing = T),]
down=down[!duplicated(down$Symbol),]
down$Score=-down$Score
up1=up
down1=down

id="104697"
up=read.table(paste(id,"_uptarget_associate_peaks.bed",sep=""),header = T)
down=read.table(paste(id,"_downtarget_associate_peaks.bed",sep = ""),header=T)
up=up[order(up$Score,decreasing = T),]
up=up[!duplicated(up$Symbol),]
down=down[order(down$Score,decreasing = T),]
down=down[!duplicated(down$Symbol),]
down$Score=-down$Score
up2=up
down2=down

id="104699"
up=read.table(paste(id,"_uptarget_associate_peaks.bed",sep=""),header = T)
down=read.table(paste(id,"_downtarget_associate_peaks.bed",sep = ""),header=T)
up=up[order(up$Score,decreasing = T),]
up=up[!duplicated(up$Symbol),]
down=down[order(down$Score,decreasing = T),]
down=down[!duplicated(down$Symbol),]
down$Score=-down$Score
up3=up
down3=down
up=merge(up1,up2,by="Symbol")
up=merge(up,up3,by="Symbol")
down=merge(down1,down2,by="Symbol")
down=merge(down,down3,by="Symbol")
all=rbind(up,down)
write.table(cbind(all$Symbol,all$Score.x,all$Score.y,all$Score),"AAAH3K_chipseq_omics.txt",quote = F,sep="\t",row.names = F,col.names  =F)


