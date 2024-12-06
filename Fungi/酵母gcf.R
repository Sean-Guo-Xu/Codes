class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
ybgc=class[which(class$V7 %in% "Saccharomycetes"),]
elsebgc=class[-which(class$V7 %in% "Saccharomycetes"),]
elsgcf=elsebgc[!duplicated(elsebgc$V4),]
yself=ybgc[-which(ybgc$V4 %in% elsebgc$V4),4]
yself=unique(yself)j
ybgc[!duplicated(ybgc$V4),4]
sor=class[which(class$V7 %in% "Sordariomycetes"),]
eur=class[which(class$V7 %in% "Eurotiomycetes"),]
sandy=sor[which(sor$V4 %in% ybgc$V4),4]
sandy=unique(sandy)
eandy=eur[which(eur$V4 %in% ybgc$V4),4]
eandy=unique(eandy)
write.table(yself,"酵母独有.txt",row.names = F,col.names = F,quote = F)

write.table(sandy,"酵母&粪壳.txt",row.names = F,col.names = F,quote = F)

write.table(without,"没有bgc的基因组.txt",row.names = F,col.names = F,quote = F,sep="\t")
