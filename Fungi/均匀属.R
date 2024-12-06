library(ggplot2)

class = read.table("D:\\bigslice\\bgc_information.txt",sep="\t",header=F,fill=T)
phy="Ascomycota"                         
class=class[which(class$V6 %in% phy),]
ord=read.table("D:\\bigslice\\gcf-rank\\Class.txt",header=T,sep="\t")

ggplot(ord) + geom_bar(aes(GCFNumber,reorder(Name,GCFNumber)),stat="identity",fill="#925E9FFF")+ylab("Class")+theme_bw()
ggsave("D:\\bigslice\\pictures\\ALLclass.tiff",height = 6.5,width = 3)

cla=ord[1,1]                        
class=class[which(class$V7 %in% cla),]
ord=read.table("D:\\bigslice\\gcf-rank\\Order.txt",header=T,sep="\t")
ord=ord[ord$Name %in% class$V8,]

ggplot(ord) + geom_bar(aes(GCFNumber,reorder(Name,GCFNumber)),stat="identity",fill="#FDAF91FF")+ylab("Order")+theme_bw()
ggsave("D:\\bigslice\\pictures\\Orders.tiff",height = 3,width = 3)

family=class[which(class$V8 %in% ord[1,1]),]
fam=read.table("D:\\bigslice\\gcf-rank\\Family.txt",header=T,sep="\t")
fam=fam[which(fam$Name %in% family$V9),]

ggplot(fam) + geom_bar(aes(GCFNumber,reorder(Name,GCFNumber)),stat="identity",fill="#42B540FF")+ylab("Family")+theme_bw()
ggsave("D:\\bigslice\\pictures\\families.tiff",height = 3,width = 3)

genus=class[which(class$V9 %in% fam[1,1]),]
gen=read.table("D:\\bigslice\\gcf-rank\\Genus.txt",header=T,sep="\t")
gen=gen[which(gen$Name %in% genus$V10),]

ggplot(gen) + geom_bar(aes(GCFNumber,reorder(Name,GCFNumber)),stat="identity",fill="#0099B4FF")+ylab("Genus")+theme_bw()
ggsave("D:\\bigslice\\pictures\\genus.tiff",height = 4,width = 3)

genus=class[which(class$V10%in% gen[1,1]),]
spe=read.table("D:\\bigslice\\gcf-rank\\Species.txt",header=T,sep="\t")
spe=spe[which(spe$Name %in% genus$V11),]
spe = spe[1:50,]

ggplot(spe) + geom_bar(aes(GCFNumber,reorder(Name,GCFNumber)),stat="identity",fill="#00468BFF")+ylab("Species")+theme_bw()
ggsave("D:\\bigslice\\pictures\\top50_species.tiff",height = 8,width = 4)

genus=class[which(class$V11%in% spe[1,1]),]
str=read.table("D:\\bigslice\\gcf-rank\\Organism.txt",header=T,sep="\t")
str=str[which(str$Name %in% genus$V12),]

ggplot(str) + geom_bar(aes(GCFNumber,reorder(Name,GCFNumber)),stat="identity",fill="#AD002AFF")+ylab("Strain")+theme_bw()
ggsave("D:\\bigslice\\pictures\\strains.tiff",height = 8,width = 4.5)

