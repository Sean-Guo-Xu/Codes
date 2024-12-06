setwd("D:\\bigslice\\Tֵ")
library(RSQLite)
library(flexclust)
allscore=NULL
allgcfscore=NULL
allMCscore=NULL
allbestscore = NULL
allARI = NULL
allpur=NULL
for (i in 1:29){
  num = 50+i*50
  file = paste("output",num,"\\result\\data.db",sep="")
sqlite    <- dbDriver("SQLite")
con <- dbConnect(sqlite,file) 
dbListTables(con)
gcf = dbReadTable(con,"gcf_membership")
bgc = dbReadTable(con,"bgc")
gcf  = merge(gcf,bgc,by.x="bgc_id",by.y="id")
mc = read.table("turecluster.txt",sep = "\t")
cdata = merge(gcf,mc,by.x="orig_folder",by.y = "V1")
cdata = cdata[,c(3,12)]
gcfscore = 0
for(i in unique(cdata$gcf_id))
{
  data = cdata[cdata$gcf_id %in% i,]
  score = length(unique(data$V2))-1
  gcfscore = gcfscore+score
}
mcscore = 0
for(i in unique(cdata$V2))
{
  data = cdata[cdata$V2 %in% i,]
  score = length(unique(data$gcf_id))-1
  mcscore = mcscore+score
}
bestscore=0
for (i in unique(cdata$V2)) {
  data = cdata[cdata$V2 %in% i,]
  if(length(unique(data$gcf_id)) ==1){
    data = cdata[cdata$gcf_id %in% unique(data$gcf_id),]
    if(length(unique(data$V2)) ==1){
      bestscore = bestscore+1
    }
  }
}
ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}
pur = ClusterPurity(cdata$V2,cdata$gcf_id)
ari = randIndex(table(cdata$gcf_id,cdata$V2),correct = T)
mcl = length(unique(cdata$V2))
gcfl = length(unique(cdata$gcf_id))
totalscore = mcscore+gcfscore
allscore = c(allscore,totalscore)
allMCscore = c(allMCscore,mcscore)
allgcfscore = c(allgcfscore,gcfscore)
allbestscore = c(allbestscore,bestscore)
allARI = c(allARI,ari)
allpur = c(allpur,pur)
}
vscore = read.table("vscore.txt")
countgcf  =read.table("countgcf.txt")
out = cbind(seq(100,1500,50),vscore,countgcf,allscore,allbestscore,allpur,allARI)
colnames(out) = c("T","vscore","gcfcount","twoscore","besthitscore","purity","ARI")
write.table(out,"result.txt",sep = "\t",quote = F,row.names = F,col.names = T)
