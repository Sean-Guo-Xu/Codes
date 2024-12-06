library(rjson)
setwd("D:\\bigslice\\化合物\\mibig_json_2.0")
name=list.files ()
cla=c()
org=c()




for (i in name){
  jsonData <- fromJSON(file =i)
  if (length(jsonData$cluster$biosyn_class)==1){
  cla=c(cla,jsonData$cluster$biosyn_class)}
  else{
    mucla=""
    for (j in jsonData$cluster$biosyn_class){
      mucla=paste(mucla,j,sep = " ")
    }
    mucla=substr(mucla,2,nchar(mucla))
    cla=c(cla,mucla)
  }
  org=c(org,jsonData$cluster$organism_name)

}
ta=cbind(name,cla,org)
ta[,1]=substr(ta[,1],1,10)


setwd("D:\\bigslice\\化合物")
bgc_gcf=read.table("bgc_id-gcf关系.txt",header=T)
bgc_mib=read.table("bgc_id-mibig关系.txt",header=T)
see=merge(bgc_gcf,bgc_mib,by="bgc_id")
see$name=substr(see$name,1,10)
see=merge(see,ta,by="name")
write.table(see,"GCFannotation.xls",row.names = F,sep="\t")
