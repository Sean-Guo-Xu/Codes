library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-CESC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-CESC_SNP.Rdata")
load("TCGA-CESC_SNP.Rdata")
