setwd("D:\\AAA\\Protein")
library(ggplot2)
aaa1 = read.table("AAA1_N.txt",sep = "\t", fileEncoding = "UTF-8",comment.char = "",fill = T,quote = "",header = T)
aaa1 = as.data.frame(aaa1)
aaa1$show = "Up" 
aaa1$show[aaa1$Ratio<1 ] = "Down" 
table(aaa1$show)

aaa2 = read.table("AAA2_N.txt",sep = "\t", fileEncoding = "UTF-8",comment.char = "",fill = T,quote = "",header = T)
aaa2 = as.data.frame(aaa2)
aaa2$show = "Up" 
aaa2$show[aaa2$Ratio<1 ] = "Down" 
table(aaa2$show)
data  = cbind(c(table(aaa1$show),table(aaa2$show)),c("Up","Down","Up","Down"),c("AAA1","AAA1","AAA2","AAA2"))
data = as.data.frame(data)
data$Num = as.numeric(data$Num)
windowsFonts(Times_New_Roman=windowsFont("Arial"))
colnames(data) = c("Num","Expression","Sample")
ggplot(data,aes(x =Sample,  y = Num,fill = Expression))+geom_bar(stat="identity",width=0.9, position='dodge')+
  geom_text(stat="identity",aes(label=Num), color="black", size=3.5,position=position_dodge(0.9),vjust=-0.5)+
  scale_fill_manual(values =c ( "#6BB7CA",  "#E07B54"))+theme_classic() +xlab("")+ylab("Number of Diff-Protein")+theme(text = element_text(family='Times_New_Roman'))
ggsave("Num_bar.png",width = 4,height = 5)
library(biomaRt)
library(ggvenn)
a = list(
    AAA1_Up = aaa1$Accession[aaa1$show %in% "Up"] ,
    AAA1_Down = aaa1$Accession[aaa1$show %in% "Down"],
    AAA2_Up = aaa2$Accession[aaa2$show %in% "Up"],
    AAA2_Down = aaa2$Accession[aaa2$show %in% "Down"]
)
pdf("protein_venn.pdf",width = 8,height = 6)
ggvenn(a,show_percentage = F,fill_color = c("#E07B54","#6BB7CA","#E07B54","#6BB7CA"),fill_alpha = 0.8,set_name_size = 5,text_size = 5)
dev.off()
library(ensembldb)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

###volcano
samepro= intersect(aaa1$Description,aaa2$Description)
pro = merge(aaa1,aaa2,by="Description")
pro$Exp= paste(pro$show.x,pro$show.y)
pro = cbind(pro$Ratio.x,pro$Ratio.y,pro$Exp)
pro = as.data.frame(pro)
pro$V1 = as.numeric(pro$V1)
pro$V2 = as.numeric(pro$V2)
colnames(pro) = c("Ratio_1","Ratio_2","Exp")
write.csv(pro,"volcano.csv",row.names = F)

max_limit <- max(max(pro$Ratio_1, na.rm = TRUE), max(pro$Ratio_2, na.rm = TRUE))
min_limit <- min(min(pro$Ratio_1, na.rm = TRUE), min(pro$Ratio_2, na.rm = TRUE))

ggplot(pro, aes(x = Ratio_1, y = Ratio_2, color = Exp)) +
  geom_point(size = 2, alpha = 0.7) +  # Adjust size and transparency
  geom_vline(xintercept = c(1/1.2, 1.2), linetype = "dashed", color = "black") +  # Vertical lines at x = 0.8 and x = 1.2
  geom_hline(yintercept = c(1/1.2, 1.2), linetype = "dashed", color = "black") +  # Horizontal lines at y = 0.8 and y = 1.2
  scale_color_manual(values = c("Up Up" = "#E07B54", "Down Down" = "#6BB7CA", "Up Down" = "#BC80BD", "Down Up" = "#BC80BD", "Other" = "grey")) +  # Custom colors
  labs(
    x = "Fold Change 1 (AAA1 vs N)",
    y = "Fold Change 2 (AAA2 vs N)",
    title = "Scatter Plot of Fold Change"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and bold the title
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12),  # Increase axis text size
    legend.position = "right",  # Move legend to the right
    legend.title = element_blank(),  # Remove legend title
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA)  # Add border
  ) +
  annotate("text", x = 1.313, y = 1.247, label = "HLA class II DRα", 
           hjust = 0, vjust = -1.5, size = 3, color = "black") +  # Add the label text
  annotate("segment", x = 1.313, xend = 1.313, y = 1.247, yend = 1.247, 
           arrow = arrow(length = unit(0.03, "npc"), type = "closed"), color = "black") +  # Add a segment (arrow) pointing to the point
  scale_x_continuous(trans = 'log10', limits = c(min_limit, max_limit)) +  # Set x-axis to log scale and consistent limits
  scale_y_continuous(trans = 'log10', limits = c(min_limit, max_limit))    # Set y-axis to log scale and consistent limits
ggsave("protein_volcano.png",width = 6,height = 5)
library(biomaRt)
# 获取蛋白质的注释信息

####相互作用###
library(STRINGdb)
samepro= intersect(aaa1$Description,aaa2$Description)
pro = merge(aaa1,aaa2,by="Description")
pro$Exp= paste(pro$show.x,pro$show.y)
pro = pro[pro$Exp %in% c("Down Down","Up Up"),]
uniprot_ids <- pro$Accession.x
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")
ppi_data <- string_db$get_interactions(unique(annot$external_gene_name))
mapped_data <- string_db$map(data.frame(uniprot_ids), "uniprot_ids", removeUnmappedRows = TRUE)
library(ggraph)
library(igraph)
string_ids <- mapped_data$STRING_id
ppi_data <- string_db$get_interactions(string_ids)
ppi_data <- ppi_data[ppi_data$from %in% "9606.ENSP00000378786" | ppi_data$to %in% "9606.ENSP00000378786",]

colnames(mapped_data)[2]  = "from"
proteins = merge(ppi_data,mapped_data,by="from")
colnames(mapped_data)[2]  = "to"
proteins = merge(proteins,mapped_data,by="to")
pro = pro[pro$Accession.x %in% c(proteins$uniprot_ids.x,proteins$uniprot_ids.y),]
# in "FGB"  "C1QB" "CTSB" "C1QC" "C1QA"
# 
colnames(pro)[c(2,15)] = colnames(proteins)[4:5]
proteins = merge(proteins,pro[,c(1,2)],by=colnames(proteins)[4])
proteins = merge(proteins,pro[,c(1,15)],by=colnames(proteins)[5])
proteins = merge()
ppi_data$to_uniprot <- string_db$get_uniprot(ppi_data$to)
graph <- graph_from_data_frame(ppi_data, directed = FALSE)
ggraph(graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = ppi_data$combined_score), arrow = arrow(length = unit(3, 'mm')), show.legend = FALSE) +
  geom_node_point(color = "steelblue", size = 5) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  theme_void()
ppi_data 
