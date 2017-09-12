
gn <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_gene.txt")

names(gn)[c(1, 22)] <- c("Peak.ID", "Smad4.Motif")
#substitute character after first Space with
#null A space (), then any character (.) any number of times (*) 
#until the end of the string ($). See ?regex to learn regular expressions. 
attach(gn)
gn$Annotation <- sub(" .N.*$", "", Annotation)
detach(gn)
attach(gn)
gn$Annotation <- sub("intron", "Intron", Annotation)
detach(gn)
attach(gn)
gn$Annotation <- sub("exon", "Exon", Annotation)
detach(gn)
attach(gn)
gn$Annotation <- sub("3' UTR", "Intergenic", Annotation)
detach(gn)
attach(gn)
gn$Annotation <- sub("5' UTR", "Intergenic", Annotation)
detach(gn)
attach(gn)
gn$Annotation <- sub("non-coding", "Intergenic", Annotation)
detach(gn)
attach(gn)
gn$Annotation <- sub("TTS", "Intergenic", Annotation)
detach(gn)
attach(gn)
gn$Annotation[Distance.to.TSS >= -5000 & Distance.to.TSS <= 5000] <- "Promoter (¡À5k)"
detach(gn)

#gn_02k <- subset(gn, Distance.to.TSS >= -2000 & Distance.to.TSS <= 2000, 
#                     select=c(2, 3, 4, 6, 8, 10, 16, 19, 22))
#gn_05k <- subset(gn, Distance.to.TSS >= -5000 & Distance.to.TSS <= 5000, 
#                 select=c(2, 3, 4, 6, 8, 10, 16, 19, 22))
#gn_10k <- subset(gn, Distance.to.TSS >= -10000 & Distance.to.TSS <= 10000, 
#                 select=c(2, 3, 4, 6, 8, 10, 16, 19, 22))

gn_05k <- unique(gn[gn$Distance.to.TSS >= -5000 & gn$Distance.to.TSS <= 5000, ]$Gene.Name)
gn_10k <- unique(gn[gn$Distance.to.TSS >= -10000 & gn$Distance.to.TSS <= 10000, ]$Gene.Name)



#gn_05k
write.table(gn_05k, "gn_05k.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_10k, "gn_10k.txt", row.names = F, col.names = F, quote = F, eol = " ")
#demo(plotmath)



rna12.5 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E12.5_M_E12.5_W/DEGseq/geneExp_E12.5_M_E12.5_W_filtered_DEGseq.txt")
rna13.5 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E13.5_M_E13.5_W/DEGseq/geneExp_E13.5_M_E13.5_W_filtered_DEGseq.txt")



#rna_e12.5_down <- rna_e12.5_down[!duplicated(rna_e12.5_down$Symbol),]
#rna_e12.5_up <- rna_e12.5_up[!duplicated(rna_e12.5_up$Symbol),]
#rna_e13.5_down <- rna_e13.5_down[!duplicated(rna_e13.5_down$Symbol),]
#rna_e13.5_up <- rna_e13.5_up[!duplicated(rna_e13.5_up$Symbol),]

#rna12.5 <- rna12.5[!duplicated(rna12.5$Symbol),]
#rna13.5 <- rna13.5[!duplicated(rna13.5$Symbol),]

rna12.5$Alteration[rna12.5$Log2Rat < 0] <-  "down"
rna12.5$Alteration[rna12.5$Log2Rat > 0] <-  "up"
rna13.5$Alteration[rna13.5$Log2Rat < 0] <-  "down"
rna13.5$Alteration[rna13.5$Log2Rat > 0] <-  "up"

rna_e12.5_down <- subset(rna12.5, rna12.5$Log2Rat < 0, )
rna_e12.5_up <- subset(rna12.5, rna12.5$Log2Rat > 0, )
rna_e13.5_down <- subset(rna13.5, rna13.5$Log2Rat < 0, )
rna_e13.5_up <- subset(rna13.5, rna13.5$Log2Rat > 0, )

down12.5 <- unique(rna_e12.5_down$Symbol)
up12.5 <- unique(rna_e12.5_up$Symbol)
down13.5 <- unique(rna_e13.5_down$Symbol)
up13.5 <- unique(rna_e13.5_up$Symbol)

down12.5
up12.5 
down13.5 
up13.5

E12.5_down_up <- intersect(down12.5, up12.5)
E13.5_down_up <- intersect(down13.5, up13.5)

E12.5_down_up
E13.5_down_up


write.table(down12.5, "rna_down12.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(up12.5, "rna_up12.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(down13.5, "rna_down13.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(up13.5, "rna_up13.5.txt", row.names = F, col.names = F, quote = F, eol = " ")





#attach(gn)
#gn$Annotation[Distance.to.TSS >= -2000 & Distance.to.TSS <= 2000] <- "Promoter (¡À2k)"
#detach(gn)
#an02_tbl <- table(gn$Annotation)
#an02_lbls <- paste(names(an02_tbl), " ~ ", an02_tbl, sep="")
#pie(an02_tbl, labels=an02_lbls, main="Annotation of 17396 Peaks")


an05_tbl <- table(gn$Annotation)
an05_lbls <- paste(names(an05_tbl), "\n", an05_tbl, sep="")
library(plotrix)
pie3D(an05_tbl, 
      labels=an05_lbls, 
      cex = 2,
      cex.main = 3,
      explode=0.02,
      col=c("black", "azure", "beige", "palegreen"),
      main="\n \n \n \n Peak Annotation\n 17396"
      )



#library(Hmisc)
#describe(gn_05k$Gene.Name)



#gn02typ_tbl <- table(gn_02k$Gene.Type)[1:5]
#gn02typ_lbls <- paste(names(gn02typ_tbl), " ~ ", gn02typ_tbl, sep=" ")
#pie3D(gn02typ_tbl, labels=gn02typ_lbls,explode=0, main="Type of 534 Unique Genes")

gn05typ_tbl <- table(gn_05k$Gene.Type)[1:5]
gn05typ_tbl
gn05typ_lbls <- paste(names(gn05typ_tbl), " ", gn05typ_tbl, sep=" ")
pie(gn05typ_tbl, labels=gn05typ_lbls, 
    cex = 1.5,
    col=c("yellow", "green", "red", "white", "blue") ,
    cex.main = 2.5,
    main="\n \n \n \n Unique Gene Type\n 1213")

#gn10typ_tbl <- table(gn_10k$Gene.Type)[1:5]
#gn10typ_lbls <- paste(names(gn10typ_tbl), " ~ ", gn10typ_tbl, sep=" ")
#pie3D(gn10typ_tbl, labels=gn10typ_lbls,explode=0, main="Type of 2167 Unique Genes")






#down12.5 <- rna_e12.5_down[!duplicated(rna_e12.5_down$Symbol),]$Symbol
#up12.5 <- rna_e12.5_up[!duplicated(rna_e12.5_up$Symbol),]$Symbol
#down13.5 <- rna_e13.5_down[!duplicated(rna_e13.5_down$Symbol),]$Symbol
#up13.5 <- rna_e13.5_up[!duplicated(rna_e13.5_up$Symbol),]$Symbol


#install.packages("grid", "VennDiagram")
#library(grid)
#library(VennDiagram)

#source("http://www.bioconductor.org/biocLite.R")

#biocLite("limma")
#biocLite("statmod")
#library("limma")



library(grid)
library(VennDiagram)

venn.diagram(list(E12.5_Rseq=down12.5,E13.5_Rseq=down13.5, ChIPseq=gn_05k), 
             fill=c("darkorchid1","olivedrab1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Down_Both.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "With Downregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             force.unique = TRUE)


venn.diagram(list(E12.5_Rseq=up12.5,E13.5_Rseq=up13.5, ChIPseq=gn_05k), 
             fill=c("darkorchid4","olivedrab4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_Both.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "With Upregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             force.unique = TRUE)

venn.diagram(list(E12.5=down12.5,E13.5=down13.5), 
             fill=c("darkorchid1","olivedrab1"),
             filename = "Venn_RNAseq_Down.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "Downregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 180,
             cat.pos = 0,
             force.unique = TRUE)

venn.diagram(list(E12.5=up12.5,E13.5=up13.5), 
             fill=c("darkorchid4","olivedrab4"),
             filename = "Venn_RNAseq_Up.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "Upregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 180,
             cat.pos = 0,
             force.unique = TRUE)

venn.diagram(list(Rseq=down12.5,ChIPseq=gn_05k), 
             fill=c("darkorchid1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_down_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Downregulated E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=up12.5,ChIPseq=gn_05k), 
             fill=c("darkorchid4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Upregulated E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=down13.5,ChIPseq=gn_05k), 
             fill=c("olivedrab1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_down_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Downregulated E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list (Rseq=up13.5,ChIPseq=gn_05k), 
             fill = c("olivedrab4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Upregulated E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)




venn.diagram(list(Down=down12.5,Up=up12.5 ), 
             fill=c("darkorchid1","turquoise1"),
             filename = "Venn_RNAseq_Down_Up_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             #rotation.degree = 90,
             #cat.pos = 0, 90, 180,
             force.unique = TRUE)

venn.diagram(list(Down=down13.5,Up=up13.5 ), 
             fill=c("darkorchid4","turquoise4"),
             filename = "Venn_RNAseq_Down_Up_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             #rotation.degree = 90,
             #cat.pos = 0, 90, 180,
             force.unique = TRUE)

gn05k_down12.5 <- intersect (down12.5, gn_05k)
gn05k_up12.5 <- intersect (up12.5, gn_05k)
gn05k_down13.5 <- intersect (down13.5, gn_05k)
gn05k_up13.5 <- intersect (up13.5, gn_05k)
gn05k_down_intersection <- intersect (gn05k_down12.5, gn05k_down13.5)
gn05k_up_intersection <- intersect (gn05k_up12.5, gn05k_up13.5)
gn05k_down_union <- union (gn05k_down12.5, gn05k_down13.5)
gn05k_up_union <- union (gn05k_up12.5, gn05k_up13.5)

gn05k_down12.5 
gn05k_up12.5 
gn05k_down13.5 
gn05k_up13.5 
gn05k_down_intersection 
gn05k_up_intersection 
gn05k_down_union 
gn05k_up_union



write.table(gn05k_down12.5, "gn05k_down12.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_up12.5, "gn05k_up12.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_down13.5, "gn05k_down13.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_up13.5, "gn05k_up13.5.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_down_intersection, "gn05k_down_intersection.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_up_intersection , "gn05k_up_intersection.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_down_union, "gn05k_down_union.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn05k_up_union, "gn05k_up_union.txt", row.names = F, col.names = F, quote = F, eol = " ")


gn05k_go <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5/gn_05k_GOTERM_BP_FAT.txt")
gn05k_go$Term <- sub(".*~", "", gn05k_go$Term)
gn05k_go$Relevance <- -log(gn05k_go$PValue)
gn05K_go9 <- gn05k_go[1:9,]
#install.packages("ggplot2")
#install.packages("gcookbook")
library(ggplot2)

ggplot(gn05K_go9, aes(x=reorder(Term, PValue), y=Relevance)) +
  geom_point(size=5) + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(angle=-45, hjust=0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="gray60", linetype="dashed"))

