#Lab
#macs14 <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_gene.txt")

#Home
macs14 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_gene.txt")
rna_e12.5_down <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E12.5_M_E12.5_W/Tophat/gene_exp_E12.5_M_E12.5_W_selected_1.5fold_symbol_downregulated.diff")
rna_e12.5_up <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E12.5_M_E12.5_W/Tophat/gene_exp_E12.5_M_E12.5_W_selected_1.5fold_symbol_upregulated.diff")
rna_e13.5_down <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E13.5_M_E13.5_W/Tophat/gene_exp_E13.5_M_E13.5_W_selected_1.5fold_symbol_downregulated.diff")
rna_e13.5_up <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E13.5_M_E13.5_W/Tophat/gene_exp_E13.5_M_E13.5_W_selected_1.5fold_symbol_upregulated.diff")

rna12.5 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E12.5_M_E12.5_W/DEGseq/geneExp_E12.5_M_E12.5_W_filtered_DEGseq.txt")
rna13.5 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E13.5_M_E13.5_W/DEGseq/geneExp_E13.5_M_E13.5_W_filtered_DEGseq.txt")

names(macs14)[c(1, 22)] <- c("Peak.ID", "Smad4.Motif")
#substitute character after first Space with
#null A space (), then any character (.) any number of times (*) 
#until the end of the string ($). See ?regex to learn regular expressions. 
attach(macs14)
macs14$Annotation <- sub(" .N.*$", "", Annotation)
detach(macs14)
attach(macs14)
macs14$Annotation <- sub("intron", "Intron", Annotation)
detach(macs14)
attach(macs14)
macs14$Annotation <- sub("exon", "Exon", Annotation)
detach(macs14)
attach(macs14)
macs14$Annotation <- sub("3' UTR", "Intergenic", Annotation)
detach(macs14)
attach(macs14)
macs14$Annotation <- sub("5' UTR", "Intergenic", Annotation)
detach(macs14)
attach(macs14)
macs14$Annotation <- sub("non-coding", "Intergenic", Annotation)
detach(macs14)
attach(macs14)
macs14$Annotation <- sub("TTS", "Intergenic", Annotation)
detach(macs14)


#attach(macs14)
#macs14$Annotation[Distance.to.TSS >= -2000 & Distance.to.TSS <= 2000] <- "Promoter (±2k)"
#detach(macs14)
#an02_tbl <- table(macs14$Annotation)
#an02_lbls <- paste(names(an02_tbl), " ~ ", an02_tbl, sep="")
#pie(an02_tbl, labels=an02_lbls, main="Annotation of 17396 Peaks")

attach(macs14)
macs14$Annotation[Distance.to.TSS >= -5000 & Distance.to.TSS <= 5000] <- "Promoter (±5k)"
detach(macs14)
an05_tbl <- table(macs14$Annotation)
an05_lbls <- paste(names(an05_tbl), " ~ ", an05_tbl, sep="")
library(plotrix)
pie3D(an05_tbl, 
      labels=an05_lbls, 
      explode=0.02,
      col=c("black", "azure", "beige", "palegreen"),
      main="Annotation of 17396 Peaks")

#attach(macs14)
#macs14$Annotation[Distance.to.TSS >= -10000 & Distance.to.TSS <= 10000] <- "Promoter (±10k)"
#detach(macs14)
##an10_tbl <- table(macs14$Annotation)
#an10_lbls <- paste(names(an10_tbl), " ~ ", an10_tbl, sep="")
#pie(an10_tbl, labels=an10_lbls, main="Annotation of 17396 Peaks")
#attach(macs14)
#macs14_02K <- macs14[which(Distance.to.TSS >= -2000 & Distance.to.TSS <= 2000) , ]
#macs14_05k <- macs14[which(Distance.to.TSS >= -5000 & Distance.to.TSS <= 5000) , ]
#macs14_10k <- macs14[which(Distance.to.TSS >= -10000 & Distance.to.TSS <= 10000) , ]
#detach(macs14)

#macs14_02k <- subset(macs14, Distance.to.TSS >= -2000 & Distance.to.TSS <= 2000, 
#                     select=c(2, 3, 4, 6, 8, 10, 16, 19, 22))
macs14_05k <- subset(macs14, Distance.to.TSS >= -5000 & Distance.to.TSS <= 5000, 
                     select=c(2, 3, 4, 6, 8, 10, 16, 19, 22))
macs14_10k <- subset(macs14, Distance.to.TSS >= -10000 & Distance.to.TSS <= 10000, 
                     select=c(2, 3, 4, 6, 8, 10, 16, 19, 22))

write.table(macs14_05k, "u.txt")
#library(Hmisc)
#describe(macs14_05k$Gene.Name)
#macs14_02k <- macs14_02k[!duplicated(macs14_02k$Gene.Name),]
macs14_05k <- macs14_05k[!duplicated(macs14_05k$Gene.Name),]
macs14_10k <- macs14_10k[!duplicated(macs14_10k$Gene.Name),]

#demo(plotmath)


#gn02typ_tbl <- table(macs14_02k$Gene.Type)[1:5]
#gn02typ_lbls <- paste(names(gn02typ_tbl), " ~ ", gn02typ_tbl, sep=" ")
#pie3D(gn02typ_tbl, labels=gn02typ_lbls,explode=0, main="Type of 534 Unique Genes")

gn05typ_tbl <- table(macs14_05k$Gene.Type)[1:5]
gn05typ_lbls <- paste(names(gn05typ_tbl), " ~ ", gn05typ_tbl, sep=" ")
pie(gn05typ_tbl, labels=gn05typ_lbls,
    col=c("yellow", "green", "orange", "blue", "red") ,
    main="1213 Unique Gene Type")

#gn10typ_tbl <- table(macs14_10k$Gene.Type)[1:5]
#gn10typ_lbls <- paste(names(gn10typ_tbl), " ~ ", gn10typ_tbl, sep=" ")
#pie3D(gn10typ_tbl, labels=gn10typ_lbls,explode=0, main="Type of 2167 Unique Genes")

#rna_e12.5_down <- rna_e12.5_down[!duplicated(rna_e12.5_down$Symbol),]
#rna_e12.5_up <- rna_e12.5_up[!duplicated(rna_e12.5_up$Symbol),]
#rna_e13.5_down <- rna_e13.5_down[!duplicated(rna_e13.5_down$Symbol),]
#rna_e13.5_up <- rna_e13.5_up[!duplicated(rna_e13.5_up$Symbol),]

rna12.5 <- rna12.5[!duplicated(rna12.5$Symbol),]
rna13.5 <- rna13.5[!duplicated(rna13.5$Symbol),]

rna_e12.5_down <- subset(rna12.5, rna12.5$Log2Rat < 0, )
rna_e12.5_up <- subset(rna12.5, rna12.5$Log2Rat > 0, )
rna_e13.5_down <- subset(rna13.5, rna13.5$Log2Rat < 0, )
rna_e13.5_up <- subset(rna13.5, rna13.5$Log2Rat > 0, )


#gn02k <- macs14_02k$Gene.Name
gn05k <- macs14_05k$Gene.Name
gn10k <- macs14_10k$Gene.Name

down12.5 <- rna_e12.5_down$Symbol 
up12.5 <- rna_e12.5_up$Symbol 
down13.5 <- rna_e13.5_down$Symbol 
up13.5 <- rna_e13.5_up$Symbol 





#install.packages("grid", "VennDiagram")
#library(grid)
#library(VennDiagram)

#source("http://www.bioconductor.org/biocLite.R")

#biocLite("limma")
#biocLite("statmod")
#library("limma")

library(VennDiagram)

library(grid)

venn.diagram(list(E12.5_Rseq=down12.5,E13.5_Rseq=down13.5, ChIPseq=gn05k), 
             fill=c("darkorchid1","olivedrab1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Down_Both.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "With Downregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             force.unique = TRUE)

venn.diagram(list(E12.5_Rseq=up12.5,E13.5_Rseq=up13.5, ChIPseq=gn05k), 
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
             force.unique = TRUE)

venn.diagram(list(E12.5=up12.5,E13.5=up13.5), 
             fill=c("darkorchid4","olivedrab4"),
             filename = "Venn_RNAseq_Up.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "Upregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             force.unique = TRUE)

venn.diagram(list(Rseq=down12.5,ChIPseq=gn05k), 
             fill=c("darkorchid1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_down_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Downregulated E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=up12.5,ChIPseq=gn05k), 
             fill=c("darkorchid4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Upregulated E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=down13.5,ChIPseq=gn05k), 
             fill=c("olivedrab1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_down_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Downregulated E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=up13.5,ChIPseq=gn05k), 
             fill=c("olivedrab4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Upregulated E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

x <- intersect (down13.5, gn10k)
x
write.table(x, "x.txt")