library(plotrix)
library(grid)
library(VennDiagram)
library(ggplot2)

#g <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_gene.txt") #lab
g <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_gene.txt") #home
     g$Annotation <- sub(" .N.*$", "", g$Annotation)
     g$Annotation <- sub("intron", "Intron", g$Annotation)
     g$Annotation <- sub("exon", "Exon", g$Annotation)
     g$Annotation <- sub("3' UTR", "Intergenic", g$Annotation)
     g$Annotation <- sub("5' UTR", "Intergenic", g$Annotation)
     g$Annotation <- sub("non-coding", "Intergenic", g$Annotation)
     g$Annotation <- sub("TTS", "Intergenic", g$Annotation)
     g$Annotation[g$Distance.to.TSS >= -5000 & g$Distance.to.TSS <= 5000] <- "Promoter (¡À5k)"

gn <- unique(g[g$Distance.to.TSS >= -5000 & g$Distance.to.TSS <= 5000, ]$Gene.Name)
write.table(gn, "gn.txt", row.names = F, col.names = F, quote = F, eol = " ")

g_an <- table(g$Annotation)
pie3D(g_an, 
      labels = paste(g_an, " ", names(g_an), sep=""),
      cex.main = 4,
      labelcex = 3,
      explode = 0.02,
      col = c("black", "azure", "beige", "palegreen"),
      main = "\n \n \n \n Peak Annotation\n 17396")

gn_tp <- table(g[gn, ]$Gene.Type)
pie(gn_tp, 
    labels = paste(names(gn_tp), "\n", gn_tp, sep=" "), 
    cex = 2.5,
    col=c("yellow", "green", "red", "white", "blue", "black") ,
    cex.main = 3.5,
    main="\n \n \n \n Unique Gene Type\n 1213")

#r12.5 <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/Jianyu/E12.5_M_E12.5_W/DEGseq/geneExp_E12.5_M_E12.5_W_filtered_DEGseq.txt") #lab
#r13.5 <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/Jianyu/E13.5_M_E13.5_W/DEGseq/geneExp_E13.5_M_E13.5_W_filtered_DEGseq.txt") #lab
r12.5 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E12.5_M_E12.5_W/DEGseq/geneExp_E12.5_M_E12.5_W_filtered_DEGseq.txt") #home
r13.5 <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/Jianyu/E13.5_M_E13.5_W/DEGseq/geneExp_E13.5_M_E13.5_W_filtered_DEGseq.txt") #home
         
         r12.5$alt[r12.5$Log2Rat < 0] <-  "down"
         r12.5$alt[r12.5$Log2Rat > 0] <-  "up"
         r13.5$alt[r13.5$Log2Rat < 0] <-  "down"
         r13.5$alt[r13.5$Log2Rat > 0] <-  "up"

r12.5d <- unique(r12.5[which(r12.5$alt == "down"), ]$Symbol)
r12.5u <- unique(r12.5[which(r12.5$alt == "up"), ]$Symbol)
r13.5d <- unique(r13.5[which(r13.5$alt == "down"), ]$Symbol)
r13.5u <- unique(r13.5[which(r13.5$alt == "up"), ]$Symbol)

gn_r12.5d <- intersect (r12.5d, gn)
gn_r12.5u <- intersect (r12.5u, gn)
gn_r13.5d <- intersect (r13.5d, gn)
gn_r13.5u <- intersect (r13.5u, gn)
gn_down_intersection <- intersect (gn_r12.5d, gn_r13.5d)
gn_up_intersection <- intersect (gn_r12.5u, gn_r13.5u)
gn_down_union <- union (gn_r12.5d, gn_r13.5d)
gn_up_union <- union (gn_r12.5u, gn_r13.5u)

gn_r12.5d 
gn_r12.5u 
gn_r13.5d 
gn_r13.5u 
gn_down_intersection 
gn_up_intersection 
gn_down_union 
gn_up_union



write.table(gn_r12.5d, "gn_r12.5d.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_r12.5u, "gn_r12.5u.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_r13.5d, "gn_r13.5d.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_r13.5u, "gn_r13.5u.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_down_intersection, "gn_down_intersection.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_up_intersection , "gn_up_intersection.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_down_union, "gn_down_union.txt", row.names = F, col.names = F, quote = F, eol = " ")
write.table(gn_up_union, "gn_up_union.txt", row.names = F, col.names = F, quote = F, eol = " ")


venn.diagram(list(E12.5_Rseq=r12.5d,E13.5_Rseq=r13.5d, ChIPseq=gn), 
             fill=c("darkorchid1","olivedrab1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Down_Both.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "With Downregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             force.unique = TRUE)


venn.diagram(list(E12.5_Rseq=r12.5u,E13.5_Rseq=r13.5u, ChIPseq=gn), 
             fill=c("darkorchid4","olivedrab4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_Both.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "With Upregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             force.unique = TRUE)

venn.diagram(list(E12.5=r12.5d,E13.5=r13.5d), 
             fill=c("darkorchid1","olivedrab1"),
             filename = "Venn_RNAseq_Down.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "Downregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 180,
             cat.pos = 0,
             force.unique = TRUE)

venn.diagram(list(E12.5=r12.5u,E13.5=r13.5u), 
             fill=c("darkorchid4","olivedrab4"),
             filename = "Venn_RNAseq_Up.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "Upregulated",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 180,
             cat.pos = 0,
             force.unique = TRUE)

venn.diagram(list(Rseq=r12.5d,ChIPseq=gn), 
             fill=c("darkorchid1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_down_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Downregulated E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=r12.5u,ChIPseq=gn), 
             fill=c("darkorchid4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Upregulated E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Rseq=r13.5d,ChIPseq=gn), 
             fill=c("olivedrab1","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_down_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Downregulated E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list (Rseq=r13.5u,ChIPseq=gn), 
             fill = c("olivedrab4","turquoise1"),
             filename = "Venn_ChIPseq_RNAseq_Up_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "ChIPseq | RNAseq", sub = "Upregulated E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             rotation.degree = 90,
             force.unique = TRUE)

venn.diagram(list(Down=r12.5d,Up=r12.5u ), 
             fill=c("darkorchid1","turquoise1"),
             filename = "Venn_RNAseq_Down_Up_12.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "E12.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             #rotation.degree = 90,
             #cat.pos = 0, 90, 180,
             force.unique = TRUE)

venn.diagram(list(Down=r13.5d,Up=r13.5u ), 
             fill=c("darkorchid4","turquoise4"),
             filename = "Venn_RNAseq_Down_Up_13.5.tiff", 
             height = 5800, width = 5800, resolution = 960, units = "px",  main = "RNAseq", sub = "E13.5",
             main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", main.col = "black", main.cex = 2, main.just = c(0.5, 1),
             sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", sub.col = "gray23", sub.cex = 1.5, sub.just = c(0.5, 1),
             #rotation.degree = 90,
             #cat.pos = 0, 90, 180,
             force.unique = TRUE)



#go_gn <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn.txt") #lab
go_gn <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn.txt") #home

go_gn$Term <- sub(".*~", "", go_gn$Term)
go_gn$Relevance <- -log(go_gn$PValue)

ggplot(go_gn, aes(x=reorder(Term, PValue), y=Relevance)) +
  geom_point(size=5) + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, size = 30),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour="gray60", linetype="dashed", size = 0.8))


#go_gn_r12.5d <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r12.5d.txt") #lab
go_gn_r12.5d <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r12.5d.txt") #home

go_gn_r12.5d$Term <- sub(".*~", "", go_gn_r12.5d$Term)
go_gn_r12.5d$Relevance <- -log(go_gn_r12.5d$PValue)

ggplot(go_gn_r12.5d, aes(y=reorder(Term, Relevance), x=Relevance)) +
  geom_point(size = 8, col = "green") + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 40),
        axis.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour="gray60", linetype="dashed", size = 0.8))


#go_gn_r12.5u <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r12.5u.txt") #lab
go_gn_r12.5u <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r12.5u.txt") #home

go_gn_r12.5u$Term <- sub(".*~", "", go_gn_r12.5u$Term)
go_gn_r12.5u$Relevance <- -log(go_gn_r12.5u$PValue)

ggplot(go_gn_r12.5u, aes(y=reorder(Term, Relevance), x=Relevance)) +
  geom_point(size = 8, col = "red") + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 40),
        axis.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour="gray60", linetype="dashed", size = 0.8))

#go_gn_r13.5d <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r13.5d.txt") #lab
go_gn_r13.5d <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r13.5d.txt") #home

go_gn_r13.5d$Term <- sub(".*~", "", go_gn_r13.5d$Term)
go_gn_r13.5d$Relevance <- -log(go_gn_r13.5d$PValue)

ggplot(go_gn_r13.5d, aes(y=reorder(Term, Relevance), x=Relevance)) +
  geom_point(size = 8, col = "green") + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 40),
        axis.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour="gray60", linetype="dashed", size = 0.8))


#go_gn_r13.5u <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r13.5u.txt") #lab
go_gn_r13.5u <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_r13.5u.txt") #home

go_gn_r13.5u$Term <- sub(".*~", "", go_gn_r13.5u$Term)
go_gn_r13.5u$Relevance <- -log(go_gn_r13.5u$PValue)

ggplot(go_gn_r13.5u, aes(y=reorder(Term, Relevance), x=Relevance)) +
  geom_point(size = 8, col = "red") + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 40),
        axis.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour="gray60", linetype="dashed", size = 0.8))

#go_gn_down_intersection <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_down_intersection.txt") #lab
go_gn_down_intersection <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_down_intersection.txt") #home

go_gn_down_intersection$Term <- sub(".*~", "", go_gn_down_intersection$Term)
go_gn_down_intersection$Relevance <- -log(go_gn_down_intersection$PValue)

ggplot(go_gn_down_intersection, aes(y=reorder(Term, Relevance), x=Relevance)) +
  geom_point(size = 8, col = "green") + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 40),
        axis.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour="gray60", linetype="dashed", size = 0.8))

#go_gn_down_union <- read.delim("C:/Users/Jun Hu/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_down_union.txt") #lab
go_gn_down_union <- read.delim("C:/Users/Vibrioh/Dropbox/Smad4/MACS14_p5_Degseq__all/go_gn_down_union.txt") #home

go_gn_down_union$Term <- sub(".*~", "", go_gn_down_union$Term)
go_gn_down_union$Relevance <- -log(go_gn_down_union$PValue)

ggplot(go_gn_down_union, aes(y=reorder(Term, Relevance), x=Relevance)) +
  geom_point(size = 8, col = "green") + # Use a larger dot
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 40),
        axis.title = element_text(size = 25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour="gray60", linetype="dashed", size = 0.8))








