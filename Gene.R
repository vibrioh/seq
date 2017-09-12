library(plotrix)
library(grid)
library(VennDiagram)
library(ggplot2)





gene <- read.delim("C:/Users/Jun Hu/Dropbox/Isl1/macs14peaks.txt")

g <- gene

g$Annotation <- sub(" .N.*$", "", g$Annotation)
g$Annotation <- sub(" .N.*$", "", g$Annotation)
g$Annotation <- sub("intron", "Intron", g$Annotation)
g$Annotation <- sub("exon", "Exon", g$Annotation)
g$Annotation <- sub("3' UTR", "Intergenic", g$Annotation)
g$Annotation <- sub("5' UTR", "Intergenic", g$Annotation)
g$Annotation <- sub("non-coding", "Intergenic", g$Annotation)
g$Annotation <- sub("TTS", "Intergenic", g$Annotation)
g$Annotation[g$Distance.to.TSS >= -10000 & g$Distance.to.TSS <= 10000] <- "Promoter (¡À10k)"

gn <- unique(g[g$Distance.to.TSS >= -10000 & g$Distance.to.TSS <= 10000, ]$Gene.Name)

write.table(gn, "gn.txt", row.names = F, col.names = F, quote = F, eol = " ")

g_an <- table(g$Annotation)

pie3D(g_an, 
      labels = paste(g_an, " ", names(g_an), sep=""),
      cex.main = 1.5,
      labelcex = 1.5,
      explode = 0.02,
      col = c("black", "azure", "beige", "palegreen"),
      main = "\n \n \n Peak Annotation \n 10812")

gp <- g[gn, ]

gn_tp <- table(factor(table(gp$Gene.Type)))


pie(gn_tp, 
    labels = paste(names(gn_tp), "\n", gn_tp, sep=" "), 
    cex = 2.5,
    #col=c("yellow", "green", "red", "white", "blue", "black") ,
    cex.main = 3.5,
    main="\n \n \n \n Unique Gene Type\n 2426")


pie(c(2078, 315, 33), labels=c("Protein-coding 2078", "miscRNA 315", "Pseudo 33"), main="\n \n  Unique Gene Type \n 2426")