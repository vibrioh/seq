###############################################################
#Install bioconductor && Package cummeRbund
###############################################################
source('http://www.bioconductor.org/biocLite.R')
biocLite('cummeRbund')

###############################################################
#Invoke cummeRbund 
###############################################################
library(cummeRbund)

###############################################################
#Read Data from Cuffdiff 
###############################################################
####Migrate Working Direction
getwd()
setwd("D:/mRNAseq/Working")
cuff_data  <- readCufflinks('diff_out')
cuff_data

###############################################################
#Quality Control for Data
###############################################################
disp <- dispersionPlot(genes(cuff_data))
disp
dispf <- dispersionPlot(cuff_data)
dispf
genes.scv <- fpkmSCVPlot(genes(cuff_data))
genes.scv
isoforms.scv <- fpkmSCVPlot(isoforms(cuff_data))
isoforms.scv
dens <- csDensity(genes(cuff_data))
dens
densRep <- csDensity(genes(cuff_data), replicates = T)
densRep
b <- csBoxplot(genes(cuff_data))
b
brep <- csBoxplot(genes(cuff_data), replicates = T)
brep
ssm <- csScatterMatrix(genes(cuff_data))
ssm
ss <- csScatter(genes(cuff_data), 'ZL', 'ZDF', smooth = T)
ss
m <- MAplot(genes(cuff_data), 'ZL', 'ZDF')
m
csv <- csVolcano(genes(cuff_data), 'ZL', 'ZDF', alpha = 0.05, xlimits = c(-30, 30), ylimits = c(0, 10),  showSignificant=T)
csv + scale_y_continuous(limits = c(0, 5))

###############################################################
#Get Data Infomation
###############################################################
runInfo(cuff_data)
replicates(cuff_data)

###############################################################
#Anatomize Data
###############################################################
gene_exp_origin <- read.delim("D:/mRNAseq/Working/diff_out/gene_exp.diff")
gene_exp <- gene_exp_origin[which(gene_exp_origin$gene != '-'), ]
gene_exp_p0.05 <- gene_exp[which(gene_exp$p_value < 0.05), ]
gene_exp_p0.05_up <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. > 0), ]
gene_exp_p0.05_down <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. < 0), ]
gene_exp_p0.05_up_f1 <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. > 1), ]
gene_exp_p0.05_up_f1.5 <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. > 1.5), ]
gene_exp_p0.05_up_f2 <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. > 2), ]
gene_exp_p0.05_down_f1 <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. < -1), ]
gene_exp_p0.05_down_f1.5 <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. < -1.5), ]
gene_exp_p0.05_down_f2 <- gene_exp_p0.05[which(gene_exp_p0.05$log2.fold_change. < -2), ]
write.csv(gene_exp, file = "gene_exp.cvs")
write.csv(gene_exp_p0.05, file = "gene_exp_p0.05.cvs")
write.csv(gene_exp_p0.05_up, file = "gene_exp_p0.05_up.cvs")
write.csv(gene_exp_p0.05_down, file = "gene_exp_p0.05_down.cvs")
write.csv(gene_exp_p0.05_up_f1, file = "gene_exp_p0.05_up_f1.cvs")
write.csv(gene_exp_p0.05_up_f1.5, file = "gene_exp_p0.05_up_f1.5.cvs")
write.csv(gene_exp_p0.05_up_f2, file = "gene_exp_p0.05_up_f2.cvs")
write.csv(gene_exp_p0.05_down_f1, file = "gene_exp_p0.05_down_f1.cvs")
write.csv(gene_exp_p0.05_down_f1.5, file = "gene_exp_p0.05_down_f1.5.cvs")
write.csv(gene_exp_p0.05_down_f2, file = "gene_exp_p0.05_down_f2.cvs")


###############################################################
#Gene Ontology
###############################################################
####Read Results
gene_exp_p0.05_down_Classification <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_Classification.txt", comment.char="#")
gene_exp_p0.05_down_GOTERM_BP_FAT <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_GOTERM_BP_FAT.txt")
gene_exp_p0.05_down_GOTERM_CC_FAT <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_GOTERM_CC_FAT.txt")
gene_exp_p0.05_down_GOTERM_MF_FAT <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_GOTERM_MF_FAT.txt")
gene_exp_p0.05_down_KEGG_PATHWAY <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_KEGG_PATHWAY.txt")
gene_exp_p0.05_down_SP_PIR_KEYWORDS <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_SP_PIR_KEYWORDS.txt")
gene_exp_p0.05_down_UP_SEQ_FEATURE <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_down_UP_SEQ_FEATURE.txt")
gene_exp_p0.05_up_Classification <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_Classification.txt", comment.char="#")
gene_exp_p0.05_up_GOTERM_BP_FAT <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_GOTERM_BP_FAT.txt")
gene_exp_p0.05_up_GOTERM_CC_FAT <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_GOTERM_CC_FAT.txt")
gene_exp_p0.05_up_GOTERM_MF_FAT <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_GOTERM_MF_FAT.txt")
gene_exp_p0.05_up_KEGG_PATHWAY <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_KEGG_PATHWAY.txt")
gene_exp_p0.05_up_SP_PIR_KEYWORDS <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_SP_PIR_KEYWORDS.txt")
gene_exp_p0.05_up_UP_SEQ_FEATURE <- read.delim("D:/mRNAseq/Working/gene_exp_p0.05_up_UP_SEQ_FEATURE.txt")
###Down
ggplot(gene_exp_p0.05_down_GOTERM_BP_FAT[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_down_GOTERM_CC_FAT[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_down_GOTERM_MF_FAT[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_down_KEGG_PATHWAY[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_down_SP_PIR_KEYWORDS[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_down_UP_SEQ_FEATURE[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
####Up
ggplot(gene_exp_p0.05_up_GOTERM_BP_FAT[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_up_GOTERM_CC_FAT[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_up_GOTERM_MF_FAT[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_up_KEGG_PATHWAY[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_up_SP_PIR_KEYWORDS[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()
ggplot(gene_exp_p0.05_up_UP_SEQ_FEATURE[1:10, ], aes(x = reorder(Term, -PValue), y = -log10(PValue))) + geom_bar(stat='identity') + coord_flip()

###############################################################
#Specify Genes of Interest
###############################################################
####Bar Plot
myGenesDown <- getGenes(cuff_data, c('Gja1', 'Opa1', 'Vdac1', 'Dnm1l', 'Kcnj11', 'Sod2', 'Hk1', 'Hk2'))
myGenesUp <- getGenes(cuff_data, c('Scn3a', 'Atp4a', 'Kcna5', 'Gja5'))
TSPO <- getGenes(cuff_data, c('Tspo', 'Actb'))
Bcl2Bax <- getGenes(cuff_data, c('Bcl2', 'Bax', 'Akt1'))
expressionBarplot(myGenesDown, showErrorbars = F, logMode = F)
expressionBarplot(myGenesUp, showErrorbars = F, logMode = F)
expressionBarplot(TSPO, showErrorbars = F, logMode = F)
expressionBarplot(Bcl2Bax, showErrorbars = F, logMode = F)
####Heatmap
HmyGenesDown <- csHeatmap(myGenesDown, replicates = T, logMode = F)
HmyGenesUp <- csHeatmap(myGenesUp, replicates = T, logMode = F)
HmyGenesDown
HmyGenesUp


