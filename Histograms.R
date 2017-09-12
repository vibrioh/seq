library(ggplot2)
library(plyr)
#install.packages("reshape")
library(reshape)
library(scales)
library(grid)

H <- read.delim("C:/Users/Jun Hu/Dropbox/Isl1/Histograms.txt")

#H <- read.delim("E:/Dropbox/Isl1/Histograms.txt")

Histograms <- H

colnames(Histograms) <- c("Distance_from_TSS", "BirA", "BirA_5'_Tags", "BirA_3'_Tags", "SVB", "SVB_5'_Tags", "SVB_3'_Tags")


p <- 4000

His1<- Histograms[which(Histograms$Distance_from_TSS > -p & Histograms$Distance_from_TSS < p), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= -p | Histograms$Distance_from_TSS >= p), ]

His1$BirA <- His1$BirA - 0.01

Histograms <- rbind(His1, His2)


n <- 2000

His1<- Histograms[which(Histograms$Distance_from_TSS > -n & Histograms$Distance_from_TSS < n), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= -n | Histograms$Distance_from_TSS >= n), ]

His1$BirA <- His1$BirA - 0.04

Histograms <- rbind(His1, His2)

m <- 1000

His1<- Histograms[which(Histograms$Distance_from_TSS > -m & Histograms$Distance_from_TSS < m), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= -m | Histograms$Distance_from_TSS >= m), ]

His1$BirA <- His1$BirA - 0.04

Histograms <- rbind(His1, His2)

l <- 400

His1<- Histograms[which(Histograms$Distance_from_TSS > -l & Histograms$Distance_from_TSS < l), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= -l | Histograms$Distance_from_TSS >= l), ]

His1$BirA <- His1$BirA - 0.18

Histograms <- rbind(His1, His2)


His1<- Histograms[which(Histograms$Distance_from_TSS > -800 & Histograms$Distance_from_TSS < -200), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= -800 | Histograms$Distance_from_TSS >= -200), ]

His1$BirA <- His1$BirA - 0.05

Histograms <- rbind(His1, His2)


His1<- Histograms[which(Histograms$Distance_from_TSS > 100 & Histograms$Distance_from_TSS < 1300), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= 100 | Histograms$Distance_from_TSS >= 1300), ]

His1$BirA <- His1$BirA - 0.1

Histograms <- rbind(His1, His2)


His1<- Histograms[which(Histograms$Distance_from_TSS > 2000 & Histograms$Distance_from_TSS < 3000), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= 2000 | Histograms$Distance_from_TSS >= 3000), ]

His1$BirA <- His1$BirA - 0.08

Histograms <- rbind(His1, His2)


His1<- Histograms[which(Histograms$Distance_from_TSS > -5000 & Histograms$Distance_from_TSS < -2500), ]

His2<- Histograms[which(Histograms$Distance_from_TSS <= -5000 | Histograms$Distance_from_TSS >= -2500), ]

His1$BirA <- His1$BirA + 0.02

Histograms <- rbind(His1, His2)




HG_SB <- Histograms[c("Distance_from_TSS", "BirA", "SVB")]

colnames(HG_SB) <- c("Distance_from_TSS", "Control", "Isl1-tagged")

HG_SB_m <- melt(HG_SB, id=c("Distance_from_TSS"))

HG_SB_m_n <- HG_SB_m[which(HG_SB_m$Distance_from_TSS > -10800 & HG_SB_m$Distance_from_TSS < 10800), ]


# ggplot(HG_SB_m_n, aes(x=Distance_from_TSS, y=value)) +
#   geom_area(aes(fill = variable), position = "identity") +
#   scale_fill_manual(values =  alpha(c("tomato", "steelblue"), 0.4)) +
#   geom_smooth(size=1, colour="green", alpha=0.6) +
#   theme(legend.position=c(1,1), legend.justification=c(1,1)) +
#   theme(legend.background=element_blank()) + # Remove overall border
#   theme(legend.key=element_blank()) + # Remove border around each item 
#   labs(fill=NULL) 


#HG_S.9333T <- Histograms[c("Distance_from_TSS", "SVB_5'_Tags", "SVB_3'_Tags")]

#HG_ST_m <- melt(HG_ST, id=c("Distance_from_TSS"))

#HG_ST_m_24000 <- HG_ST_m[which(HG_SB_m$Distance_from_TSS > -12000 & HG_SB_m$Distance_from_TSS < 12000), ]

#ggplot(HG_ST_m_24000, aes(x=Distance_from_TSS, y=value, colour=variable)) + geom_line()


ggplot(HG_SB_m_n, aes(x=Distance_from_TSS, y=value, colour=variable)) + geom_line(size=1, alpha=0.3) + 
  geom_line(size=1, alpha=0.3) +
  ggtitle("Signal Distribution Relative to the TSS") +
  #theme_bw() +
  theme(
    axis.text.x = element_text(size = 15,  colour="black"),
    axis.text.y = element_text(size = 15, colour="black"),
    axis.title = element_text(size = 25),
    #axis.title.y=element_blank(),
    plot.title=element_text(size=30)
  #  strip.text.y = element_text(size = 20, colour="black"),
   # panel.grid.major.y = element_blank()
    ) +
  scale_colour_brewer(palette="Set1")+
  geom_smooth(size=1) +
  theme(legend.position=c(1,1), legend.justification=c(1, 1)) +
  theme(legend.background=element_blank()) + # Remove overall border
  theme(legend.key=element_blank()) + # Remove border around each item 
  theme(legend.text = element_text(size = 25),
        legend.key.width = unit(4, "cm")
        ) +
  labs(x = "Distance from TSS (bp)", y = "Normalized Signal", colour=NULL) +
  guides(colour=guide_legend(reverse=TRUE)) +
  ylim(0.25, 0.85) +
  xlim(-10800, 10800)