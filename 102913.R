library(plotrix)
library(grid)
library(ggplot2)
library(tcltk)
library(VennDiagram)
library(bvenn)
library(eVenn)
library(extrafont)



#install.packages("extrafont")

font_import()
loadfonts(device = "win")


go_gn$Term <- sub(".*~", "", go_gn$Term)
go_gn$Regulation <- "chip"

go_gn_termorder <- go_gn$Term[order(go_gn$Regulation, -go_gn$PValue)]
go_gn$Term <- factor(go_gn$Term, levels=go_gn_termorder)
ggplot(go_gn, aes(x=-log10(PValue), y=Term)) +
  ggtitle("Gene Ontology") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 10^-03, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go_gn$Term, lable=go_gn$Count, data=go_gn) +
  geom_segment(aes(yend=Term), xend=0, size=6, colour = "cyan4") +
  geom_point(size=18, colour = "cyan4") +
  #scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  #theme_bw() +
  theme(
    axis.text.x = element_text(size = 38, face="bold", colour="black"),
    axis.text.y = element_text(size = 38, face="plain", colour="black"),
    axis.title = element_text(size = 38),
    axis.title.y=element_blank(),
    plot.title=element_text(size=46, face="bold"),
    strip.text.y = element_text(size = 38, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=10, fontface="bold.italic", colour="cornsilk") 

Go_F <- Go_F[1:5, ]
Go_F$Term <- sub(".*~", "", Go_F$Term)
Go_F$Regulation <- "chip"

Go_F_termorder <- Go_F$Term[order(Go_F$Regulation, -Go_F$PValue)]
Go_F$Term <- factor(Go_F$Term, levels=Go_F_termorder)
ggplot(Go_F, aes(x=-log10(PValue), y=Term)) +
  ggtitle("Functional Categories") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 10^-03, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=Go_F$Term, lable=Go_F$Count, data=Go_F) +
  geom_segment(aes(yend=Term), xend=0, size=6, colour = "palegreen3") +
  geom_point(size=18, colour = "palegreen3") +
  #scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  #theme_bw() +
  theme(
    axis.text.x = element_text(size = 38, face="bold", colour="black"),
    axis.text.y = element_text(size = 38, face="plain", colour="black"),
    axis.title = element_text(size = 38),
    axis.title.y=element_blank(),
    plot.title=element_text(size=46, face="bold"),
    strip.text.y = element_text(size = 38, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=10, fontface="bold.italic", colour="cornsilk") 



###############################################################################################################################

go_r12.5_termorder <- go_r12.5$Term[order(go_r12.5$Regulation, -go_r12.5$PValue)]
go_r12.5$Term <- factor(go_r12.5$Term, levels=go_r12.5_termorder)
ggplot(go_r12.5, aes(x=-log10(PValue), y=Term)) +
  ggtitle("E12.5") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 0.001, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go_r12.5$Term, lable=go_r12.5$Count, data=go_r12.5) +
  geom_segment(aes(yend=Term, colour = Regulation), xend=0, size=6) +
  geom_point(size=18, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  #theme_bw() +
  theme(
    axis.text.x = element_text(size = 40, face="bold", colour="black"),
    axis.text.y = element_text(size = 46, face="plain", colour="black"),
    axis.title = element_text(size = 40),
    axis.title.y=element_blank(),
    plot.title=element_text(size=46, face="bold"),
    strip.text.y = element_text(size = 46, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=9, fontface="bold.italic", colour="white") +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y") 

############################################################################################################################

go_r13.5_termorder <- go_r13.5$Term[order(go_r13.5$Regulation, -go_r13.5$PValue)]
go_r13.5$Term <- factor(go_r13.5$Term, levels=go_r13.5_termorder)
ggplot(go_r13.5, aes(x=-log10(PValue), y=Term)) +
  ggtitle("E13.5") +
  scale_x_continuous(breaks=c(5*(0:10)),
                     labels=expression(1, 10^-05, 10^-10, 10^-15, 10^-20, 10^-25, 10^-30, 10^-35, 10^-40, 10^-45, 10^-50),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go_r12.5$Term, lable=go_r12.5$Count, data=go_r12.5) +
  geom_segment(aes(yend=Term, colour = Regulation), xend=0, size=6) +
  geom_point(size=18, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  #theme_bw() +
  theme(
    axis.text.x = element_text(size = 40, face="bold", colour="black"),
    axis.text.y = element_text(size = 46, face="plain", colour="black"),
    axis.title = element_text(size = 40),
    axis.title.y=element_blank(),
    plot.title=element_text(size=46, face="bold"),
    strip.text.y = element_text(size = 46, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=9, fontface="bold.italic", colour="white") +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y") 

##########################################################################################################################

go_r14.5_termorder <- go_r14.5$Term[order(go_r14.5$Regulation, -go_r14.5$PValue)]
go_r14.5$Term <- factor(go_r14.5$Term, levels=go_r14.5_termorder)
ggplot(go_r14.5, aes(x=-log10(PValue), y=Term)) +
  ggtitle("E14.5") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 0.001, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go_r14.5$Term, lable=go_r14.5$Count, data=go_r14.5) +
  geom_segment(aes(yend=Term), xend=0, colour="grey50") +
  geom_point(size=11, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face="bold"),
    axis.text.y = element_text(size = 19, face="bold"),
    axis.title = element_text(size = 19),
    axis.title.y=element_blank(),
    plot.title=element_text(size=25),
    strip.text.y = element_text(size = 18, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=5, fontface="bold.italic", colour="white") +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y")


##################################################################################################################################

go_in_r12.5_r13.5_termorder <- go_in_r12.5_r13.5$Term[order(go_in_r12.5_r13.5$Regulation, -go_in_r12.5_r13.5$PValue)]
go_in_r12.5_r13.5$Term <- factor(go_in_r12.5_r13.5$Term, levels=go_in_r12.5_r13.5_termorder)
ggplot(go_in_r12.5_r13.5, aes(x=-log10(PValue), y=Term)) +
  ggtitle("E12.5 กษ E13.5") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 0.001, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go_in_r12.5_r13.5$Term, lable=go_in_r12.5_r13.5$Count, data=go_in_r12.5_r13.5) +
  geom_segment(aes(yend=Term), xend=0, colour="grey50") +
  geom_point(size=11, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face="bold"),
    axis.text.y = element_text(size = 19, face="bold"),
    axis.title = element_text(size = 19),
    axis.title.y=element_blank(),
    plot.title=element_text(size=25),
    strip.text.y = element_text(size = 18, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=5, fontface="bold.italic", colour="white") +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y")

###################################################################################################################33

go_in_r13.5_r14.5_termorder <- go_in_r13.5_r14.5$Term[order(go_in_r13.5_r14.5$Regulation, -go_in_r13.5_r14.5$PValue)]
go_in_r13.5_r14.5$Term <- factor(go_in_r13.5_r14.5$Term, levels=go_in_r13.5_r14.5_termorder)
ggplot(go_in_r13.5_r14.5, aes(x=-log10(PValue), y=Term)) +
  ggtitle("E13.5 กษ E14.5") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 0.001, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go_in_r13.5_r14.5$Term, lable=go_in_r13.5_r14.5$Count, data=go_in_r13.5_r14.5) +
  geom_segment(aes(yend=Term), xend=0, colour="grey50") +
  geom_point(size=11, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Up-regulated","Down-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face="bold"),
    axis.text.y = element_text(size = 19, face="bold"),
    axis.title = element_text(size = 19),
    axis.title.y=element_blank(),
    plot.title=element_text(size=25),
    strip.text.y = element_text(size = 18, colour="black"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=5, fontface="bold.italic", colour="white") +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y")

















go <- rbind(go_r12.5, go_r13.5, go_r14.5, go_in_r12.5_r13.5, go_in_r13.5_r14.5)

go$order <- factor(c(80:1))

#go_r12.5_termorder <- go_r12.5$Term[order(go_r12.5$Regulation, -go$PValue)]
#go_r12.5$Term <- factor(go_r12.5$Term, levels=go_r12.5_termorder)

ggplot(go, aes(x=-log10(PValue), y=order)) +
  ggtitle("Gene Ontology Term") +
  scale_x_continuous(breaks=c(0:50),
                     labels=expression(1, 0.1, 0.01, 0.001, 10^-04, 10^-05, 10^-06, 10^-07, 10^-08, 10^-09, 
                                       10^-10, 10^-11, 10^-12, 10^-13, 10^-14, 10^-15, 10^-16, 10^-17, 10^-18, 
                                       10^-19, 10^-20, 10^-21, 10^-22, 10^-23, 10^-24, 10^-25, 10^-26, 10^-27, 
                                       10^-28, 10^-29, 10^-30, 10^-31, 10^-32, 10^-33, 10^-34, 10^-35, 10^-36, 
                                       10^-37, 10^-38, 10^-39, 10^-40, 10^-41, 10^-42, 10^-43, 10^-44, 10^-45, 
                                       10^-46, 10^-47, 10^-48, 10^-49, 10^-50),
                     name="P-Value") +
  scale_y_discrete(breaks=go$order,
                   labels=go$Term)+
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go$Term, lable=go$Count, data=go) +
  geom_segment(aes(yend=order), xend=0, colour="grey50") +
  geom_point(size=11, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Down-regulated","Up-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face="bold"),
    axis.text.y = element_text(size = 18, face="bold"),
    axis.title = element_text(size = 23),
    axis.title.y=element_blank(),
    plot.title=element_text(size=25),
    strip.text.y = element_text(size = 18, colour="white", face="bold"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=5, fontface="bold.italic", colour="white") +
  facet_grid(Stage ~ ., scales="free_y", space="free_y")

go <- rbind(go_r12.5, go_r13.5, go_in_r12.5_r13.5)

go$order <- factor(c(48:1))

#go_r12.5_termorder <- go_r12.5$Term[order(go_r12.5$Regulation, -go$PValue)]
#go_r12.5$Term <- factor(go_r12.5$Term, levels=go_r12.5_termorder)

ggplot(go, aes(x=-log10(PValue), y=order)) +
  ggtitle("Gene Ontology Term") +
  scale_x_continuous(breaks=c(0:50),
                     labels=expression(1, 0.1, 0.01, 0.001, 10^-04, 10^-05, 10^-06, 10^-07, 10^-08, 10^-09, 
                                       10^-10, 10^-11, 10^-12, 10^-13, 10^-14, 10^-15, 10^-16, 10^-17, 10^-18, 
                                       10^-19, 10^-20, 10^-21, 10^-22, 10^-23, 10^-24, 10^-25, 10^-26, 10^-27, 
                                       10^-28, 10^-29, 10^-30, 10^-31, 10^-32, 10^-33, 10^-34, 10^-35, 10^-36, 
                                       10^-37, 10^-38, 10^-39, 10^-40, 10^-41, 10^-42, 10^-43, 10^-44, 10^-45, 
                                       10^-46, 10^-47, 10^-48, 10^-49, 10^-50),
                     name="P-Value") +
  scale_y_discrete(breaks=go$order,
                   labels=go$Term)+
  #annotation_logticks() + 
  #xlim(0, 20) +
  #scale_x_reverse() +
  #geom_text(x=50, y=go$Term, lable=go$Count, data=go) +
  geom_segment(aes(yend=order), xend=0, colour="grey50") +
  geom_point(size=11, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Down-regulated","Up-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face="bold"),
    axis.text.y = element_text(size = 18, face="bold"),
    axis.title = element_text(size = 23),
    axis.title.y=element_blank(),
    plot.title=element_text(size=25),
    strip.text.y = element_text(size = 18, colour="white", face="bold"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=5, fontface="bold.italic", colour="white") +
  facet_grid(Stage ~ ., scales="free_y", space="free_y")





################################################################################################################


go_r12.5 <- rbind(go_r12.5u[1:5, ], go_r12.5d[1:5, ])
go_r13.5 <- rbind(go_r13.5u[1:5, ], go_r13.5d[1:5, ])
go_r14.5 <- rbind(go_r14.5u[1:5, ], go_r14.5d[1:5, ])
go_in_r12.5_r13.5 <- rbind(go_in_r12.5u_r13.5u[1:5, ], go_in_r12.5d_r13.5d[1:5, ])
go_in_r13.5_r14.5 <- rbind(go_in_r13.5u_r14.5u[1:5, ], go_in_r13.5d_r14.5d[1:5, ])

go_r12.5$Term <- sub(".*~", "", go_r12.5$Term)
go_r13.5$Term <- sub(".*~", "", go_r13.5$Term)
go_r14.5$Term <- sub(".*~", "", go_r14.5$Term)
go_in_r12.5_r13.5$Term <- sub(".*~", "", go_in_r12.5_r13.5$Term)
go_in_r13.5_r14.5$Term <- sub(".*~", "", go_in_r13.5_r14.5$Term)

go <- rbind(go_r12.5, go_r13.5, go_in_r12.5_r13.5)

go$order <- factor(c(30:1))

ggplot(go, aes(x=-log10(PValue), y=order)) +
  ggtitle("Gene Ontology Term") +
  scale_x_continuous(breaks=c(3*(0:10)),
                     labels=expression(1, 0.001, 10^-06, 10^-09, 10^-12, 10^-15, 10^-18, 10^-21, 10^-24, 10^-27, 10^-30),
                     name="P-Value") +
  scale_y_discrete(breaks=go$order,
                   labels=go$Term)+
  geom_segment(aes(yend=order), xend=0, colour="grey50") +
  geom_point(size=10, aes(colour=Regulation, shape=Stage)) +
  #scale_colour_brewer(palette="Set1", limits=c("Down-regulated","Up-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, face="italic"),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18, face="italic"),
    axis.title.y=element_blank(),
    plot.title=element_text(size=23),
    strip.text.y = element_text(size = 18, colour="white", face="bold"),
    panel.grid.major.y = element_blank()) +
  geom_text(aes(label=Count), size=5, fontface="bold.italic", colour="white") +
  facet_grid(Stage ~ ., scales="free_y", space="free_y")




##############################################################################################################


























go_r13.5_termorder <- go_r13.5$Term[order(go_r13.5$Regulation, -go_r13.5$PValue)]
go_r13.5$Term <- factor(go_r13.5$Term, levels=go_r13.5_termorder)
ggplot(go_r13.5, aes(x=-log(PValue), y=Term)) +
  geom_segment(aes(yend=Term), xend=0, colour="grey50") +
  geom_point(size=4, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Down-regulated","Up-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 17),
    panel.grid.major.y = element_blank()) +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y")

go_in_r_termorder <- go_in_r$Term[order(go_in_r$Regulation, -go_in_r$PValue)]
go_in_r$Term <- factor(go_in_r$Term, levels=go_in_r_termorder)
ggplot(go_in_r, aes(x=-log(PValue), y=Term)) +
  geom_segment(aes(yend=Term), xend=0, colour="grey50") +
  geom_point(size=4, aes(colour=Regulation)) +
  scale_colour_brewer(palette="Set1", limits=c("Down-regulated","Up-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 17),
    panel.grid.major.y = element_blank()) +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y")

go_r_all <- rbind (go_r12.5, go_r13.5, go_in_r)
go_r_all_termorder <- go_r12.5$Term[order(go_r_all$Stage, go_r_all$Regulation, -go_r_all$PValue)]
go_r_all$Term <- factor(go_r_all$Term, levels=go_r_all_termorder)

go_r_all_termorder <- go_r_all$Term[order(go_r_all$Regulation, -go_r_all$PValue)]
go_r_all$Term <- factor(go_r_all$Term, levels=go_r_all_termorder)
ggplot(go_r_all, aes(x=-log(PValue), y=Term)) +
  geom_segment(aes(yend=Term), xend=0, colour="grey50") +
  geom_point(size=4, aes(colour=Regulation, shape=Stage)) +
  scale_colour_brewer(palette="Set1", limits=c("Down-regulated","Up-regulated"), guide=FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 17),
    panel.grid.major.y = element_blank()) +
 # facet_grid(Stage ~ ., scales="free_y", space="free_y") +
  facet_grid(Regulation ~ ., scales="free_y", space="free_y")

