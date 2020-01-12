# sufig1
options(stringsAsFactors = FALSE)
options(scipen = 100)


library(basicTrendline)
library(broom)
library(BuenColors)
library(CAGEr)
library(data.table)
library(DESeq2)
library(dplyr)
library(genefilter)
library(ggalt)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gplots)
library(gridExtra)
library(heatmap.plus)
library(nlstools)
library(pheatmap)
library(purrr)
library(RColorBrewer)
library(reshape)
library(rlist)
library(Rmisc)
library(rsample)
library(rstatix)
library(scales)
library(splines)
library(statmod)
library(stringr)



setwd("~/zjw/scCAT-seq-master/Analysing/")
load("sufig1.RData")



### sufig1 c




#########call peak
for(i in grep("myCAGEset",objects(),value = T)) {
  a <- paste('normalizeTagCount(',i,', method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
  print(a)
  eval(parse(text=a))
}














tc_C1_CAGE_ERCC <- tagClusters(myCAGEsetC1_CAGE_ERCC)[[1]]
tc_STRT_ERCC <- tagClusters(myCAGEsetSTRT_ERCC)[[1]]
tc_NAR_ERCC <- tagClusters(myCAGEsetNAR_ERCC)[[1]]
tc_CAT_ERCC <- tagClusters(myCAGEsetERCC_5cap_filteryes)[[1]]


tc_C1_CAGE_ERCC[,5] <- as.character(tc_C1_CAGE_ERCC[,5])
tc_STRT_ERCC[,5] <- as.character(tc_STRT_ERCC[,5])
tc_NAR_ERCC[,5] <- as.character(tc_NAR_ERCC[,5])
tc_CAT_ERCC[,5] <- as.character(tc_CAT_ERCC[,5])

tc_C1_CAGE_ERCC <- tc_C1_CAGE_ERCC[tc_C1_CAGE_ERCC[,5]=="+",]
tc_STRT_ERCC <- tc_STRT_ERCC[tc_STRT_ERCC[,5]=="+",]
tc_NAR_ERCC <- tc_NAR_ERCC[tc_NAR_ERCC[,5]=="+",]
tc_CAT_ERCC <- tc_CAT_ERCC[tc_CAT_ERCC[,5]=="+",]











#3

temp <- data.frame()
for(i in 1:20) {
  temp <- rbind(temp,
                data.frame(position=i,peaknumber=length(tc_CAT_ERCC[tc_CAT_ERCC[,8]==i,6])))
}
temp <- rbind(temp,data.frame(position=21,peaknumber=length(tc_CAT_ERCC[tc_CAT_ERCC[,8]>=21,6]))) 



ggplot(temp,aes(x=position,y=peaknumber)) +
  geom_bar(position=position_dodge(0.8),width=0.8,stat="identity") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 12),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(1,21),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",">20"))+
  #scale_y_continuous(limits = c(0,1000000),breaks = seq(0,1000000,250000),
  #                   labels = c("0","25","50","75","100"))+
  labs(x=element_blank())








temp <- data.frame()
for(i in 1:20) {
  temp <- rbind(temp,
                data.frame(position=i,peaknumber=length(tc_NAR_ERCC[tc_NAR_ERCC[,8]==i,6])))
}
temp <- rbind(temp,data.frame(position=21,peaknumber=length(tc_NAR_ERCC[tc_NAR_ERCC[,8]>=21,6]))) 



ggplot(temp,aes(x=position,y=peaknumber)) +
  geom_bar(position=position_dodge(0.8),width=0.8,stat="identity") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 12),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(1,21),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",">20"))+
  #scale_y_continuous(limits = c(0,1000000),breaks = seq(0,1000000,250000),
  #                   labels = c("0","25","50","75","100"))+
  labs(x=element_blank())






temp <- data.frame()
for(i in 1:20) {
  temp <- rbind(temp,
                data.frame(position=i,peaknumber=length(tc_C1_CAGE_ERCC[tc_C1_CAGE_ERCC[,8]==i,6])))
}
temp <- rbind(temp,data.frame(position=21,peaknumber=length(tc_C1_CAGE_ERCC[tc_C1_CAGE_ERCC[,8]>=21,6]))) 



ggplot(temp,aes(x=position,y=peaknumber)) +
  geom_bar(position=position_dodge(0.8),width=0.8,stat="identity") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 12),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(1,21),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",">20"))+
  #scale_y_continuous(limits = c(0,1000000),breaks = seq(0,1000000,250000),
  #                   labels = c("0","25","50","75","100"))+
  labs(x=element_blank())






temp <- data.frame()
for(i in 1:20) {
  temp <- rbind(temp,
                data.frame(position=i,peaknumber=length(tc_NAR_ERCC[tc_NAR_ERCC[,8]==i,6])))
}
temp <- rbind(temp,data.frame(position=21,peaknumber=length(tc_NAR_ERCC[tc_NAR_ERCC[,8]>=21,6]))) 



ggplot(temp,aes(x=position,y=peaknumber)) +
  geom_bar(position=position_dodge(0.8),width=0.8,stat="identity") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 12),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(1,21),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",">20"))+
  #scale_y_continuous(limits = c(0,1000000),breaks = seq(0,1000000,250000),
  #                   labels = c("0","25","50","75","100"))+
  labs(x=element_blank())









### sufig1 h; sifig1 i


for(i in 1:nrow(DRG_tss)) {
  a <- DRG_tss_region[DRG_tss_region[,1]==DRG_tss[i,1] & DRG_tss_region[,3]==DRG_tss[i,3] & DRG_tss_region[,6]==DRG_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  DRG_tss[i,7] <- b
}


DRG_tss[is.na(DRG_tss[,7]),7] <- "intergenic"
DRG_tss[DRG_tss[,7]=="TES1000",7] <- "intergenic"
DRG_tss[DRG_tss[,7]=="last_exon",7] <- "other_exon"
DRG_tss[DRG_tss[,7]=="last_intron",7] <- "other intron"
DRG_tss[DRG_tss[,7]=="one_exon",7] <- "first_exon"
DRG_tss[DRG_tss[,7]=="one_intron",7] <- "first_intron"
DRG_tss[DRG_tss[,7]=="first_exon",7] <- "first exon"
DRG_tss[DRG_tss[,7]=="first_intron",7] <- "first intron"
DRG_tss[DRG_tss[,7]=="other_exon",7] <- "other exon"
DRG_tss[DRG_tss[,7]=="other_intron",7] <- "other intron"
DRG_tss[DRG_tss[,7]=="5UTR",7] <- "5'UTR"
DRG_tss[DRG_tss[,7]=="3UTR",7] <- "3'UTR"
DRG_tss[DRG_tss[,7]=="TSS1000",7] <- "TSS±1000"

DRG_tss_df <- as.data.frame(table(DRG_tss[,7]))
DRG_tss_df$Var1 <- factor(DRG_tss_df$Var1, levels=c("TSS±1000", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
DRG_tss_df <- DRG_tss_df[order(DRG_tss_df[,1]),]
DRG_tss_df[,1] <- paste(DRG_tss_df[,1]," (",round(DRG_tss_df[,2]/sum(DRG_tss_df[,2])*100,2),"%)",sep = "")
DRG_tss_df$Var1 <- factor(DRG_tss_df$Var1, levels=c(grep("TSS",DRG_tss_df[,1],value = T),
                                                    grep("5'UTR",DRG_tss_df[,1],value = T), 
                                                    grep("first exon",DRG_tss_df[,1],value = T),
                                                    grep("first intron",DRG_tss_df[,1],value = T),
                                                    grep("other exon",DRG_tss_df[,1],value = T),
                                                    grep("other intron",DRG_tss_df[,1],value = T),
                                                    grep("3'UTR",DRG_tss_df[,1],value = T),
                                                    grep("intergenic",DRG_tss_df[,1],value = T)), ordered=TRUE)
DRG_tss_df <- DRG_tss_df[order(DRG_tss_df[,1]),]

for(i in 1:nrow(DRG_tes)) {
  a <- DRG_tes_region[DRG_tes_region[,1]==DRG_tes[i,1] & DRG_tes_region[,3]==DRG_tes[i,3] & DRG_tes_region[,6]==DRG_tes[i,6],]
  b <- ordered(a[,13])
  b <- ordered(b, levels = c("TES1000", "3UTR", "last_exon","one_exon","last_intron","one_intron","Other_exon","first_exon","other_intron","first_intron","5UTR","TSS1000") )
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  DRG_tes[i,7] <- b
}

DRG_tes[is.na(DRG_tes[,7]),7] <- "intergenic"
DRG_tes[DRG_tes[,7]=="TSS1000",7] <- "intergenic"
DRG_tes[DRG_tes[,7]=="first_exon",7] <- "other exon"
DRG_tes[DRG_tes[,7]=="first_intron",7] <- "other intron"
DRG_tes[DRG_tes[,7]=="one_exon",7] <- "last exon"
DRG_tes[DRG_tes[,7]=="one_intron",7] <- "last intron"
DRG_tes[DRG_tes[,7]=="last_exon",7] <- "last exon"
DRG_tes[DRG_tes[,7]=="last_intron",7] <- "last intron"
DRG_tes[DRG_tes[,7]=="other_exon",7] <- "other exon"
DRG_tes[DRG_tes[,7]=="other_intron",7] <- "other intron"
DRG_tes[DRG_tes[,7]=="5UTR",7] <- "5'UTR"
DRG_tes[DRG_tes[,7]=="3UTR",7] <- "3'UTR"
DRG_tes[DRG_tes[,7]=="TES1000",7] <- "TES±1000"



DRG_tes_df <- as.data.frame(table(DRG_tes[,7]))
DRG_tes_df$Var1 <- factor(DRG_tes_df$Var1, levels=c("TES±1000", "3'UTR", "last exon","last intron","other exon","other intron","5'UTR","intergenic"), ordered=TRUE)
DRG_tes_df <- DRG_tes_df[order(DRG_tes_df[,1]),]
DRG_tes_df[,1] <- paste(DRG_tes_df[,1]," (",round(DRG_tes_df[,2]/sum(DRG_tes_df[,2])*100,2),"%)",sep = "")
DRG_tes_df$Var1 <- factor(DRG_tes_df$Var1, levels=c(grep("TES",DRG_tes_df[,1],value = T),
                                                    grep("3'UTR",DRG_tes_df[,1],value = T),
                                                    grep("last exon",DRG_tes_df[,1],value = T),
                                                    grep("last intron",DRG_tes_df[,1],value = T),
                                                    grep("other exon",DRG_tes_df[,1],value = T),
                                                    grep("other intron",DRG_tes_df[,1],value = T),
                                                    grep("5'UTR",DRG_tes_df[,1],value = T),
                                                    grep("intergenic",DRG_tes_df[,1],value = T)), ordered=TRUE)
DRG_tes_df <- DRG_tes_df[order(DRG_tes_df[,1]),]



ggplot(DRG_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+ 
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL)

ggplot(DRG_tes_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+ ## 
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL)

