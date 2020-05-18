# supplementaryfig1
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





### supplementaryfig1 c

load("supplementaryfig1.RData")


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



