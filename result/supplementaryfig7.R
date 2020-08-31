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




load("supplementaryfig10.RData")


####supplementaryfig7 g 10x tss

cell <- c("ARPE","organoid","HEK293T","hESC","mESC")

for (i in cell) {
  if (i=="ARPE") {novel_tss <- data.frame()}
  a <- paste('temp1 <- nrow(',i,'[',i,'$isindsc==0 & ',i,'$model.prediction==1,])',sep="")
  eval(parse(text=a))
  print(a)
  a <- paste('temp2 <- nrow(',i,'_intergenic[',i,'_intergenic$isindsc==0 & ',i,'_intergenic$model.prediction==1,])',sep="")
  eval(parse(text=a))
  print(a)
  novel_tss <- rbind(novel_tss,
                     data.frame(cell=i,number=sum(temp1,temp2)))
}


novel_tss[,1] <- factor(novel_tss[,1],levels = c("HEK293T","hESC","ARPE","organoid","mESC"))




ggplot(novel_tss,aes(x=cell,y=number,fill=cell)) +
  geom_bar(position=position_dodge(0.7),width=0.5,stat="identity") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  #theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="none") +
  labs(x=element_blank(),y="Number")+
  scale_y_continuous(breaks = seq(0,16000,2000),limits = c(0,16000),expand=c(0,0))

