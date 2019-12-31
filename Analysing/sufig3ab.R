options(stringsAsFactors = FALSE)
options(scipen = 100)
library(CAGEr)
library(rlist)
library(ggplot2)
library(Rmisc)
library(scales)
library(basicTrendline)
library(broom)
library(dplyr)
library(nlstools)

setwd("G:/CAGEr/CAGEr20190704sensitivity_new/")   #####做对比图就选这个


##############一般从这里开始

rm(list = ls())



#for(i in grep("5bp",list.files(),value = T)) {
#  a <- paste(strsplit(i,split = "5bp.txt")[[1]][1],' <- read.table("',i,'")',sep = "")
#  print(a)
#  eval(parse(text=a))
#}


#for(i in grep("e2560000|e10000|e20000|e40000|e80000",objects(),value = T)) {
#  a <- paste('rm(',i,')',sep = "")
#  print(a)
#  eval(parse(text=a))
#}



###############这一步很重要，把转录本信息和基因信息加在一起
#for(i in grep("seed",objects(),value = T)) {
#  a <- paste(i,'[,2] <- paste(',i,'[,2],',i,'[,3])',sep = "")
#  print(a)
#  eval(parse(text=a))
#  a <- paste(i,' <- ',i,'[,c(1,2)]',sep = "")
#  print(a)
#  eval(parse(text=a))
#}

#####正负合并
#for(i in grep("plus",objects(),value = T)) {
#  for(j in grep("minus",objects(),value = T)) {
#    if(strsplit(i,split = "plus_")[[1]][2]==strsplit(j,split = "minus_")[[1]][2] & strsplit(i,split = "plus")[[1]][1]==strsplit(j,split = "minus")[[1]][1]) {
#      print(c(i,j))
#      a <- paste(strsplit(i,split = "_plus")[[1]][1],strsplit(i,split = "_plus")[[1]][2],' <- rbind(',i,',',j,')',sep = "")
#      print(a)
#      eval(parse(text=a))
#      a <- paste('rm(',i,')',sep = "")
#      print(a)
#      eval(parse(text=a))
#      a <- paste('rm(',j,')',sep = "")
#      print(a)
#      eval(parse(text=a))
#    }
#  }
#}



load("sufig3a.RData")
load("sufig3b.RData")






###############对单端处理求方便
for(i in grep("tss|tes",objects(),value = T)) {
  a <- paste(i,' <- ',i,'[,c(2,1,1)]',sep = "")
  print(a)
  eval(parse(text=a))
}




df <- data.frame()
for(i in grep("_sample_",objects(),value = T)) {
  a <- paste('b <- ',i,sep = "")
  eval(parse(text=a))
  df <- rbind(df,data.frame(number=nrow(b[b[,2]>=2 & b[,3]>=2,]),
                            size=strsplit(strsplit(i,split = "size")[[1]][2],split = "_")[[1]][1],
                            methods=paste(strsplit(i,split = "_sample_")[[1]][1],strsplit(i,split = "00_")[[1]][2],sep = "_")
  ))
}




df[,2] <- as.numeric(df[,2])

df[,2] <- log2(df[,2])



dfc <- summarySE(df, measurevar="number", groupvars=c("size","methods"))
dfc[,2] <- as.character(dfc[,2])

dfc <- dfc[dfc[,2] %in% c("CAT_TSS_tss","STRT_tss","SMART_TSS_tss","C1_CAGE_tss","NAR_tss"),]
dfc[,2] <- factor(dfc[,2],levels = c("CAT_TSS_tss","STRT_tss","SMART_TSS_tss","C1_CAGE_tss","NAR_tss"),ordered = T)

ggplot(dfc, aes(x=size, y=number, colour=methods)) +
  geom_errorbar(aes(ymin=number-se, ymax=number+se), width=0.1,size=0.8) +
  geom_line(size=1) +
  #geom_point()+
  scale_x_continuous(breaks = seq(log2(5000*2^1),log2(5000*2^8),1),
                     labels = c("10","20","40","80","160","320","640","1,280"))+
  scale_y_continuous(limits = c(0,20000),breaks = seq(0,20000,5000),
                     labels = c("0","5,000","10,000","15,000","20,000"))+
  scale_colour_manual(values=c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      breaks=c("CAT_TSS_tss","STRT_tss","SMART_TSS_tss","C1_CAGE_tss","NAR_tss"),
                      labels=c("CAT-seq","C1 STRT","Smart-seq2","C1 CAGE","NAR"))+ #改折线颜色
  labs(x="Number of reads uniquely mapped",y="Number of transcripts")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  #theme(panel.background = element_blank())+
  theme(panel.border = element_blank())+ #去除边框
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), #加上x轴
        axis.line.y=element_line(linetype=1,color="black",size=1), #加上y轴
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 10.5), #改y轴字体
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))



dfc <- summarySE(df, measurevar="number", groupvars=c("size","methods"))
dfc[,2] <- as.character(dfc[,2])
#dfc <- dfc[dfc[,2] %in% c("BAT_ESC_tes","CAT_TES_tes","SMART_TES_tes"),]
#dfc[,2] <- factor(dfc[,2],levels = c("CAT_TES_tes","BAT_ESC_tes","SMART_TES_tes"),ordered = T)

dfc <- dfc[dfc[,2] %in% c("BAT_tes","CAT_TES_tes","SMART_TES_tes"),]
dfc[,2] <- factor(dfc[,2],levels = c("CAT_TES_tes","BAT_tes","SMART_TES_tes"),ordered = T)

ggplot(dfc, aes(x=size, y=number, colour=methods)) +
  geom_errorbar(aes(ymin=number-se, ymax=number+se), width=0.1,size=0.8) +
  geom_line(size=1) +
  #geom_point()+
  scale_x_continuous(breaks = seq(log2(5000*2^1),log2(5000*2^8),1),
                     labels = c("10","20","40","80","160","320","640","1,280"))+
  scale_y_continuous(limits = c(0,20000),breaks = seq(0,20000,5000),
                     labels = c("0","5,000","10,000","15,000","20,000"))+
  scale_colour_manual(values=c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      breaks=c("CAT_TES_tes","BAT_tes","SMART_TES_tes"),
                      labels=c("CAT-seq","BAT-seq","Smart-seq2"))+ #改折线颜色  
  labs(x="Number of reads uniquely mapped",y="Number of transcripts")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  #theme(panel.background = element_blank())+
  theme(panel.border = element_blank())+ #去除边框
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), #加上x轴
        axis.line.y=element_line(linetype=1,color="black",size=1), #加上y轴
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 10.5), #改y轴字体
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))


#dfc <- summarySE(df, measurevar="number", groupvars=c("size","methods"))
#dfc[,2] <- as.character(dfc[,2])
dfc <- dfc[41:80,]
dfc <- dfc[dfc[,2] %in% c("CAT-seq simultaneous","Smart-seq2 simultaneous"),]
ggplot(dfc, aes(x=size, y=number, colour=methods)) +
  geom_errorbar(aes(ymin=number-se, ymax=number+se), width=0.1,size=0.8) +
  geom_line(size=1) +
  #geom_point()+
  scale_x_continuous(breaks = seq(log2(5000*2^1),log2(5000*2^8),1),
                     labels = c("0.01M","0.02M","0.04M","0.08M","0.16M","0.32M","0.64M","1.28M"))+
  #scale_y_continuous(limits = c(0,7000))+
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+ #改折线颜色
  labs(x="reads of uniquely map",y="number of end detected")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  #theme(panel.background = element_blank())+
  theme(panel.border = element_blank())+ #去除边框
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), #加上x轴
        axis.line.y=element_line(linetype=1,color="black",size=1), #加上y轴
        axis.text.x = element_text(face = "plain",size = 10.5,angle=45,hjust = 1,vjust = 1), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 10.5), #改y轴字体
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.75,0.25),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))
