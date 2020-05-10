# supplementary2
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





## supplementary2a supplementary2b

load("supplementary2.RData")

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

dfc <- dfc[dfc[,2] %in% c("CAT_TSS_tss","STRT_tss","SMART_TSS_tss","C1_CAGE_tss","hESCsccattss_tss"),]
dfc[,2] <- factor(dfc[,2],levels = c("CAT_TSS_tss","STRT_tss","SMART_TSS_tss","C1_CAGE_tss","hESCsccattss_tss"),ordered = T)

ggplot(dfc, aes(x=size, y=number, colour=methods)) +
  geom_errorbar(aes(ymin=number-se, ymax=number+se), width=0.1,size=0.8) +
  geom_line(size=1) +
  #geom_point()+
  scale_x_continuous(breaks = seq(log2(5000*2^1),log2(5000*2^8),1),
                     labels = c("10","20","40","80","160","320","640","1,280"))+
  scale_y_continuous(limits = c(0,25000),breaks = seq(0,25000,5000),)+
  scale_colour_manual(values=c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      breaks=c("CAT_TSS_tss","STRT_tss","SMART_TSS_tss","C1_CAGE_tss","NAR_tss","hESCsccattss_tss"),
                      labels=c("CAT-seq_DRG","C1 STRT","Smart-seq2","C1 CAGE","NAR","CAT-seq_hESC"))+
  labs(x="Number of reads uniquely mapped",y="Number of transcripts")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  #theme(panel.background = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 10.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))



dfc <- summarySE(df, measurevar="number", groupvars=c("size","methods"))
dfc[,2] <- as.character(dfc[,2])
#dfc <- dfc[dfc[,2] %in% c("BAT_ESC_tes","CAT_TES_tes","SMART_TES_tes"),]
#dfc[,2] <- factor(dfc[,2],levels = c("CAT_TES_tes","BAT_ESC_tes","SMART_TES_tes"),ordered = T)

dfc <- dfc[dfc[,2] %in% c("BAT_tes","CAT_TES_tes","SMART_TES_tes","hESCsccattes_tes"),]
dfc[,2] <- factor(dfc[,2],levels = c("CAT_TES_tes","BAT_tes","SMART_TES_tes","hESCsccattes_tes"),ordered = T)

ggplot(dfc, aes(x=size, y=number, colour=methods)) +
  geom_errorbar(aes(ymin=number-se, ymax=number+se), width=0.1,size=0.8) +
  geom_line(size=1) +
  #geom_point()+
  scale_x_continuous(breaks = seq(log2(5000*2^1),log2(5000*2^8),1),
                     labels = c("10","20","40","80","160","320","640","1,280"))+
  scale_y_continuous(limits = c(0,25000),breaks = seq(0,25000,5000),)+
  scale_colour_manual(values=c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      breaks=c("CAT_TES_tes","BAT_tes","SMART_TES_tes","hESCsccattes_tes"),
                      labels=c("CAT-seq_DRG","BAT-seq","Smart-seq2","CAT-seq_hESC"))+
  labs(x="Number of reads uniquely mapped",y="Number of transcripts")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  #theme(panel.background = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 10.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))
















### supplementary2 c

rm(list=ls())

load("supplementary2.RData")


#P42_51_dominant_tss_in_gene <- cbind(read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/20190831P42_51_5cap_peak_with_prediction.csv"),
#                                     read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/tc_P42_51_5cap_peak.csvnew.csv")[,7])

#P44_51_dominant_tss_in_gene <- cbind(read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/20190831P44_51_5cap_peak_with_prediction.csv"),
#                                     read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/tc_P44_51_5cap_peak.csvnew.csv")[,7])

#P44_71_dominant_tss_in_gene <- cbind(read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/20190831P44_71_5cap_peak_with_prediction.csv"),
#                                     read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/tc_P44_71_5cap_peak.csvnew.csv")[,7])


#P42_51_dominant_tes_in_gene <- cbind(read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/20190831P42_51_3tail_peak_with_prediction.csv"),
#                                     read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/tc_P42_51_3tail_peak.csvnew.csv")[,7])

#P44_51_dominant_tes_in_gene <- cbind(read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/20190831P44_51_3tail_peak_with_prediction.csv"),
#                                     read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/tc_P44_51_3tail_peak.csvnew.csv")[,7])

#P44_71_dominant_tes_in_gene <- cbind(read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/20190831P44_71_3tail_peak_with_prediction.csv"),
#                                     read.csv("C:/Users/zhong/Desktop/20190831peak_outputs-20190908T082631Z-001/20190831peak_outputs/tc_P44_71_3tail_peak.csvnew.csv")[,7])

#P42_51_dominant_tss_in_gene <- P42_51_dominant_tss_in_gene[P42_51_dominant_tss_in_gene[,16]==1,]
#P44_51_dominant_tss_in_gene <- P44_51_dominant_tss_in_gene[P44_51_dominant_tss_in_gene[,16]==1,]
#P44_71_dominant_tss_in_gene <- P44_71_dominant_tss_in_gene[P44_71_dominant_tss_in_gene[,16]==1,]

#P42_51_dominant_tes_in_gene <- P42_51_dominant_tes_in_gene[P42_51_dominant_tes_in_gene[,16]==1,]
#P44_51_dominant_tes_in_gene <- P44_51_dominant_tes_in_gene[P44_51_dominant_tes_in_gene[,16]==1,]
#P44_71_dominant_tes_in_gene <- P44_71_dominant_tes_in_gene[P44_71_dominant_tes_in_gene[,16]==1,]



#P42_51_dominant_tss_in_gene <- P42_51_dominant_tss_in_gene[,c(3,17,2,7,6)]
#P44_51_dominant_tss_in_gene <- P44_51_dominant_tss_in_gene[,c(3,17,2,7,6)]
#P44_71_dominant_tss_in_gene <- P44_71_dominant_tss_in_gene[,c(3,17,2,7,6)]

#P42_51_dominant_tes_in_gene <- P42_51_dominant_tes_in_gene[,c(3,17,2,7,6)]
#P44_51_dominant_tes_in_gene <- P44_51_dominant_tes_in_gene[,c(3,17,2,7,6)]
#P44_71_dominant_tes_in_gene <- P44_71_dominant_tes_in_gene[,c(3,17,2,7,6)]

P42_51_dominant_tss_in_gene <- data.frame(chr=P42_51_dominant_tss_in_gene[,1],
                                          start=P42_51_dominant_tss_in_gene[,2]-1,
                                          end=P42_51_dominant_tss_in_gene[,2],
                                          P42_51_dominant_tss_in_gene[,c(3,4,5)])

P44_51_dominant_tss_in_gene <- data.frame(chr=P44_51_dominant_tss_in_gene[,1],
                                          start=P44_51_dominant_tss_in_gene[,2]-1,
                                          end=P44_51_dominant_tss_in_gene[,2],
                                          P44_51_dominant_tss_in_gene[,c(3,4,5)])

P44_71_dominant_tss_in_gene <- data.frame(chr=P44_71_dominant_tss_in_gene[,1],
                                          start=P44_71_dominant_tss_in_gene[,2]-1,
                                          end=P44_71_dominant_tss_in_gene[,2],
                                          P44_71_dominant_tss_in_gene[,c(3,4,5)])


P42_51_dominant_tes_in_gene <- data.frame(chr=P42_51_dominant_tes_in_gene[,1],
                                          start=P42_51_dominant_tes_in_gene[,2]-1,
                                          end=P42_51_dominant_tes_in_gene[,2],
                                          P42_51_dominant_tes_in_gene[,c(3,4,5)])

P44_51_dominant_tes_in_gene <- data.frame(chr=P44_51_dominant_tes_in_gene[,1],
                                          start=P44_51_dominant_tes_in_gene[,2]-1,
                                          end=P44_51_dominant_tes_in_gene[,2],
                                          P44_51_dominant_tes_in_gene[,c(3,4,5)])

P44_71_dominant_tes_in_gene <- data.frame(chr=P44_71_dominant_tes_in_gene[,1],
                                          start=P44_71_dominant_tes_in_gene[,2]-1,
                                          end=P44_71_dominant_tes_in_gene[,2],
                                          P44_71_dominant_tes_in_gene[,c(3,4,5)])


P42_51 <- merge(P42_51_dominant_tss_in_gene,P42_51_dominant_tes_in_gene,by = "gene")
P44_51 <- merge(P44_51_dominant_tss_in_gene,P44_51_dominant_tes_in_gene,by = "gene")
P44_71 <- merge(P44_71_dominant_tss_in_gene,P44_71_dominant_tes_in_gene,by = "gene")



for(i in 1:nrow(P42_51)) {
  if(P42_51[i,6]=="+" & P42_51[i,3]>P42_51[i,9]) P42_51[i,11] <- "."
  if(P42_51[i,6]=="-" & P42_51[i,3]<P42_51[i,9]) P42_51[i,11] <- "."
  if(length(grep("inter",P42_51[i,1]))==1 & abs(P42_51[i,4]-P42_51[i,9])>100000) P42_51[i,11] <- "."
}
P42_51 <- P42_51[!P42_51[,11]==".",]

for(i in 1:nrow(P44_51)) {
  if(P44_51[i,6]=="+" & P44_51[i,4]>P44_51[i,9]) P44_51[i,11] <- "."
  if(P44_51[i,6]=="-" & P44_51[i,4]<P44_51[i,9]) P44_51[i,11] <- "."
  if(length(grep("inter",P44_51[i,1]))==1 & abs(P44_51[i,4]-P44_51[i,9])>100000) P44_51[i,11] <- "."
}
P44_51 <- P44_51[!P44_51[,11]==".",]

for(i in 1:nrow(P44_71)) {
  if(P44_71[i,6]=="+" & P44_71[i,3]>P44_71[i,9]) P44_71[i,11] <- "."
  if(P44_71[i,6]=="-" & P44_71[i,3]<P44_71[i,9]) P44_71[i,11] <- "."
  if(length(grep("inter",P44_71[i,1]))==1 & abs(P44_71[i,4]-P44_71[i,9])>100000) P44_71[i,11] <- "."
}
P44_71 <- P44_71[!P44_71[,11]==".",]




#### gene number
length(unique(P42_51[,1]))
length(unique(P44_51[,1]))
length(unique(P44_71[,1]))



P42_51_temp <- data.frame()
for(i in unique(P42_51[,1])) {
  temp <- P42_51[P42_51[,1]==i,]
  P42_51_temp <- rbind(P42_51_temp,
                       data.frame(V1=i,v2=max(c(length(unique(temp[,4])),length(unique(temp[,9]))))))
}


P44_51_temp <- data.frame()
for(i in unique(P44_51[,1])) {
  temp <- P44_51[P44_51[,1]==i,]
  P44_51_temp <- rbind(P44_51_temp,
                       data.frame(V1=i,v2=max(c(length(unique(temp[,4])),length(unique(temp[,9]))))))
}



P44_71_temp <- data.frame()
for(i in unique(P44_71[,1])) {
  temp <- P44_71[P44_71[,1]==i,]
  P44_71_temp <- rbind(P44_71_temp,
                       data.frame(V1=i,v2=max(c(length(unique(temp[,4])),length(unique(temp[,9]))))))
}









OC_1_in_gene <- unique(OC_1_in_gene)
OC_2_in_gene <- unique(OC_2_in_gene)
OC_3_in_gene <- unique(OC_3_in_gene)


OC_1_in_gene_temp <- as.data.frame(table(OC_1_in_gene[,10]),stringsAsFactors = F)
OC_2_in_gene_temp <- as.data.frame(table(OC_2_in_gene[,10]),stringsAsFactors = F)
OC_3_in_gene_temp <- as.data.frame(table(OC_3_in_gene[,10]),stringsAsFactors = F)


colnames(P42_51_temp) <- c("Var1","Freq")
colnames(P44_51_temp) <- c("Var1","Freq")
colnames(P44_71_temp) <- c("Var1","Freq")


df <- rbind(data.frame(P42_51_temp,type="scCAT-seq_sc#1"),
            data.frame(P44_51_temp,type="scCAT-seq_sc#2"),
            data.frame(P44_71_temp,type="scCAT-seq_sc#3"),
            data.frame(OC_1_in_gene_temp,type="ISO-seq_sc#1"),
            data.frame(OC_2_in_gene_temp,type="ISO-seq_sc#2"),
            data.frame(OC_3_in_gene_temp,type="ISO-seq_sc#3"))


df <- df[df[,1] %in% as.data.frame(table(df[,1]),stringsAsFactors = F)[as.data.frame(table(df[,1]),stringsAsFactors = F)[,2]==6,1],]






df <- rbind(data.frame(type="gene",group=c("scCAT-seq","scCAT-seq","scCAT-seq","Iso-seq","Iso-seq","Iso-seq"),number=c(length(unique(P42_51[,1])),length(unique(P44_51[,1])),length(unique(P44_71[,1])),length(unique(OC_1_in_gene[,10])),length(unique(OC_2_in_gene[,10])),length(unique(OC_3_in_gene[,10])))),
            data.frame(type="transcript",group=c("scCAT-seq","scCAT-seq","scCAT-seq","Iso-seq","Iso-seq","Iso-seq"),number=c(nrow(P42_51),nrow(P44_51),nrow(P44_71),length(OC_1_in_gene[,10]),length(OC_2_in_gene[,10]),length(OC_3_in_gene[,10]))))
df[,2] <- factor(df[,2],levels = c("scCAT-seq","Iso-seq"))

ggbarplot(df, x="type", y="number", add = "mean_se", fill = "group",color = "group",add.params = list(color = "black"),
          palette = "jco", position = position_dodge(0.8))+ 
  stat_compare_means(aes(group=group), label = "..p.format..", label.y = 9000,method = "t.test") +
  scale_y_continuous(limits = c(0,10000),expand=c(0,0),breaks = seq(0,10000,2000),labels = c("0","2,000","4,000","6,000","8,000","10,000")) +
  scale_fill_manual(values=c("#E69F00", "#999999"))+
  scale_color_manual(values=c("#E69F00", "#999999"))


















### supplementary2d
rm(list=ls())

load("supplementary2.RData")


#D3_tss <- cbind(read.csv("C:/Users/zhong/Desktop/201907211815/20190711D3tss_peak_new.csv"),
#                read.csv("C:/Users/zhong/Desktop/201907211815/majority_vote_tss_test_D3.csv")[,c(12,13)])
#D3_tss <- D3_tss[D3_tss[,17]==1,]
#D3_tss <- D3_tss[,c(2,3,4,1,6,5)]

#D3_tes <- cbind(read.csv("C:/Users/zhong/Desktop/201907211815/20190711D3tes_peak_new.csv"),
#                read.csv("C:/Users/zhong/Desktop/201907211815/majority_vote_tes_test_D3.csv")[,c(12,13)])
#D3_tes <- D3_tes[D3_tes[,17]==1,]
#D3_tes <- D3_tes[,c(2,3,4,1,6,5)]

######## OC 1 gene list 
#OC_in_gene <- read.table("iso_seq_OC_1_in_gene.bed")    
#iso_seq_OC_gene_list <- unique(OC_in_gene[,10])

####combine all peak


for(i in grep("D3_tss|D3_tes",objects(),value = T)) {
  a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    b[1,5] <- sum(b[,5])
    a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- rbind(',strsplit(i,split = "_3tail")[[1]][1],'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}

D3_tss_major <- D3_tss_major[order(D3_tss_major[,5],decreasing = T),]
D3_tss_major[1:3208,7] <- "high"
D3_tss_major[3209:6416,7] <- "medium"
D3_tss_major[6417:9626,7] <- "low"


D3_tes_major <- D3_tes_major[order(D3_tes_major[,5],decreasing = T),]
D3_tes_major[1:4200,7] <- "high"
D3_tes_major[4201:8400,7] <- "medium"
D3_tes_major[8400:12602,7] <- "low"




length(unique(intersect(iso_seq_OC_gene_list,D3_tss_major[D3_tss_major[,7]=="high",4])))
length(unique(intersect(iso_seq_OC_gene_list,D3_tss_major[D3_tss_major[,7]=="medium",4])))
length(unique(intersect(iso_seq_OC_gene_list,D3_tss_major[D3_tss_major[,7]=="low",4])))



length(unique(intersect(iso_seq_OC_gene_list,D3_tes_major[D3_tes_major[,7]=="high",4])))
length(unique(intersect(iso_seq_OC_gene_list,D3_tes_major[D3_tes_major[,7]=="medium",4])))
length(unique(intersect(iso_seq_OC_gene_list,D3_tes_major[D3_tes_major[,7]=="low",4])))


D3_tss_major[D3_tss_major[,4] %in% iso_seq_OC_gene_list,8] <- "detected_by_isoseq"
D3_tss_major[!D3_tss_major[,4] %in% iso_seq_OC_gene_list,8] <- "not_detected_by_isoseq"


D3_tes_major[D3_tes_major[,4] %in% iso_seq_OC_gene_list,8] <- "detected_by_isoseq"
D3_tes_major[!D3_tes_major[,4] %in% iso_seq_OC_gene_list,8] <- "not_detected_by_isoseq"



#TES  quantification
df <- data.frame()

df <- rbind(df,data.frame(RPM="1-10",type="not_detected",number=nrow(D3_tes_major[D3_tes_major[,5]>1 & D3_tes_major[,5]<10 & D3_tes_major[,8]=="not_detected_by_isoseq",])))
df <- rbind(df,data.frame(RPM="10-100",type="not_detected",number=nrow(D3_tes_major[D3_tes_major[,5]>10 & D3_tes_major[,5]<100 & D3_tes_major[,8]=="not_detected_by_isoseq",])))
df <- rbind(df,data.frame(RPM="100+",type="not_detected",number=nrow(D3_tes_major[D3_tes_major[,5]>100 & D3_tes_major[,8]=="not_detected_by_isoseq",])))

df <- rbind(df,data.frame(RPM="1-10",type="detected",number=nrow(D3_tes_major[D3_tes_major[,5]>1 & D3_tes_major[,5]<10 & D3_tes_major[,8]=="detected_by_isoseq",])))
df <- rbind(df,data.frame(RPM="10-100",type="detected",number=nrow(D3_tes_major[D3_tes_major[,5]>10 & D3_tes_major[,5]<100 & D3_tes_major[,8]=="detected_by_isoseq",])))
df <- rbind(df,data.frame(RPM="100+",type="detected",number=nrow(D3_tes_major[D3_tes_major[,5]>100  & D3_tes_major[,8]=="detected_by_isoseq",])))




df$RPM = factor(df$RPM, levels=c('1-10','10-100','100+'))
df$type = factor(df$type, levels=c('not_detected','detected'))

ggplot(data=df,aes(x=RPM,y=number,fill=type))+
  geom_bar(stat="identity",position="fill",width = 0.85,color = "black",size=0.7)+  
  labs(x="Gene expression level in scCAT-seq (TES RPM)",y="Fraction of gene",title="")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 10,angle=35,hjust = 1,vjust = 1), 
        axis.text.y = element_text(face = "plain",size = 10),
        legend.title = element_blank())+
  
  scale_fill_manual(values=c("#C7CACA","#7A3430"),
                    labels=c("Not detected in Iso-seq","Detected in Iso-seq"),
                    breaks=c("not_detected","detected")) 










### supplementary2 g
#load("supplementary2.RData")




#OC_1_in_gene <- read.table("iso_seq_OC_1_1_in_gene.bed")
#OC_2_in_gene <- read.table("iso_seq_OC_1_2_in_gene.bed")
#OC_3_in_gene <- read.table("iso_seq_OC_1_3_in_gene.bed")




#P42_51_5cap_in_gene <- read.table("P42_51_5cap_dominant_tes_ustream2k_and_genebody.bed")
#P44_51_5cap_in_gene <- read.table("P44_51_5cap_dominant_tes_ustream2k_and_genebody.bed")
#P44_71_5cap_in_gene <- read.table("P44_71_5cap_dominant_tes_ustream2k_and_genebody.bed")
#P42_51_3tail_in_gene <- read.table("P42_51_3tail_dominant_tes_genebody_and_downstream2k.bed")
#P44_51_3tail_in_gene <- read.table("P44_51_3tail_dominant_tes_genebody_and_downstream2k.bed")
#P44_71_3tail_in_gene <- read.table("P44_71_3tail_dominant_tes_genebody_and_downstream2k.bed")






P42_51_5cap_in_gene <- P42_51_5cap_in_gene[,c(1,2,3,10,5,6)]
P44_51_5cap_in_gene <- P44_51_5cap_in_gene[,c(1,2,3,10,5,6)]
P44_71_5cap_in_gene <- P44_71_5cap_in_gene[,c(1,2,3,10,5,6)]
P42_51_3tail_in_gene <- P42_51_3tail_in_gene[,c(1,2,3,10,5,6)]
P44_51_3tail_in_gene <- P44_51_3tail_in_gene[,c(1,2,3,10,5,6)]
P44_71_3tail_in_gene <- P44_71_3tail_in_gene[,c(1,2,3,10,5,6)]




P42_51_5cap <- P42_51_5cap_in_gene
P44_51_5cap <- P44_51_5cap_in_gene
P44_71_5cap <- P44_71_5cap_in_gene

P42_51_3tail <- P42_51_3tail_in_gene
P44_51_3tail <- P44_51_3tail_in_gene
P44_71_3tail <- P44_71_3tail_in_gene




P42_51 <- merge(P42_51_5cap,P42_51_3tail,by = "V10")
P44_51 <- merge(P44_51_5cap,P44_51_3tail,by = "V10")
P44_71 <- merge(P44_71_5cap,P44_71_3tail,by = "V10")

for(i in 1:nrow(P42_51)) {
  if(P42_51[i,6]=="+" & P42_51[i,3]>P42_51[i,9]) P42_51[i,11] <- "."
  if(P42_51[i,6]=="-" & P42_51[i,3]<P42_51[i,9]) P42_51[i,11] <- "."
  if(length(grep("inter",P42_51[i,1]))==1 & abs(P42_51[i,4]-P42_51[i,9])>100000) P42_51[i,11] <- "."
}
P42_51 <- P42_51[!P42_51[,11]==".",]

for(i in 1:nrow(P44_51)) {
  if(P44_51[i,6]=="+" & P44_51[i,4]>P44_51[i,9]) P44_51[i,11] <- "."
  if(P44_51[i,6]=="-" & P44_51[i,4]<P44_51[i,9]) P44_51[i,11] <- "."
  if(length(grep("inter",P44_51[i,1]))==1 & abs(P44_51[i,4]-P44_51[i,9])>100000) P44_51[i,11] <- "."
}
P44_51 <- P44_51[!P44_51[,11]==".",]

for(i in 1:nrow(P44_71)) {
  if(P44_71[i,6]=="+" & P44_71[i,3]>P44_71[i,9]) P44_71[i,11] <- "."
  if(P44_71[i,6]=="-" & P44_71[i,3]<P44_71[i,9]) P44_71[i,11] <- "."
  if(length(grep("inter",P44_71[i,1]))==1 & abs(P44_71[i,4]-P44_71[i,9])>100000) P44_71[i,11] <- "."
}
P44_71 <- P44_71[!P44_71[,11]==".",]




P42_51_new <- data.frame()
for(i in 1:nrow(P42_51)) {
  if(P42_51[i,6]=="+") {
    P42_51_new <- rbind(P42_51_new,
                        data.frame(V1=P42_51[i,2],
                                   V2=P42_51[i,3],
                                   V3=P42_51[i,9],
                                   P42_51[i,c(1,1,6)]))}
  if(P42_51[i,6]=="-") {
    P42_51_new <- rbind(P42_51_new,
                        data.frame(V1=P42_51[i,2],
                                   V2=P42_51[i,8],
                                   V3=P42_51[i,4],
                                   P42_51[i,c(1,1,6)]))}
}



P44_51_new <- data.frame()
for(i in 1:nrow(P44_51)) {
  if(P44_51[i,6]=="+") {
    P44_51_new <- rbind(P44_51_new,
                        data.frame(V1=P44_51[i,2],
                                   V2=P44_51[i,3],
                                   V3=P44_51[i,9],
                                   P44_51[i,c(1,1,6)]))}
  if(P44_51[i,6]=="-") {
    P44_51_new <- rbind(P44_51_new,
                        data.frame(V1=P44_51[i,2],
                                   V2=P44_51[i,8],
                                   V3=P44_51[i,4],
                                   P44_51[i,c(1,1,6)]))}
}


P44_71_new <- data.frame()
for(i in 1:nrow(P44_71)) {
  if(P44_71[i,6]=="+") {
    P44_71_new <- rbind(P44_71_new,
                        data.frame(V1=P44_71[i,2],
                                   V2=P44_71[i,3],
                                   V3=P44_71[i,9],
                                   P44_71[i,c(1,1,6)]))}
  if(P44_71[i,6]=="-") {
    P44_71_new <- rbind(P44_71_new,
                        data.frame(V1=P44_71[i,2],
                                   V2=P44_71[i,8],
                                   V3=P44_71[i,4],
                                   P44_71[i,c(1,1,6)]))}
}






df <- rbind(data.frame(V1=c(length(intersect(P42_51[,1],P44_51[,1]))/length(union(P42_51[,1],P44_51[,1])),
                            length(intersect(P42_51[,1],P44_71[,1]))/length(union(P42_51[,1],P44_71[,1])),
                            length(intersect(P44_51[,1],P44_71[,1]))/length(union(P44_51[,1],P44_71[,1]))),
                       V2="cat"),
            data.frame(V1=c(length(intersect(OC_1_in_gene[,10],OC_2_in_gene[,10]))/length(union(OC_1_in_gene[,10],OC_2_in_gene[,10])),
                            length(intersect(OC_1_in_gene[,10],OC_3_in_gene[,10]))/length(union(OC_1_in_gene[,10],OC_3_in_gene[,10])),
                            length(intersect(OC_2_in_gene[,10],OC_3_in_gene[,10]))/length(union(OC_2_in_gene[,10],OC_3_in_gene[,10]))),
                       V2="iso")
            
)


df[,2] <- as.character(df[,2])

my_comparisons <- list(c("iso", "cat"))

ggbarplot(df, x="V2", y="V1", add = "mean_se", fill = "V2",color = "V2",add.params = list(color = "black"),
          palette = "jco", position = position_dodge(0.8),ylab = "overlap rate")+ 
  stat_compare_means(comparisons=my_comparisons, label = "p-value",method = "t.test") +
  scale_y_continuous(limits = c(0,0.8),expand=c(0,0),breaks = seq(0,0.8,0.2),labels = c("0","0.2","0.4","0.6","0.8")) +
  scale_fill_manual(values=c("#E69F00", "#999999"))+
  scale_color_manual(values=c("#E69F00", "#999999"))

