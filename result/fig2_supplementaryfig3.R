# cost
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




load("fig2_supplementaryfig3.RData")

#sccat <- read.table("/home/zjw/zjw/nc/20200415cost/sccat.tsv")
#pacbio <- read.table("/home/zjw/zjw/nc/20200415cost/pacbio.tsv")
### fig2a


catseq <- data.frame(number=c(647.5,1225.5,1647.0,2020.5,2365.0,2650.5,2941.5,3191.5,3400.0,3577.0,3790.5,3982.0,4148.0,4374.5,4486.5,4690.5,4794.0,4973.0,5117.5,5194.5,5395.5,5488.5,5609.5,5726.0,5762.5,5912.0,6065.5,6185.5,6261.0,6352.5,6422.0,6509.5,6567.5,6644.0),
                     money=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4)/0.04384678*1.82)   #######新一批测序 0.04384678的转化率



y1 <- catseq[,1]
x1 <- catseq[,2]






###############pacbio


money <- c(0.32,0.35,0.89,0.24,1.03,3.52,0.22,1.31)*9000/20/6.87


pacbio <- data.frame(number=c(798,423,785,282,811,1200,356,958), 
                     money=c(0.32,0.35,0.89,0.24,1.03,3.52,0.22,1.31)*9000/20/6.87)

y2 <- pacbio[,1]
x2 <- pacbio[,2]



##############smart-seq2


smartseq <- data.frame(number=c(5.5,10,19.5,45,146,313,478.5,653.5,800.5,966.5,1087.5,1213),
                       money=c(0.64,1.28,2.56,5,10,20,30,40,50,60,70,80)/0.36*1.82)  


y3 <- smartseq[,1]
x3 <- smartseq[,2]

#### http://blog.sciencenet.cn/blog-651374-1126673.html
x <- c(1, 3, 6,  9,  13,   17)
y <- c(5, 8, 11, 13, 13.2, 13.5)
trendline(x, y, model="power2P", eDigit = 3, eSize = 1.4, text.col = "blue")

####x1x2x3 y1y2y3
trendline(x1, y1, model="power2P", ePos.x = "topleft", summary=TRUE, 
          eDigit=5,show.Rpvalue = FALSE,xlab = "cost",ylab = "number",xlim=c(0,500),ylim=c(0,7000))

par(new=TRUE)


trendline(x2, y2, model="power2P", ePos.x = "topleft", summary=TRUE, 
          eDigit=5,show.Rpvalue = FALSE,xlab = "cost",ylab = "number",xlim=c(0,500),ylim=c(0,7000))

par(new=TRUE)

trendline(x3, y3, model="power2P", ePos.x = "topleft", summary=TRUE, 
          eDigit=5,show.Rpvalue = FALSE,xlab = "cost",ylab = "number",xlim=c(0,500),ylim=c(0,7000))





### fig2b


df <- rbind(data.frame(group="scCAT-seq",value=sccat[,2]/sccat[,1]*1000),
            data.frame(group="ScISOr-seq",value=pacbio[,2]/pacbio[,1]*1000))

my_comparisons <- list(c("scCAT-seq", "ScISOr-seq"))

ggboxplot(df, x="group", y="value", fill = "group",ylab="Cost($) to cover 1000 transcript",xlab="",legend = "right",add.params = list(binwidth = 0.02)) +
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13)) +
  stat_compare_means(comparisons=my_comparisons,label = "p-value", method = "wilcox.test",paired = F) +
  scale_y_continuous(limits = c(0,250),expand = c(0,0))








#### fig2c


D3_4_5cap_dominant_tss_in_gene <- D3_4_5cap_dominant_tss_in_gene[,c(1,2,3,10,5,6)]
D3_4_3tail_dominant_tss_in_gene <- D3_4_3tail_dominant_tss_in_gene[,c(1,2,3,10,5,6)]



for(i in grep("in_gene",objects(),value = T)) {
  a <- paste(i,'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    b[1,5] <- sum(b[,5])
    a <- paste(i,'_major <- rbind(',i,'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}



#OC1_1 <- read.table("G:/CAGEr/CAGEr20190821/OC_1_1_OC_1_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
#OC1_2 <- read.table("G:/CAGEr/CAGEr20190821/OC_1_2_OC_1_2.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
#OC1_3 <- read.table("G:/CAGEr/CAGEr20190821/OC_1_3_OC_1_3.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
#DRG_1 <- read.table("G:/CAGEr/CAGEr20190821/DRG_1_DRG_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
#DRG_2 <- read.table("G:/CAGEr/CAGEr20190821/DRG_2_DRG_2.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")


counts <- data.frame(gene=OC1_1[,1],oc1_1=OC1_1[,2],oc1_2=OC1_2[,2],oc1_3=OC1_3[,2],drg_1=DRG_1[,2],drg_2=DRG_2[,2])
counts <- counts[-c(54147:54151),]
counts <- counts[!(counts[,2]==0 & counts[,3]==0 & counts[,4]==0 & counts[,5]==0 & counts[,6]==0),]




for(i in 1:nrow(counts)) {
  if(counts[i,1] %in% ENSM2symbol[,1]) counts[i,1] <- ENSM2symbol[ENSM2symbol[,1]==counts[i,1],2]
}

colnames(counts)[1] <- "gene"


#smart <- read.table("G:/CAGEr/CAGEr20190616/smart_seq2_single_cell_martix.txt",header = T)

smart <- data.frame(gene=rownames(smart),read_count=smart[,24])
smart[,2] <- smart[,2]/sum(smart[,2])*10^6
smart <- smart[1:54028,]
smart <- smart[!smart[,2]==0,]


tss_genelist <- D3_4_5cap_dominant_tss_in_gene_major[,4]
tes_genelist <- D3_4_3tail_dominant_tss_in_gene_major[,4]
iso_genelist <- counts[!counts[,2]==0,1]

df3 <- rbind(data.frame(smart[smart[,1] %in% tss_genelist,],V3="tss"),
             data.frame(smart[smart[,1] %in% tes_genelist,],V3="tes"),
             data.frame(smart[smart[,1] %in% iso_genelist,],V3="iso"))

df3[,2] <- log10(df3[,2])

colnames(df3) <- c("gene","V1","V2")

my_comparisons <- list(c("iso", "tes"),c("iso", "tss"))
ggviolin(df3, x="V2", y="V1", fill = "V2", palette = "jco",ylab="smart_RPM_logscale",xlab="",legend = "right") +   ###boxplot可改成violin
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  stat_compare_means(comparisons=my_comparisons,label = "p-value", method = "wilcox.test") 




my_comparisons <- list(c("iso", "tes"),c("iso", "tss"))
ggviolin(df3, x="V2", y="V1", fill = "V2", palette = "jco",ylab="smart_RPM_logscale",xlab="",legend = "right") +   ###boxplot可改成violin
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  stat_compare_means(comparisons=my_comparisons,label = "p-value", method = "wilcox.test") 









####fig2d oocyte

number <- c()


number <- c(number,nrow(D4tss[D4tss$isindsc==0 & D4tss$model.prediction==1,])+nrow(D4tss_intergenic[D4tss_intergenic$isindsc==0 & D4tss_intergenic$model.prediction==1,]))


number <- c(number,nrow(D4tes[D4tes$isinPAS==0 & D4tes$model.prediction==1,])+nrow(D4tes_intergenic[D4tes_intergenic$isinPAS==0 & D4tes_intergenic$model.prediction==1,]))



number <- c(number,nrow(D4_novel_gene))

number

novel <- data.frame(type=c("Novel TSS","Novel TES","Novel gene"),
                    number)

novel[,1] <- factor(novel[,1],levels = c("Novel TSS","Novel TES","Novel gene"))

ggplot(novel,aes(x=type,y=number,fill=type)) +
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
  scale_y_continuous(breaks = seq(0,6000,1000),limits = c(0,6200),expand=c(0,0))






###### fig2e DRG

number <- c()


number <- c(number,nrow(DRGtss[DRGtss$isindsc==0 & DRGtss$model.prediction==1,])+nrow(DRGtss_intergenic[DRGtss_intergenic$isindsc==0 & DRGtss_intergenic$model.prediction==1,]))


number <- c(number,nrow(DRGtes[DRGtes$isinPAS==0 & DRGtes$model.prediction==1,])+nrow(DRGtes_intergenic[DRGtes_intergenic$isinPAS==0 & DRGtes_intergenic$model.prediction==1,]))



number <- c(number,nrow(DRG_novel_gene))

number

novel <- data.frame(type=c("Novel TSS","Novel TES","Novel gene"),
                    number)

novel[,1] <- factor(novel[,1],levels = c("Novel TSS","Novel TES","Novel gene"))

ggplot(novel,aes(x=type,y=number,fill=type)) +
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
  scale_y_continuous(expand=c(0,0))














####supplementaryfig3a HEK293T


number <- c()


number <- c(number,nrow(HEK293Ttss[HEK293Ttss$isindsc==0 & HEK293Ttss$model.prediction==1,])+nrow(HEK293Ttss_intergenic[HEK293Ttss_intergenic$isindsc==0 & HEK293Ttss_intergenic$model.prediction==1,]))


number <- c(number,nrow(HEK293Ttes[HEK293Ttes$isinPAS==0 & HEK293Ttes$model.prediction==1,])+nrow(HEK293Ttes_intergenic[HEK293Ttes_intergenic$isinPAS==0 & HEK293Ttes_intergenic$model.prediction==1,]))



number <- c(number,nrow(HEK293T_novel_gene))

number

novel <- data.frame(type=c("Novel TSS","Novel TES","Novel gene"),
                    number)

novel[,1] <- factor(novel[,1],levels = c("Novel TSS","Novel TES","Novel gene"))

ggplot(novel,aes(x=type,y=number,fill=type)) +
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
  scale_y_continuous(expand=c(0,0))

















#### supplementaryfig3b hESC



number <- c()


number <- c(number,nrow(hESCtss[hESCtss$isindsc==0 & hESCtss$model.prediction==1,])+nrow(hESCtss_intergenic[hESCtss_intergenic$isindsc==0 & hESCtss_intergenic$model.prediction==1,]))


number <- c(number,nrow(hESCtes[hESCtes$isinPAS==0 & hESCtes$model.prediction==1,])+nrow(hESCtes_intergenic[hESCtes_intergenic$isinPAS==0 & hESCtes_intergenic$model.prediction==1,]))






number <- c(number,nrow(hESC_novel_gene))

number

novel <- data.frame(type=c("Novel TSS","Novel TES","Novel gene"),
                    number)

novel[,1] <- factor(novel[,1],levels = c("Novel TSS","Novel TES","Novel gene"))

ggplot(novel,aes(x=type,y=number,fill=type)) +
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
  scale_y_continuous(expand=c(0,0))




















####supplementaryfig7h organiod





number <- c()


number <- c(number,nrow(organiodtss[organiodtss$isindsc==0 & organiodtss$model.prediction==1,])+nrow(organiodtss_intergenic[organiodtss_intergenic$isindsc==0 & organiodtss_intergenic$model.prediction==1,]))


number <- c(number,nrow(organiodtes[organiodtes$isinPAS==0 & organiodtes$model.prediction==1,])+nrow(organiodtes_intergenic[organiodtes_intergenic$isinPAS==0 & organiodtes_intergenic$model.prediction==1,]))





number <- c(number,nrow(organiod_novel_gene))

number

novel <- data.frame(type=c("Novel TSS","Novel TES","Novel gene"),
                    number)

novel[,1] <- factor(novel[,1],levels = c("Novel TSS","Novel TES","Novel gene"))

ggplot(novel,aes(x=type,y=number,fill=type)) +
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
  scale_y_continuous(breaks = seq(0,6000,1000),expand=c(0,0))







