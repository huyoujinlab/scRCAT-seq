setwd("C:/Users/zhong/Desktop/20190712peak/")

#DRG_tss <- read.table("C:/Users/zhong/Desktop/D_dominant_tss_in_gene_1.bed")
#DRG_tes <- read.table("C:/Users/zhong/Desktop/D_dominant_tes_in_gene_1.bed")


DRG_tss_region <- read.table("C:/Users/zhong/Desktop/D_dominant_tss_in_gene_1_region.bed",sep = "\t")
DRG_tes_region <- read.table("C:/Users/zhong/Desktop/D_dominant_tes_in_gene_1_region.bed",sep = "\t")

load("pie_plot.RData")

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
  coord_polar(theta = "y")+ ## 以y轴建立极坐标
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +#去除边框
  theme(axis.text.x = element_blank(), #改x轴字体
        axis.text.y = element_blank(), #改y轴字体
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL)

ggplot(DRG_tes_df, aes(x = 1, weight = Freq, fill = Var1)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+ ## 以y轴建立极坐标
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +#去除边框
  theme(axis.text.x = element_blank(), #改x轴字体
        axis.text.y = element_blank(), #改y轴字体
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL)
