options(stringsAsFactors = FALSE)
options(scipen = 100)

library(basicTrendline)
library(broom)
library(BuenColors)
library(data.table)
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

args = commandArgs(T)

##human

hESC_tss <- read.csv(args[1])
temp <-  read.csv(args[2])

dim(hESC_tss)
dim(temp)

hESC_tss <- rbind(hESC_tss,temp)
hESC_tss <- hESC_tss[hESC_tss$model.prediction==1,]
hESC_tss <- data.frame(hESC_tss[,2],hESC_tss[,7]-1,hESC_tss[,7],hESC_tss[c(1,6,5)])

hESC_tes <- read.csv(args[3])
temp <-  read.csv(args[4])
outdir <- args[5]
if (!dir.exists(outdir)) dir.create(outdir)

dim(hESC_tes)
dim(temp)

hESC_tes <- rbind(hESC_tes,temp)
hESC_tes <- hESC_tes[hESC_tes$model.prediction==1,]
hESC_tes <- data.frame(hESC_tes[,2],hESC_tes[,7]-1,hESC_tes[,7],hESC_tes[c(1,6,5)])


write.table(hESC_tss,"outdir/novel/hESC_tss.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(hESC_tes,"outdir/novel/hESC_tes.bed",quote = F,row.names = F,col.names = F,sep = "\t")


###bedtools intersect




system("bedtools intersect -s -a outdir/novel/hESC_tss.bed -b reference/gencode_hg38_all_gene_all_transcript_region.bed -wa -wb > outdir/novel/hESC_tss_region.bed")
system("bedtools intersect -s -a outdir/novel/hESC_tes.bed -b reference/gencode_hg38_all_gene_all_transcript_region.bed -wa -wb > outdir/novel/hESC_tes_region.bed")


hESC_tss_region <- read.table("outdir/novel/hESC_tss_region.bed",sep = "\t")

hESC_tes_region <- read.table("outdir/novel/hESC_tes_region.bed",sep = "\t")

for(i in 1:nrow(hESC_tss)) {
  a <- hESC_tss_region[hESC_tss_region[,1]==hESC_tss[i,1] & hESC_tss_region[,3]==hESC_tss[i,3] & hESC_tss_region[,6]==hESC_tss[i,6],]
  b <- unique(ordered(a[,13]))
  b <- ordered(b, levels = c("TSS1000", "5UTR", "first_exon","one_exon","first_intron","one_intron","Other_exon","last_exon","other_intron","last_intron","3UTR","TES1000"))
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  hESC_tss[i,7] <- b
}
hESC_tss[is.na(hESC_tss[,7]),7] <- "intergenic"
hESC_tss[hESC_tss[,7]=="TES1000",7] <- "intergenic"
hESC_tss[hESC_tss[,7]=="last_exon",7] <- "other_exon"
hESC_tss[hESC_tss[,7]=="last_intron",7] <- "other intron"
hESC_tss[hESC_tss[,7]=="one_exon",7] <- "first_exon"
hESC_tss[hESC_tss[,7]=="one_intron",7] <- "first_intron"
hESC_tss[hESC_tss[,7]=="first_exon",7] <- "first exon"
hESC_tss[hESC_tss[,7]=="first_intron",7] <- "first intron"
hESC_tss[hESC_tss[,7]=="other_exon",7] <- "other exon"
hESC_tss[hESC_tss[,7]=="other_intron",7] <- "other intron"
hESC_tss[hESC_tss[,7]=="5UTR",7] <- "5'UTR"
hESC_tss[hESC_tss[,7]=="3UTR",7] <- "3'UTR"
hESC_tss[hESC_tss[,7]=="TSS1000",7] <- "TSS±1000"
hESC_tss_df <- as.data.frame(table(hESC_tss[,7]))
hESC_tss_df$Var1 <- factor(hESC_tss_df$Var1, levels=c("TSS±1000", "5'UTR", "first exon","first intron","other exon","other intron","3'UTR","intergenic"), ordered=TRUE)
hESC_tss_df <- hESC_tss_df[order(hESC_tss_df[,1]),]
hESC_tss_df[,1] <- paste(hESC_tss_df[,1]," (",round(hESC_tss_df[,2]/sum(hESC_tss_df[,2])*100,2),"%)",sep = "")
hESC_tss_df$Var1 <- factor(hESC_tss_df$Var1, levels=c(grep("TSS",hESC_tss_df[,1],value = T),
                                                    grep("5'UTR",hESC_tss_df[,1],value = T),
                                                    grep("first exon",hESC_tss_df[,1],value = T),
                                                    grep("first intron",hESC_tss_df[,1],value = T),
                                                    grep("other exon",hESC_tss_df[,1],value = T),
                                                    grep("other intron",hESC_tss_df[,1],value = T),
                                                    grep("3'UTR",hESC_tss_df[,1],value = T),
                                                    grep("intergenic",hESC_tss_df[,1],value = T)), ordered=TRUE)
hESC_tss_df <- hESC_tss_df[order(hESC_tss_df[,1]),]
for(i in 1:nrow(hESC_tes)) {
  a <- hESC_tes_region[hESC_tes_region[,1]==hESC_tes[i,1] & hESC_tes_region[,3]==hESC_tes[i,3] & hESC_tes_region[,6]==hESC_tes[i,6],]
  b <- ordered(a[,13])
  b <- ordered(b, levels = c("TES1000", "3UTR", "last_exon","one_exon","last_intron","one_intron","Other_exon","first_exon","other_intron","first_intron","5UTR","TSS1000") )
  b <- as.character(b[order(b)][1])
  #print(b)
  #b <- b[1]
  #print(b)
  hESC_tes[i,7] <- b
}
hESC_tes[is.na(hESC_tes[,7]),7] <- "intergenic"
hESC_tes[hESC_tes[,7]=="TSS1000",7] <- "intergenic"
hESC_tes[hESC_tes[,7]=="first_exon",7] <- "other exon"
hESC_tes[hESC_tes[,7]=="first_intron",7] <- "other intron"
hESC_tes[hESC_tes[,7]=="one_exon",7] <- "last exon"
hESC_tes[hESC_tes[,7]=="one_intron",7] <- "last intron"
hESC_tes[hESC_tes[,7]=="last_exon",7] <- "last exon"
hESC_tes[hESC_tes[,7]=="last_intron",7] <- "last intron"
hESC_tes[hESC_tes[,7]=="other_exon",7] <- "other exon"
hESC_tes[hESC_tes[,7]=="other_intron",7] <- "other intron"
hESC_tes[hESC_tes[,7]=="5UTR",7] <- "5'UTR"
hESC_tes[hESC_tes[,7]=="3UTR",7] <- "3'UTR"
hESC_tes[hESC_tes[,7]=="TES1000",7] <- "TES±1000"
hESC_tes_df <- as.data.frame(table(hESC_tes[,7]))
hESC_tes_df$Var1 <- factor(hESC_tes_df$Var1, levels=c("TES±1000", "3'UTR", "last exon","last intron","other exon","other intron","5'UTR","intergenic"), ordered=TRUE)
hESC_tes_df <- hESC_tes_df[order(hESC_tes_df[,1]),]
hESC_tes_df[,1] <- paste(hESC_tes_df[,1]," (",round(hESC_tes_df[,2]/sum(hESC_tes_df[,2])*100,2),"%)",sep = "")
hESC_tes_df$Var1 <- factor(hESC_tes_df$Var1, levels=c(grep("TES",hESC_tes_df[,1],value = T),
                                                    grep("3'UTR",hESC_tes_df[,1],value = T),
                                                    grep("last exon",hESC_tes_df[,1],value = T),
                                                    grep("last intron",hESC_tes_df[,1],value = T),
                                                    grep("other exon",hESC_tes_df[,1],value = T),
                                                    grep("other intron",hESC_tes_df[,1],value = T),
                                                    grep("5'UTR",hESC_tes_df[,1],value = T),
                                                    grep("intergenic",hESC_tes_df[,1],value = T)), ordered=TRUE)
hESC_tes_df <- hESC_tes_df[order(hESC_tes_df[,1]),]


p <- ggplot(hESC_tss_df, aes(x = 1, weight = Freq, fill = Var1)) +
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

ggsave(p, filename = file.path(outdir, "pieTSS.pdf"), height = 4, width = 4)

p <- ggplot(hESC_tes_df, aes(x = 1, weight = Freq, fill = Var1)) +
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

ggsave(p, filename = file.path(outdir, "pieTES.pdf"), height = 4, width = 4)
