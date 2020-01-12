# fig2
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
load("fig2.RData")



### fig2 b; fig2 c

options(stringsAsFactors = FALSE)
options(scipen = 100)



#gencode_tss <- read.table("G:/annotation/gencode_mm10_all_gene_all_transcript_tss.bed",sep = "\t")
#gencode_tes <- read.table("G:/annotation/gencode_mm10_all_gene_all_transcript_tes.bed",sep = "\t")

#D_dominant_tss_in_gene_1 <- cbind(read.csv("20190711DRGtss_peak_new.csv"),
#                                  read.csv("majority_vote_tss_test_DRG.csv")[,c(12,13)])
#D_dominant_tss_in_gene_1 <- D_dominant_tss_in_gene_1[D_dominant_tss_in_gene_1[,17]==1,]
#D_dominant_tss_in_gene_1 <- data.frame(V1=D_dominant_tss_in_gene_1[,2],V2=D_dominant_tss_in_gene_1[,7]-1,V3=D_dominant_tss_in_gene_1[,7],
#                                       V10=D_dominant_tss_in_gene_1[,1],V5=D_dominant_tss_in_gene_1[,6],V6=D_dominant_tss_in_gene_1[,5])


#O_dominant_tss_in_gene_1 <- cbind(read.csv("20190717D4_tss_peak_new.csv"),
#                                  read.csv("majority_vote_tss_test_D4.csv")[,c(12,13)])
#O_dominant_tss_in_gene_1 <- O_dominant_tss_in_gene_1[O_dominant_tss_in_gene_1[,17]==1,]
#O_dominant_tss_in_gene_1 <- data.frame(V1=O_dominant_tss_in_gene_1[,2],V2=O_dominant_tss_in_gene_1[,7]-1,V3=O_dominant_tss_in_gene_1[,7],
#                                       V10=O_dominant_tss_in_gene_1[,1],V5=O_dominant_tss_in_gene_1[,6],V6=O_dominant_tss_in_gene_1[,5])









#D_dominant_tes_in_gene_1 <- cbind(read.csv("20190711DRGtes_peak_new.csv"),
#                                  read.csv("majority_vote_tes_test_DRG.csv")[,c(12,13)])
#D_dominant_tes_in_gene_1 <- D_dominant_tes_in_gene_1[D_dominant_tes_in_gene_1[,17]==1,]
#D_dominant_tes_in_gene_1 <- data.frame(V1=D_dominant_tes_in_gene_1[,2],V2=D_dominant_tes_in_gene_1[,7]-1,V3=D_dominant_tes_in_gene_1[,7],
#                                       V10=D_dominant_tes_in_gene_1[,1],V5=D_dominant_tes_in_gene_1[,6],V6=D_dominant_tes_in_gene_1[,5])


#O_dominant_tes_in_gene_1 <- cbind(read.csv("20190717D4tes_peak_new.csv"),
#                                  read.csv("majority_vote_tes_test_D4.csv")[,c(12,13)])
#O_dominant_tes_in_gene_1 <- O_dominant_tes_in_gene_1[O_dominant_tes_in_gene_1[,17]==1,]
#O_dominant_tes_in_gene_1 <- data.frame(V1=O_dominant_tes_in_gene_1[,2],V2=O_dominant_tes_in_gene_1[,7]-1,V3=O_dominant_tes_in_gene_1[,7],
#                                       V10=O_dominant_tes_in_gene_1[,1],V5=O_dominant_tes_in_gene_1[,6],V6=O_dominant_tes_in_gene_1[,5])









D_dominant_tss_novel <- data.frame()
for(i in 1:nrow(D_dominant_tss_in_gene_1)) {
  a <- gencode_tss[gencode_tss[,1]==D_dominant_tss_in_gene_1[i,1] & gencode_tss[,6]==D_dominant_tss_in_gene_1[i,6] & gencode_tss[,3]<D_dominant_tss_in_gene_1[i,3]+101 & gencode_tss[,3]>D_dominant_tss_in_gene_1[i,3]-101,]
  if(nrow(a)==0) D_dominant_tss_novel <- rbind(D_dominant_tss_novel,D_dominant_tss_in_gene_1[i,])
}



O_dominant_tss_novel <- data.frame()
for(i in 1:nrow(O_dominant_tss_in_gene_1)) {
  a <- gencode_tss[gencode_tss[,1]==O_dominant_tss_in_gene_1[i,1] & gencode_tss[,6]==O_dominant_tss_in_gene_1[i,6] & gencode_tss[,3]<O_dominant_tss_in_gene_1[i,3]+101 & gencode_tss[,3]>O_dominant_tss_in_gene_1[i,3]-101,]
  if(nrow(a)==0) O_dominant_tss_novel <- rbind(O_dominant_tss_novel,O_dominant_tss_in_gene_1[i,])
}








D_dominant_tes_novel <- data.frame()
for(i in 1:nrow(D_dominant_tes_in_gene_1)) {
  a <- gencode_tes[gencode_tes[,1]==D_dominant_tes_in_gene_1[i,1] & gencode_tes[,6]==D_dominant_tes_in_gene_1[i,6] & gencode_tes[,3]<D_dominant_tes_in_gene_1[i,3]+101 & gencode_tes[,3]>D_dominant_tes_in_gene_1[i,3]-101,]
  if(nrow(a)==0) D_dominant_tes_novel <- rbind(D_dominant_tes_novel,D_dominant_tes_in_gene_1[i,])
}



O_dominant_tes_novel <- data.frame()
for(i in 1:nrow(O_dominant_tes_in_gene_1)) {
  a <- gencode_tes[gencode_tes[,1]==O_dominant_tes_in_gene_1[i,1] & gencode_tes[,6]==O_dominant_tes_in_gene_1[i,6] & gencode_tes[,3]<O_dominant_tes_in_gene_1[i,3]+101 & gencode_tes[,3]>O_dominant_tes_in_gene_1[i,3]-101,]
  if(nrow(a)==0) O_dominant_tes_novel <- rbind(O_dominant_tes_novel,O_dominant_tes_in_gene_1[i,])
}










##### novel gene



#normalizeTagCount(myCAGEsetDtss, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
#normalizeTagCount(myCAGEsetoocytetss, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)

#normalizeTagCount(myCAGEsetDtes, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
#normalizeTagCount(myCAGEsetoocytetes, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)



#i=10
#clusterCTSS(object = myCAGEsetDtss, method = "distclu", threshold = i, 
#            nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,
#            keepSingletonsAbove = 1, maxDist = 20,
#            useMulticore = FALSE, nrCores = 8)
#clusterCTSS(object = myCAGEsetoocytetss, method = "distclu", threshold = i, 
#            nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,
#            keepSingletonsAbove = 1, maxDist = 20,
#            useMulticore = FALSE, nrCores = 8)
#clusterCTSS(object = myCAGEsetDtes, method = "distclu", threshold = i, 
#            nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,
#            keepSingletonsAbove = 1, maxDist = 20,
#            useMulticore = FALSE, nrCores = 8)
#clusterCTSS(object = myCAGEsetoocytetes, method = "distclu", threshold = i, 
#            nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,
#            keepSingletonsAbove = 1, maxDist = 20,
#            useMulticore = FALSE, nrCores = 8)


#tc_Dtss <- tagClusters(myCAGEsetDtss)[[1]]
#tc_oocytetss <- tagClusters(myCAGEsetoocytetss)[[1]]
#tc_Dtes <- tagClusters(myCAGEsetDtes)[[1]]
#tc_oocytetes <- tagClusters(myCAGEsetoocytetes)[[1]]



#D_dominant_tss <- data.frame(tc_Dtss[,2],tc_Dtss[,8]-1,tc_Dtss[,8],tc_Dtss[,c(6,6,5)],stringsAsFactors = FALSE)
#oocyte_dominant_tss <- data.frame(tc_oocytetss[,2],tc_oocytetss[,8]-1,tc_oocytetss[,8],tc_oocytetss[,c(6,6,5)],stringsAsFactors = FALSE)
#D_dominant_tes <- data.frame(tc_Dtes[,2],tc_Dtes[,8]-1,tc_Dtes[,8],tc_Dtes[,c(6,6,5)],stringsAsFactors = FALSE)
#oocyte_dominant_tes <- data.frame(tc_oocytetes[,2],tc_oocytetes[,8]-1,tc_oocytetes[,8],tc_oocytetes[,c(6,6,5)],stringsAsFactors = FALSE)


#D_dominant_tss_plus <- D_dominant_tss[D_dominant_tss[,6]=="+",]
#D_dominant_tss_minus <- D_dominant_tss[D_dominant_tss[,6]=="-",]
#D_dominant_tes_plus <- D_dominant_tes[D_dominant_tes[,6]=="+",]
#D_dominant_tes_minus <- D_dominant_tes[D_dominant_tes[,6]=="-",]
#oocyte_dominant_tss_plus <- oocyte_dominant_tss[oocyte_dominant_tss[,6]=="+",]
#oocyte_dominant_tss_minus <- oocyte_dominant_tss[oocyte_dominant_tss[,6]=="-",]
#oocyte_dominant_tes_plus <- oocyte_dominant_tes[oocyte_dominant_tes[,6]=="+",]
#oocyte_dominant_tes_minus <- oocyte_dominant_tes[oocyte_dominant_tes[,6]=="-",]




#write.table(D_dominant_tss_plus,file = "D_dominant_tss_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(D_dominant_tss_minus,file = "D_dominant_tss_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(D_dominant_tes_plus,file = "D_dominant_tes_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(D_dominant_tes_minus,file = "D_dominant_tes_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(oocyte_dominant_tss_plus,file = "oocyte_dominant_tss_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(oocyte_dominant_tss_minus,file = "oocyte_dominant_tss_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(oocyte_dominant_tes_plus,file = "oocyte_dominant_tes_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
#write.table(oocyte_dominant_tes_minus,file = "oocyte_dominant_tes_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")





#gencode_tss <- read.table("G:/annotation/gencode_mm10_all_gene_all_transcript_tss.bed",header = FALSE,sep = "\t")


#gencode_tes <- read.table("G:/annotation/gencode_mm10_all_gene_all_transcript_tes.bed",header = FALSE,sep = "\t")








#D_dominant_tss_novel <- data.frame()
#for(i in 1:nrow(D_dominant_tss)) {
#  a <- gencode_tss[gencode_tss[,1]==D_dominant_tss[i,1] & 
#                     gencode_tss[,6]==D_dominant_tss[i,6] & 
#                     abs(gencode_tss[,3]-D_dominant_tss[i,3])<101,]
#  if(nrow(a)==0) D_dominant_tss_novel <- rbind(D_dominant_tss_novel,D_dominant_tss[i,])
#}


#D_dominant_tes_novel <- data.frame()
#for(i in 1:nrow(D_dominant_tes)) {
#  a <- gencode_tes[gencode_tes[,1]==D_dominant_tes[i,1] & 
#                     gencode_tes[,6]==D_dominant_tes[i,6] & 
#                     abs(gencode_tes[,3]-D_dominant_tes[i,3])<101,]
#  if(nrow(a)==0) D_dominant_tes_novel <- rbind(D_dominant_tes_novel,D_dominant_tes[i,])
#}


#oocyte_dominant_tss_novel <- data.frame()
#for(i in 1:nrow(oocyte_dominant_tss)) {
#  a <- gencode_tss[gencode_tss[,1]==oocyte_dominant_tss[i,1] & 
#                     gencode_tss[,6]==oocyte_dominant_tss[i,6] & 
#                     abs(gencode_tss[,3]-oocyte_dominant_tss[i,3])<101,]
#  if(nrow(a)==0) oocyte_dominant_tss_novel <- rbind(oocyte_dominant_tss_novel,oocyte_dominant_tss[i,])
#}


#oocyte_dominant_tes_novel <- data.frame()
#for(i in 1:nrow(oocyte_dominant_tes)) {
#  a <- gencode_tes[gencode_tes[,1]==oocyte_dominant_tes[i,1] & 
#                     gencode_tes[,6]==oocyte_dominant_tes[i,6] & 
#                     abs(gencode_tes[,3]-oocyte_dominant_tes[i,3])<101,]
#  if(nrow(a)==0) oocyte_dominant_tes_novel <- rbind(oocyte_dominant_tes_novel,oocyte_dominant_tes[i,])
#}









##### bedtools intersect

#bedtools intersect -a D_dominant_tss_plus.bed -b gencode_mm10_all_gene_intergenic_plus.bed -wa -wb >D_dominant_tss_plus_intergenic.bed
#bedtools intersect -a oocyte_dominant_tss_plus.bed -b gencode_mm10_all_gene_intergenic_plus.bed -wa -wb >oocyte_dominant_tss_plus_intergenic.bed
#bedtools intersect -a D_dominant_tes_plus.bed -b gencode_mm10_all_gene_intergenic_plus.bed -wa -wb >D_dominant_tes_plus_intergenic.bed
#bedtools intersect -a oocyte_dominant_tes_plus.bed -b gencode_mm10_all_gene_intergenic_plus.bed -wa -wb >oocyte_dominant_tes_plus_intergenic.bed
#bedtools intersect -a D_dominant_tss_minus.bed -b gencode_mm10_all_gene_intergenic_minus.bed -wa -wb >D_dominant_tss_minus_intergenic.bed
#bedtools intersect -a oocyte_dominant_tss_minus.bed -b gencode_mm10_all_gene_intergenic_minus.bed -wa -wb >oocyte_dominant_tss_minus_intergenic.bed
#bedtools intersect -a D_dominant_tes_minus.bed -b gencode_mm10_all_gene_intergenic_minus.bed -wa -wb >D_dominant_tes_minus_intergenic.bed
#bedtools intersect -a oocyte_dominant_tes_minus.bed -b gencode_mm10_all_gene_intergenic_minus.bed -wa -wb >oocyte_dominant_tes_minus_intergenic.bed






#D_tss_intergenic_plus <- read.table("D_dominant_tss_plus_intergenic.bed",header = FALSE)
#oocyte_tss_intergenic_plus <- read.table("oocyte_dominant_tss_plus_intergenic.bed",header = FALSE)
#D_tes_intergenic_plus <- read.table("D_dominant_tes_plus_intergenic.bed",header = FALSE)
#oocyte_tes_intergenic_plus <- read.table("oocyte_dominant_tes_plus_intergenic.bed",header = FALSE)
#D_tss_intergenic_minus <- read.table("D_dominant_tss_minus_intergenic.bed",header = FALSE)
#oocyte_tss_intergenic_minus <- read.table("oocyte_dominant_tss_minus_intergenic.bed",header = FALSE)
#D_tes_intergenic_minus <- read.table("D_dominant_tes_minus_intergenic.bed",header = FALSE)
#oocyte_tes_intergenic_minus <- read.table("oocyte_dominant_tes_minus_intergenic.bed",header = FALSE)



## novel gene here for easy

D_tss_intergenic_plus <- D_tss_intergenic_plus[,c(1,2,3,10,5,6)]
oocyte_tss_intergenic_plus <- oocyte_tss_intergenic_plus[,c(1,2,3,10,5,6)]
D_tes_intergenic_plus <- D_tes_intergenic_plus[,c(1,2,3,10,5,6)]
oocyte_tes_intergenic_plus <- oocyte_tes_intergenic_plus[,c(1,2,3,10,5,6)]
D_tss_intergenic_minus <- D_tss_intergenic_minus[,c(1,2,3,10,5,6)]
oocyte_tss_intergenic_minus <- oocyte_tss_intergenic_minus[,c(1,2,3,10,5,6)]
D_tes_intergenic_minus <- D_tes_intergenic_minus[,c(1,2,3,10,5,6)]
oocyte_tes_intergenic_minus <- oocyte_tes_intergenic_minus[,c(1,2,3,10,5,6)]




D_plus_novel_gene <- merge(D_tss_intergenic_plus,D_tes_intergenic_plus,by = "V10")
oocyte_plus_novel_gene <- merge(oocyte_tss_intergenic_plus,oocyte_tes_intergenic_plus,by = "V10")
D_minus_novel_gene <- merge(D_tss_intergenic_minus,D_tes_intergenic_minus,by = "V10")
oocyte_minus_novel_gene <- merge(oocyte_tss_intergenic_minus,oocyte_tes_intergenic_minus,by = "V10")

D_novel_gene <- rbind(D_plus_novel_gene,D_minus_novel_gene)
oocyte_novel_gene <- rbind(oocyte_plus_novel_gene,oocyte_minus_novel_gene)






D_plus_novel_gene <- D_plus_novel_gene[D_plus_novel_gene[,4]<D_plus_novel_gene[,9],]
oocyte_plus_novel_gene <- oocyte_plus_novel_gene[oocyte_plus_novel_gene[,4]<oocyte_plus_novel_gene[,9],]
D_minus_novel_gene <- D_minus_novel_gene[D_minus_novel_gene[,4]>D_minus_novel_gene[,9],]
oocyte_minus_novel_gene <- oocyte_minus_novel_gene[oocyte_minus_novel_gene[,4]>oocyte_minus_novel_gene[,9],]



D_plus_novel_gene <- D_plus_novel_gene[abs(D_plus_novel_gene[,4]-D_plus_novel_gene[,9])<100000,]
oocyte_plus_novel_gene <- oocyte_plus_novel_gene[abs(oocyte_plus_novel_gene[,4]-oocyte_plus_novel_gene[,9])<100000,]
D_minus_novel_gene <- D_minus_novel_gene[abs(D_minus_novel_gene[,4]-D_minus_novel_gene[,9])<100000,]
oocyte_minus_novel_gene <- oocyte_minus_novel_gene[abs(oocyte_minus_novel_gene[,4]-oocyte_minus_novel_gene[,9])<100000,]























































a <- data.frame(V1=c(2031,4693,71),V2=c("Novel TSS","Novel TES","Novel gene"))
a <- data.frame(V1=c(3102,5746,125),V2=c("Novel TSS","Novel TES","Novel gene"))


a[,2] <- factor(a[,2],levels = c("Novel TSS","Novel TES","Novel gene"),ordered = T)
a <- a[order(a[,2]),]
ggplot(a,aes(x=V2,y=V1,fill=V2)) +
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
  #scale_y_continuous(breaks = seq(0,3000,1000),limits = c(0,3500),expand=c(0,0)) 
  scale_y_continuous(breaks = seq(0,5000,1000),limits = c(0,5000),expand=c(0,0),labels = c("0","1,000","2,000","3,000","4,000","5,000"))

