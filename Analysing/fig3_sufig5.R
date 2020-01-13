# fig3 sufig5
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
load("fig3_sufig5.RData")

### fig3 a; sufig5 a




#make sure that no *ctss files are in working directory


####5' single cell callpeak
#for(i in objects()[grep("5cap",objects())]) {
#  a <- paste('normalizeTagCount(',i,', method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)',sep = "\t")
#  eval(parse(text=a))
#  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 10, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
#  eval(parse(text=a))
#}


####3' single cell callpeak
#for(i in objects()[grep("3tail",objects())]) {
#  a <- paste('normalizeTagCount(',i,', method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)',sep = "\t")
#  eval(parse(text=a))
#  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 10, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
#  eval(parse(text=a))
#}

#rm(i)
#rm(a)
#objectlist <- objects()



#for(i in grep("myCAGEset",objects(),value = T)) {
#  a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- tagClusters(',i,')[[1]]',sep = "")
#  eval(parse(text=a))
#}

#for(i in grep("5cap|tss",grep("tc",objects(),value = T),value = T)) {
#  a <- paste('colnames(',i,') <- c("cluster","chr","start","end","strand","rpm","nr_tss","dominant_tss","rpm.dominant_tss")',sep = "")
#  eval(parse(text=a))
#}

#for(i in grep("3tail|tes",grep("tc",objects(),value = T),value = T)) {
#  a <- paste('colnames(',i,') <- c("cluster","chr","start","end","strand","rpm","nr_tes","dominant_tes","rpm.dominant_tes")',sep = "")
#  eval(parse(text=a))
#}


#for(i in grep("tc",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.txt',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = T,sep = "\t")',sep = "")
#  eval(parse(text=a))
#}

#for(i in grep("tc",objects(),value = TRUE)) {
#  a <- paste(paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tss',sep = ""),' <- data.frame(',i,'[,2],',i,'[,8]-1,',i,'[,c(8,6,6,5)],stringsAsFactors = FALSE)',sep = "")
#  eval(parse(text=a))
#}


#for(i in grep("dominant",objects(),value = TRUE)) {
#  a <- paste(paste(i,'_plus',sep = ""),' <- ',i,'[',i,'[,6]=="+",]',sep = "")
#  eval(parse(text=a))
#  a <- paste(paste(i,'_minus',sep = ""),' <- ',i,'[',i,'[,6]=="-",]',sep = "")
#  eval(parse(text=a))
#}

##output dominant plus
#for(i in grep("plus",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
#  eval(parse(text=a))
#}

##output dominant minus
#for(i in grep("minus",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
#  eval(parse(text=a))
#}












####bedtools intersect

#rm(list = ls())


#for(i in grep("bed",grep("stream",list.files(),value = T),value = T)) {
#  a <- paste(strsplit(i,split = "[.]")[[1]][1],' <- read.table("',i,'",header = FALSE)',sep = "")
#  eval(parse(text=a))
#  print(a)
#}

#load("correlation_heatmap.RData")
#load("correlation_heatmap1.RData")
#load("correlation_heatmap2.RData")

for(i in grep("stream",objects(),value = T)) {
  a <- paste(i,' <- ',i,'[,c(1,3,6,4)]',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste(i,' <- unique(',i,')',sep = "")
  eval(parse(text=a))
  print(a)
}






###combine plus and minus
for(i in grep("plus",grep("stream",objects(),value = T),value = T)) {
  a <- paste(strsplit(i,split = "plus")[[1]][1],strsplit(i,split = "plus_")[[1]][2],' <- rbind(',i,',',strsplit(i,split = "plus")[[1]][1],'minus',strsplit(i,split = "plus")[[1]][2],')',sep = "")
  eval(parse(text=a))
  print(a)
}



for(i in grep("cap",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste('b <- CTSStagCount(',i,')',sep = "")
  eval(parse(text=a))
  print(a)
  a <- 'b <- sum(b[,4])'
  eval(parse(text=a))
  print(a)
  a <- paste(strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_upstream2k_and_genebody[,4] <- round(',strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_upstream2k_and_genebody[,4]/1000000*',b,')',sep = "")
  eval(parse(text=a))
  print(a)
}


for(i in grep("tail",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste('b <- CTSStagCount(',i,')',sep = "")
  eval(parse(text=a))
  print(a)
  a <- 'b <- sum(b[,4])'
  eval(parse(text=a))
  print(a)
  a <- paste(strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_gene_genebody_and_downstream20k[,4] <- round(',strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_gene_genebody_and_downstream20k[,4]/1000000*',b,')',sep = "")
  eval(parse(text=a))
  print(a)
}






##output
for(i in grep("tss_gene",objects(),value = T)) {
  a <- paste('write.table(',i,',"',paste(i,'.ctss',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  eval(parse(text=a))
}

for(i in grep("tss_upstream",objects(),value = T)) {
  a <- paste('write.table(',i,',"',paste(i,'.ctss',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  eval(parse(text=a))
}





tss_samples <- c()
tss_inputFiles <- c()
for(i in grep("5cap",grep("ctss",list.files(),value = T),value = T)) {
  a <- paste('tss_samples <- c(tss_samples,"',strsplit(i,split = "_dominant")[[1]][1],'")',sep = "")
  eval(parse(text=a))
  a <- paste('tss_inputFiles <- c(tss_inputFiles,"',i,'")',sep = "")
  eval(parse(text=a))
}


tes_samples <- c()
tes_inputFiles <- c()
for(i in grep("3tail",grep("ctss",list.files(),value = T),value = T)) {
  a <- paste('tes_samples <- c(tes_samples,"',strsplit(i,split = "_dominant")[[1]][1],'")',sep = "")
  eval(parse(text=a))
  a <- paste('tes_inputFiles <- c(tes_inputFiles,"',i,'")',sep = "")
  eval(parse(text=a))
}

myCAGEsettss <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = tss_inputFiles, inputFilesType = "ctss",sampleLabels = tss_samples)
myCAGEsettes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = tes_inputFiles, inputFilesType = "ctss",sampleLabels = tes_samples)

getCTSS(myCAGEsettss,removeFirstG = FALSE,correctSystematicG = FALSE)
getCTSS(myCAGEsettes,removeFirstG = FALSE,correctSystematicG = FALSE)

tss_df <- CTSStagCount(myCAGEsettss)
tes_df <- CTSStagCount(myCAGEsettes)


cor_tss_df <- cor(tss_df[,c(-1,-2,-3)])
cor_tes_df <- cor(tes_df[,c(-1,-2,-3)])


for(i in 1:nrow(cor_tss_df)) {
  for(j in 1:ncol(cor_tss_df)) {
    rowname <- rownames(cor_tss_df)[i]
    colname <- colnames(cor_tss_df)[j]
    a <- paste('temp <- data.frame(V1=tss_df$',rowname,',V2=tss_df$',colname,')',sep = "")
    eval(parse(text=a))
    print(a)
    temp <- temp[!(temp[,1]==0 & temp[,2]==0),]
    cor_tss_df[i,j] <- cor(temp[,1],temp[,2],method = "pearson")
  }
}


for(i in 1:nrow(cor_tes_df)) {
  for(j in 1:ncol(cor_tes_df)) {
    rowname <- rownames(cor_tes_df)[i]
    colname <- colnames(cor_tes_df)[j]
    a <- paste('temp <- data.frame(V1=tes_df$',rowname,',V2=tes_df$',colname,')',sep = "")
    eval(parse(text=a))
    print(a)
    temp <- temp[!(temp[,1]==0 & temp[,2]==0),]
    cor_tes_df[i,j] <- cor(temp[,1],temp[,2],method = "pearson")
  }
}


#pdf(file=paste("heat_cor_5cap.pdf", sep=""), height = 10, width = 10)

heatmap.2(cor_tss_df, 
          Rowv=TRUE, Colv="Rowv", dendrogram='row',
          trace='none',
          scale = "none", 
          col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),
          symbreaks = FALSE,
          breaks = seq(0,1,0.004),
          srtCol=45, adjCol=c(1,1),
          #key.title = "Pearson",
          keysize = 0.8,
          RowSideColors = c(rep("#0000CD",25),rep("#A52A2A",8),rep("#1C1C1C",22)),
          ColSideColors = c(rep("#0000CD",25),rep("#A52A2A",8),rep("#1C1C1C",22)),
          labCol = "",
          labRow = "",
          density.info="none") 
#legend("topright",      # location of the legend on the heatmap plot
#       legend = c("DRG", "Oocyte D3", "Oocyte D4"), # category labels
#       fill = c("#0000CD", "#1C1C1C", "#A52A2A"),  # color key
#       border = c("#0000CD", "#1C1C1C", "#A52A2A"))
dev.off() 

#pdf(file=paste("heat_cor_3tail.pdf", sep=""), height = 10, width = 10)
heatmap.2(cor_tes_df, 
          Rowv=TRUE, Colv="Rowv", dendrogram='row',
          trace='none',
          scale = "none",
          col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),
          symbreaks = FALSE,
          breaks = seq(0,1,0.004),
          srtCol=45, adjCol=c(1,1),
          #key.title = "Pearson",
          keysize = 0.8,
          RowSideColors = c(rep("#0000CD",25),rep("#A52A2A",8),rep("#1C1C1C",22)),
          ColSideColors = c(rep("#0000CD",25),rep("#A52A2A",8),rep("#1C1C1C",22)),
          labCol = "",
          labRow = "",
          density.info="none") 

dev.off() 




### fig3 b; sufig5 b; sufig5 c; sufig5 d; sufig5 e; sufig5 f













condition <- factor(c(rep("D",19),rep("P",12)))

dds <- DESeqDataSetFromMatrix(smart_seq2_readcount, DataFrame(condition), design= ~ condition )

dds <- DESeq(dds)

res <- results(dds)

mcols(res)

summary(res)

a <- as.data.frame(res)

no_sign <- subset(res, padj >= 0.05)

no_sign <- as.data.frame(no_sign)





boxplot <- melt(smart_seq2_RPM[rownames(smart_seq2_RPM) %in% "Grpel1",])
boxplot[,1] <- as.character(boxplot[,1])
boxplot[1:19,1] <- "DRG"
boxplot[20:31,1] <- "oocyte D3"
boxplot[,2] <- log10(boxplot[,2]+1)

my_comparison <- list(c("DRG","oocyte D3"))
ggboxplot(boxplot, x="variable", y="value", fill = "variable", palette = "jco",xlab="",legend = "right") + 
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  labs(y="log10 (RPM+1)")+
  #labs(y="log10 (absolute gene expression +1)")+
  stat_compare_means(comparisons = my_comparison,aes(label=p.sign..),label = "p-value", method = "wilcox.test")









hm_df_5cap <- log10(hm_df_5cap+1)




################
rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))

hm_df_5cap <- hm_df_5cap[rowname %in% rownames(no_sign),]


################






###############

double <- as.data.frame(table(rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))))
double <- as.character(double[double[,2]==2,1])
hm_df_5cap <- hm_df_5cap[unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1])) %in% double,]

##############







boxplot3 <- melt(hm_df_5cap[rownames(hm_df_5cap) %in% "Grpel1 _D",])
boxplot4 <- melt(hm_df_5cap[rownames(hm_df_5cap) %in% "Grpel1 _P",])

boxplot <- data.frame(variable=rownames(boxplot),value=boxplot[,1])
boxplot3[,1] <- as.character(boxplot3[,1])
boxplot3[1:25,1] <- "DRG"
boxplot3[26:47,1] <- "oocyte D3"

boxplot4[,1] <- as.character(boxplot4[,1])
boxplot4[1:25,1] <- "DRG"
boxplot4[26:47,1] <- "oocyte D3"


my_comparison <- list(c("DRG","oocyte D3"))
ggboxplot(boxplot3, x="variable", y="value", fill = "variable", palette = "jco",xlab="",legend = "right") + 
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  labs(y="log10 (RPM+1)")+
  stat_compare_means(comparisons = my_comparison,aes(label=p.sign..),label = "p-value", method = "wilcox.test")

ggboxplot(boxplot4, x="variable", y="value", fill = "variable", palette = "jco",xlab="",legend = "right") + 
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  labs(y="log10 (RPM+1)")+
  stat_compare_means(comparisons = my_comparison,aes(label=p.sign..),label = "p-value", method = "wilcox.test")







a <- pheatmap(hm_df_5cap,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
a <- pheatmap(hm_df_5cap)

summary(a)



order_row= c(seq(1,nrow(hm_df_5cap),2),seq(2,nrow(hm_df_5cap),2))
hm_df_5cap <- hm_df_5cap[order_row,]

annotation_col = data.frame(celltype=factor(rep(c("Oocyte_D3", "DRG"), c(22,25))))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]



temp <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = "_"), function(x) x[2]))  ###根据isoform分开

temp[temp=="D"] <- "isoform_1"  ###

temp[temp=="P"] <- "isoform_2"  ###




annotation_row = data.frame(isoform=factor(temp))

ann_colors <- list(celltype = c(Oocyte_D3 = "#1C1C1C", DRG = "#0000CD"),
                   isoform = c(isoform_1 = "#1B9E77", isoform_2 = "#D95F02"))



rownames(annotation_row) <- rownames(hm_df_5cap)   ###



pheatmap(hm_df_5cap,
         scale="row",
         col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(200)),
         cluster_rows = F,
         fontsize_row = 5,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_rownames = F,
         show_colnames = F,
         border_color = F)



hm_df_5cap_normal <- as.matrix(hm_df_5cap_normal)
hm_df_5cap_normal <- log10(hm_df_5cap_normal+1)


##############

hm_df_5cap_normal <- hm_df_5cap_normal[rownames(hm_df_5cap_normal) %in% as.character(as.data.frame(table(unlist(strsplit(rownames(hm_df_5cap),split = " _"))))[,1]),]

##############

a <- pheatmap(hm_df_5cap_normal,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
a <- pheatmap(hm_df_5cap_normal)
summary(a)
order_row = a$tree_row$order
hm_df_5cap_normal <- hm_df_5cap_normal[order_row,]


annotation_col = data.frame(
  celltype=factor(rep(c("DRG","Oocyte_D3", "DRG","Oocyte_D3"), c(13,1,12,21))))

ann_colors <- list(celltype = c(Oocyte_D3 = "#1C1C1C", DRG = "#0000CD"))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]


pheatmap(hm_df_5cap_normal,
         scale="row",
         #color = colorRampPalette(colors = c("blue","white","red"))(100),
         col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),
         cluster_rows = F,
         fontsize_row = 5,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         show_rownames = F,
         show_colnames = F,
         border_color = F)















#3'






hm_df_3tail <- log10(hm_df_3tail+1)




################
rowname <- unlist(lapply(strsplit(rownames(hm_df_3tail),split = " _"), function(x) x[1]))

hm_df_3tail <- hm_df_3tail[rowname %in% rownames(no_sign),]


################

##############

double <- as.data.frame(table(rowname <- unlist(lapply(strsplit(rownames(hm_df_3tail),split = " _"), function(x) x[1]))))
double <- as.character(double[double[,2]==2,1])

hm_df_3tail <- hm_df_3tail[unlist(lapply(strsplit(rownames(hm_df_3tail),split = " _"), function(x) x[1])) %in% double,]

##############





a <- pheatmap(hm_df_3tail,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)

summary(a)



order_row= c(seq(1,nrow(hm_df_3tail),2),seq(2,nrow(hm_df_3tail),2))
hm_df_3tail <- hm_df_3tail[order_row,]



annotation_col = data.frame(celltype=factor(rep(c("DRG", "Oocyte_D3"), c(25,22))))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]


temp <- unlist(lapply(strsplit(rownames(hm_df_3tail),split = "_"), function(x) x[2]))

temp[temp=="D"] <- "isoform_1"

temp[temp=="P"] <- "isoform_2"



annotation_row = data.frame(isoform=factor(temp))



rownames(annotation_row) <- rownames(hm_df_3tail)


ann_colors <- list(celltype = c(Oocyte_D3 = "#1C1C1C", DRG = "#0000CD"),
                   isoform = c(isoform_1 = "#1B9E77", isoform_2 = "#D95F02"))



pheatmap(hm_df_3tail,
         scale="row",
         col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),
         cluster_rows = F,
         fontsize_row = 5,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_rownames = F,
         show_colnames = F,
         border_color = F)




##############

hm_df_3tail_normal <- hm_df_3tail_normal[hm_df_3tail_normal[,1] %in% as.character(as.data.frame(table(unlist(strsplit(rownames(hm_df_3tail),split = " _"))))[,1]),]

##############

rownames(hm_df_3tail_normal) <- hm_df_3tail_normal[,1]
hm_df_3tail_normal <- hm_df_3tail_normal[,-1]
hm_df_3tail_normal <- as.matrix(hm_df_3tail_normal)
hm_df_3tail_normal <- log10(hm_df_3tail_normal+1)

a <- pheatmap(hm_df_3tail_normal,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
summary(a)
order_row = a$tree_row$order
hm_df_3tail_normal <- hm_df_3tail_normal[order_row,]


annotation_col = data.frame(celltype=factor(rep(c("Oocyte_D3","DRG"), c(22,25))))

ann_colors <- list(celltype = c(Oocyte_D3 = "#1C1C1C", DRG = "#0000CD"))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]


pheatmap(hm_df_3tail_normal,
         scale="row",
         col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),
         cluster_rows = F,
         fontsize_row = 5,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         show_rownames = F,
         show_colnames = F,
         border_color = F)



### fig3 e; fig3 f

O_merge_1 <- merge(O_dominant_tss_in_gene_1,O_dominant_tes_in_gene_1,by="V10")
P_merge_1 <- merge(P_dominant_tss_in_gene_1,P_dominant_tes_in_gene_1,by="V10")



colnames(O_merge_1) <- c("gene_ID","chr","V3","tss_cor","tss_tpm","strand","V7","V8","tes_cor","tes_tpm","V11")
colnames(P_merge_1) <- c("gene_ID","chr","V3","tss_cor","tss_tpm","strand","V7","V8","tes_cor","tes_tpm","V11")



O_merge_1 <- O_merge_1[,c(1,2,4,5,9,10,6)]
P_merge_1 <- P_merge_1[,c(1,2,4,5,9,10,6)]






for(i in row(O_merge_1)) {
  if(O_merge_1[i,7]=="+") {if(O_merge_1[i,3]>O_merge_1[i,5]) O_merge_1[i,7] <- "."}
  else {if(O_merge_1[i,3]<O_merge_1[i,5]) O_merge_1[i,7] <- "."}
}
for(i in row(P_merge_1)) {
  if(P_merge_1[i,7]=="+") {if(P_merge_1[i,3]>P_merge_1[i,5]) P_merge_1[i,7] <- "."}
  else {if(P_merge_1[i,3]<P_merge_1[i,5]) P_merge_1[i,7] <- "."}
}



O_merge_1 <- O_merge_1[!O_merge_1[,7]==".",]
P_merge_1 <- P_merge_1[!P_merge_1[,7]==".",]










O_major_isoform <- data.frame()
for(i in unique(O_merge_1[,1])) {
  a <- O_merge_1[O_merge_1[,1]==i,]
  a <- a[order(a[,4],a[,6],decreasing = T),]
  O_major_isoform <- rbind(O_major_isoform,a[1,])
  
}


P_major_isoform <- data.frame()
for(i in unique(P_merge_1[,1])) {
  a <- P_merge_1[P_merge_1[,1]==i,]
  a <- a[order(a[,4],a[,6],decreasing = T),]
  P_major_isoform <- rbind(P_major_isoform,a[1,])
  
}





O_major_isoform[,8] <- rep("O",nrow(O_major_isoform))
P_major_isoform[,8] <- rep("P",nrow(P_major_isoform))




union <- rbind(O_major_isoform,P_major_isoform)


union_2_diff <- data.frame()
union_2_same <- data.frame()
union_1 <- data.frame()
for(i in unique(union[,1])) {
  a <- union[union[,1] %in% i,]
  if(nrow(a)==2) {
    if(abs(a[1,3]-a[2,3])<21 & abs(a[1,5]-a[2,5])<21) union_2_same <- rbind(union_2_same,a)
    else union_2_diff <- rbind(union_2_diff,a)
  }
  if(nrow(a)==1) union_1 <- rbind(union_1,a)
}




union_1_true <- data.frame()
for(i in 1:nrow(union_1)) {
  one <- union_1[i,]
  if(one[1,8]=="O") {
    two_tss <- P_dominant_tss_in_gene_1[P_dominant_tss_in_gene_1[,4]==one[1,1],]
    two_tes <- P_dominant_tes_in_gene_1[P_dominant_tes_in_gene_1[,4]==one[1,1],]
    two_tss <- two_tss[order(two_tss[,5],decreasing = TRUE),]
    two_tes <- two_tes[order(two_tes[,5],decreasing = TRUE),]
    b <- two_tss[abs(two_tss[,3]-one[1,3])<21,]
    d <- two_tes[abs(two_tes[,3]-one[1,5])<21,]
    two <- data.frame(gene_ID=one[1,1],chr=one[1,2],tss_cor=b[1,3],tss_tpm=b[1,5],
                      tes_cor=d[1,3],tes_tpm=d[1,5],strand=one[1,7],V8="P")
    merge <- merge(two_tss,two_tes,by="V10")
    colnames(merge) <- c("gene_ID","chr","V3","tss_cor","tss_tpm","strand","V7","V8","tes_cor","tes_tpm","V11")
    merge <- merge[,c(1,2,4,5,9,10,6)]
    for(i in row(merge)) {
      if(merge[i,7]=="+") {if(merge[i,3]>merge[i,5]) merge[i,7] <- "."}
      else {if(merge[i,3]<merge[i,5]) merge[i,7] <- "."}
    }
    merge <- merge[!merge[,7]==".",]
    merge <- merge[order(merge[,4]+merge[,6],decreasing = T),]
    three <- data.frame(merge[1,],V8="P")
    four_tss <- O_dominant_tss_in_gene_1[O_dominant_tss_in_gene_1[,4]==one[1,1],]
    four_tes <- O_dominant_tes_in_gene_1[O_dominant_tes_in_gene_1[,4]==one[1,1],]
    b <- four_tss[abs(four_tss[,3]-three[1,3])<21,]
    d <- four_tes[abs(four_tes[,3]-three[1,5])<21,]
    four <- data.frame(gene_ID=one[1,1],chr=one[1,2],tss_cor=b[1,3],tss_tpm=b[1,5],
                       tes_cor=d[1,3],tes_tpm=d[1,5],strand=one[1,7],V8="O")
  }
  else {
    two_tss <- O_dominant_tss_in_gene_1[O_dominant_tss_in_gene_1[,4]==one[1,1],]
    two_tes <- O_dominant_tes_in_gene_1[O_dominant_tes_in_gene_1[,4]==one[1,1],]
    two_tss <- two_tss[order(two_tss[,5],decreasing = TRUE),]
    two_tes <- two_tes[order(two_tes[,5],decreasing = TRUE),]
    b <- two_tss[abs(two_tss[,3]-one[1,3])<21,]
    d <- two_tes[abs(two_tes[,3]-one[1,5])<21,]
    two <- data.frame(gene_ID=one[1,1],chr=one[1,2],tss_cor=b[1,3],tss_tpm=b[1,5],
                      tes_cor=d[1,3],tes_tpm=d[1,5],strand=one[1,7],V8="O")
    merge <- merge(two_tss,two_tes,by="V10")
    colnames(merge) <- c("gene_ID","chr","V3","tss_cor","tss_tpm","strand","V7","V8","tes_cor","tes_tpm","V11")
    merge <- merge[,c(1,2,4,5,9,10,6)]
    for(i in row(merge)) {
      if(merge[i,7]=="+") {if(merge[i,3]>merge[i,5]) merge[i,7] <- "."}
      else {if(merge[i,3]<merge[i,5]) merge[i,7] <- "."}
    }
    merge <- merge[!merge[,7]==".",]
    merge <- merge[order(merge[,4]+merge[,6],decreasing = T),]
    three <- data.frame(merge[1,],V8="O")    
    four_tss <- P_dominant_tss_in_gene_1[P_dominant_tss_in_gene_1[,4]==one[1,1],]
    four_tes <- P_dominant_tes_in_gene_1[P_dominant_tes_in_gene_1[,4]==one[1,1],]
    b <- four_tss[abs(four_tss[,3]-three[1,3])<21,]
    d <- four_tes[abs(four_tes[,3]-three[1,5])<21,]
    four <- data.frame(gene_ID=one[1,1],chr=one[1,2],tss_cor=b[1,3],tss_tpm=b[1,5],
                       tes_cor=d[1,3],tes_tpm=d[1,5],strand=one[1,7],V8="P")
  }
  union_1_true <- rbind(union_1_true,one,two,three,four)
}











union_2_diff_true <- data.frame()
for(i in 1:(nrow(union_2_diff)/2)) {
  one <- union_2_diff[2*i-1,]
  three <- union_2_diff[2*i,]
  two_tss <- P_dominant_tss_in_gene_1[P_dominant_tss_in_gene_1[,4]==one[1,1],]
  two_tes <- P_dominant_tes_in_gene_1[P_dominant_tes_in_gene_1[,4]==one[1,1],]
  two_tss <- two_tss[order(two_tss[,5],decreasing = TRUE),]
  two_tes <- two_tes[order(two_tes[,5],decreasing = TRUE),]
  b <- two_tss[abs(two_tss[,3]-one[1,3])<21,]
  d <- two_tes[abs(two_tes[,3]-one[1,5])<21,]
  two <- data.frame(gene_ID=one[1,1],chr=one[1,2],tss_cor=b[1,3],tss_tpm=b[1,5],
                    tes_cor=d[1,3],tes_tpm=d[1,5],strand=one[1,7],V8="P")
  four_tss <- O_dominant_tss_in_gene_1[O_dominant_tss_in_gene_1[,4]==three[1,1],]
  four_tes <- O_dominant_tes_in_gene_1[O_dominant_tes_in_gene_1[,4]==three[1,1],]
  four_tss <- four_tss[order(four_tss[,5],decreasing = TRUE),]
  four_tes <- four_tes[order(four_tes[,5],decreasing = TRUE),]
  b <- four_tss[abs(four_tss[,3]-three[1,3])<21,]
  d <- four_tes[abs(four_tes[,3]-three[1,5])<21,]
  four <- data.frame(gene_ID=three[1,1],chr=three[1,2],tss_cor=b[1,3],tss_tpm=b[1,5],
                     tes_cor=d[1,3],tes_tpm=d[1,5],strand=three[1,7],V8="O")
  
  union_2_diff_true <- rbind(union_2_diff_true,one,two,three,four)
}







for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(is.na(union_2_diff_true[4*i-2,3])) {
    union_2_diff_true[4*i-2,3] <- union_2_diff_true[4*i-3,3]
    union_2_diff_true[4*i-2,4] <- 0.01
  }
  if(is.na(union_2_diff_true[4*i-2,5])) {
    union_2_diff_true[4*i-2,5] <- union_2_diff_true[4*i-3,5]
    union_2_diff_true[4*i-2,6] <- 0.5
  }
  if(is.na(union_2_diff_true[4*i,3])) {
    union_2_diff_true[4*i,3] <- union_2_diff_true[4*i-1,3]
    union_2_diff_true[4*i,4] <- 0.01
  }
  if(is.na(union_2_diff_true[4*i,5])) {
    union_2_diff_true[4*i,5] <- union_2_diff_true[4*i-1,5]
    union_2_diff_true[4*i,6] <- 0.5
  }
}






for(i in 1:(nrow(union_2_diff_true)/4)) {
  a <- abs(union_2_diff_true[4*i-3,3]-union_2_diff_true[4*i-3,5])
  b <- abs(union_2_diff_true[4*i-1,3]-union_2_diff_true[4*i-1,5])
  print(a)
  print(b)
  union_2_diff_true[4*i-3,9] <- a
  union_2_diff_true[4*i-1,9] <- b
  if(a>b) {
    union_2_diff_true[4*i-3,10] <- "long"
    union_2_diff_true[4*i-2,10] <- "long"
    union_2_diff_true[4*i-1,10] <- "short"
    union_2_diff_true[4*i,10] <- "short"}
  if(a<b) {
    union_2_diff_true[4*i-3,10] <- "short"
    union_2_diff_true[4*i-2,10] <- "short"
    union_2_diff_true[4*i-1,10] <- "long"
    union_2_diff_true[4*i,10] <- "long"}
}



union_2_diff_true_same_tss <- data.frame()
union_2_diff_true_same_tes <- data.frame()
union_2_diff_true_total_diff <- data.frame()
for(i in 1:(nrow(union_2_diff_true)/4)) {
  a <- union_2_diff_true[(4*i-3):(4*i),]
  if((abs(a[1,3]-a[3,3]))<21) union_2_diff_true_same_tss <- rbind(union_2_diff_true_same_tss,a)
  if((abs(a[1,5]-a[3,5]))<21) union_2_diff_true_same_tes <- rbind(union_2_diff_true_same_tes,a)
  if((abs(a[1,5]-a[3,5]))>20 & ((abs(a[1,3]-a[3,3]))>20)) union_2_diff_true_total_diff <- rbind(union_2_diff_true_total_diff,a)
}


union_2_diff_0515true <- rbind(data.frame(union_2_diff_true_same_tss,type="same_tss"),
                               data.frame(union_2_diff_true_same_tes,type="same_tes"),
                               data.frame(union_2_diff_true_total_diff,type="total_diff"))
union_2_diff_0515true[,12] <- "same"
union_2_diff_0515true[,13] <- "same"
union_2_diff_0515true[union_2_diff_0515true[,11]=="same_tss",13] <- union_2_diff_0515true[union_2_diff_0515true[,11]=="same_tss",10]
union_2_diff_0515true[union_2_diff_0515true[,11]=="same_tes",12] <- union_2_diff_0515true[union_2_diff_0515true[,11]=="same_tes",10]
colnames(union_2_diff_0515true)[12:13] <- c("TSS_long?","TES_long?")

for(i in 1:(nrow(union_2_diff_0515true)/4)) {
  if(union_2_diff_0515true[4*i-3,11]=="total_diff") {
    if(union_2_diff_0515true[4*i-3,7]=="+") {
      if(union_2_diff_0515true[4*i-3,3]>union_2_diff_0515true[4*i-1,3]) {
        union_2_diff_0515true[4*i-3,12] <- "short"
        union_2_diff_0515true[4*i-2,12] <- "short"
        union_2_diff_0515true[4*i-1,12] <- "long"
        union_2_diff_0515true[4*i,12] <- "long"
      }
      if(union_2_diff_0515true[4*i-3,3]<union_2_diff_0515true[4*i-1,3]) {
        union_2_diff_0515true[4*i-3,12] <- "long"
        union_2_diff_0515true[4*i-2,12] <- "long"
        union_2_diff_0515true[4*i-1,12] <- "short"
        union_2_diff_0515true[4*i,12] <- "short"
      }
      if(union_2_diff_0515true[4*i-3,5]>union_2_diff_0515true[4*i-1,5]) {
        union_2_diff_0515true[4*i-3,13] <- "long"
        union_2_diff_0515true[4*i-2,13] <- "long"
        union_2_diff_0515true[4*i-1,13] <- "short"
        union_2_diff_0515true[4*i,13] <- "short"
      }
      if(union_2_diff_0515true[4*i-3,5]<union_2_diff_0515true[4*i-1,5]) {
        union_2_diff_0515true[4*i-3,13] <- "short"
        union_2_diff_0515true[4*i-2,13] <- "short"
        union_2_diff_0515true[4*i-1,13] <- "long"
        union_2_diff_0515true[4*i,13] <- "long"
      }
    }
    if(union_2_diff_0515true[4*i-3,7]=="-") {
      if(union_2_diff_0515true[4*i-3,3]>union_2_diff_0515true[4*i-1,3]) {
        union_2_diff_0515true[4*i-3,12] <- "long"
        union_2_diff_0515true[4*i-2,12] <- "long"
        union_2_diff_0515true[4*i-1,12] <- "short"
        union_2_diff_0515true[4*i,12] <- "short"
      }
      if(union_2_diff_0515true[4*i-3,3]<union_2_diff_0515true[4*i-1,3]) {
        union_2_diff_0515true[4*i-3,12] <- "short"
        union_2_diff_0515true[4*i-2,12] <- "short"
        union_2_diff_0515true[4*i-1,12] <- "long"
        union_2_diff_0515true[4*i,12] <- "long"
      }
      if(union_2_diff_0515true[4*i-3,5]>union_2_diff_0515true[4*i-1,5]) {
        union_2_diff_0515true[4*i-3,13] <- "short"
        union_2_diff_0515true[4*i-2,13] <- "short"
        union_2_diff_0515true[4*i-1,13] <- "long"
        union_2_diff_0515true[4*i,13] <- "long"
      }
      if(union_2_diff_0515true[4*i-3,5]<union_2_diff_0515true[4*i-1,5]) {
        union_2_diff_0515true[4*i-3,13] <- "long"
        union_2_diff_0515true[4*i-2,13] <- "long"
        union_2_diff_0515true[4*i-1,13] <- "short"
        union_2_diff_0515true[4*i,13] <- "short"
      }
    }
  }
}



for(i in 1:(nrow(union_2_diff_0515true)/4)) {
  a <- union_2_diff_0515true[(4*i-3):(4*i),]
  if(union_2_diff_0515true[4*i-3,12]=="same") {
    union_2_diff_0515true[(4*i-3):(4*i),14] <- 1
    union_2_diff_0515true[(4*i-3):(4*i),15] <- "same"
    union_2_diff_0515true[(4*i-3),16] <- log10(ceiling(a[a[,8]=="O" ,4])[1]/ceiling(a[a[,8]=="P" ,4])[1])
    union_2_diff_0515true[(4*i-3),17] <- union_2_diff_0515true[(4*i-3),16]
  }
  
  if(!union_2_diff_0515true[4*i-3,12]=="same") {
    union_2_diff_0515true[(4*i-3):(4*i),14] <- fisher.test(matrix(c(ceiling(a[a[,8]=="O" & a[,12]=="long",4]),
                                                                    ceiling(a[a[,8]=="O" & a[,12]=="short",4]),
                                                                    ceiling(a[a[,8]=="P" & a[,12]=="long",4]),
                                                                    ceiling(a[a[,8]=="P" & a[,12]=="short",4])),byrow = T,nrow = 2))$p.value
    if(ceiling(a[a[,8]=="O" & a[,12]=="long",4])/ceiling(a[a[,8]=="O" & a[,12]=="short",4])>ceiling(a[a[,8]=="P" & a[,12]=="long",4])/ceiling(a[a[,8]=="P" & a[,12]=="short",4])) {
      union_2_diff_0515true[(4*i-3):(4*i),15] <- "longer"
    }
    if(ceiling(a[a[,8]=="O" & a[,12]=="long",4])/ceiling(a[a[,8]=="O" & a[,12]=="short",4])<ceiling(a[a[,8]=="P" & a[,12]=="long",4])/ceiling(a[a[,8]=="P" & a[,12]=="short",4])) {
      union_2_diff_0515true[(4*i-3):(4*i),15] <- "shorter"
    }
    union_2_diff_0515true[(4*i-3),16] <- log10(ceiling(a[a[,8]=="O" & a[,12]=="long",4])/ceiling(a[a[,8]=="P" & a[,12]=="long",4]))
    union_2_diff_0515true[(4*i-3),17] <- log10(ceiling(a[a[,8]=="O" & a[,12]=="short",4])/ceiling(a[a[,8]=="P" & a[,12]=="short",4]))
  }
}


for(i in 1:(nrow(union_2_diff_0515true)/4)) {
  if(union_2_diff_0515true[4*i-3,14]<0.05) {
    if(union_2_diff_0515true[4*i-3,15]=="longer") {
      union_2_diff_0515true[4*i-3,18] <- "significant longer"
    }
    else union_2_diff_0515true[4*i-3,18] <- "significant shorter"
  }
  else union_2_diff_0515true[4*i-3,18] <- "no significance"
}


colnames(union_2_diff_0515true)[14:18] <- c("TSS_pvalue","TSS_choose_longer?","TSS_long_iso_O_devide_P","TSS_short_iso_O_devide_P","TSS_event")






#####TES


for(i in 1:(nrow(union_2_diff_0515true)/4)) {
  a <- union_2_diff_0515true[(4*i-3):(4*i),]
  if(union_2_diff_0515true[4*i-3,13]=="same") {
    union_2_diff_0515true[(4*i-3):(4*i),19] <- 1
    union_2_diff_0515true[(4*i-3):(4*i),20] <- "same"
    union_2_diff_0515true[(4*i-3),21] <- log10(ceiling(a[a[,8]=="O" ,6])[1]/ceiling(a[a[,8]=="P" ,6])[1])
    union_2_diff_0515true[(4*i-3),22] <- union_2_diff_0515true[(4*i-3),21]
  }
  
  if(!union_2_diff_0515true[4*i-3,13]=="same") {
    union_2_diff_0515true[(4*i-3):(4*i),19] <- fisher.test(matrix(c(ceiling(a[a[,8]=="O" & a[,13]=="long",6]),
                                                                    ceiling(a[a[,8]=="O" & a[,13]=="short",6]),
                                                                    ceiling(a[a[,8]=="P" & a[,13]=="long",6]),
                                                                    ceiling(a[a[,8]=="P" & a[,13]=="short",6])),byrow = T,nrow = 2))$p.value
    if(ceiling(a[a[,8]=="O" & a[,13]=="long",6])/ceiling(a[a[,8]=="O" & a[,13]=="short",6])>ceiling(a[a[,8]=="P" & a[,13]=="long",6])/ceiling(a[a[,8]=="P" & a[,13]=="short",6])) {
      union_2_diff_0515true[(4*i-3):(4*i),20] <- "longer"
    }
    if(ceiling(a[a[,8]=="O" & a[,13]=="long",6])/ceiling(a[a[,8]=="O" & a[,13]=="short",6])<ceiling(a[a[,8]=="P" & a[,13]=="long",6])/ceiling(a[a[,8]=="P" & a[,13]=="short",6])) {
      union_2_diff_0515true[(4*i-3):(4*i),20] <- "shorter"
    }
    union_2_diff_0515true[(4*i-3),21] <- log10(ceiling(a[a[,8]=="O" & a[,13]=="long",6])/ceiling(a[a[,8]=="P" & a[,13]=="long",6]))
    union_2_diff_0515true[(4*i-3),22] <- log10(ceiling(a[a[,8]=="O" & a[,13]=="short",6])/ceiling(a[a[,8]=="P" & a[,13]=="short",6]))
  }
}


for(i in 1:(nrow(union_2_diff_0515true)/4)) {
  if(union_2_diff_0515true[4*i-3,19]<0.05) {
    if(union_2_diff_0515true[4*i-3,20]=="longer") {
      union_2_diff_0515true[4*i-3,23] <- "significant longer"
    }
    else union_2_diff_0515true[4*i-3,23] <- "significant shorter"
  }
  else union_2_diff_0515true[4*i-3,23] <- "no significance"
}


colnames(union_2_diff_0515true)[19:23] <- c("TES_pvalue","TES_choose_longer?","TES_long_iso_O_devide_P","TES_short_iso_O_devide_P","TES_event")

union_2_diff_0515true[union_2_diff_0515true[,1]=="Litaf",20][1] <- "same"


a <- na.omit(union_2_diff_0515true)

union_2_diff_0515true[!union_2_diff_0515true[,1] %in% a[,1],]


nrow(a[a[,18]=="significant shorter" & a[,23]=="significant shorter",])

nrow(a[a[,18]=="significant shorter" & a[,23]=="no significance",])

nrow(a[a[,18]=="significant shorter" & a[,23]=="significant longer",])


nrow(a[a[,18]=="no significance" & a[,23]=="significant shorter",])

nrow(a[a[,18]=="no significance" & a[,23]=="no significance",])

nrow(a[a[,18]=="no significance" & a[,23]=="significant longer",])


nrow(a[a[,18]=="significant longer" & a[,23]=="significant shorter",])

nrow(a[a[,18]=="significant longer" & a[,23]=="no significance",])

nrow(a[a[,18]=="significant longer" & a[,23]=="significant longer",])



ggplot(a[a[,23]=="no significance",],aes(x=TSS_long_iso_O_devide_P,y=TSS_short_iso_O_devide_P,color=TSS_event))+
  geom_point()+
  labs(x="Ratio of isoform abundance(log10),longer",y="Ratio of isoform abundance(log10),shorter") +
  scale_color_manual(values = c("#BEBEBE","#00CDCD","#FF4500"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1,linetype = 1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5)) +
  geom_abline(slope = 1, intercept=0, na.rm = FALSE, show.legend = NA,linetype="dashed",size=1)+
  theme(legend.position=c(10,10))+   
  xlim(c(-3.5,3.5))+
  ylim(c(-3.5,3.5))


ggplot(a[a[,18]=="no significance",],aes(x=TES_long_iso_O_devide_P,y=TES_short_iso_O_devide_P,color=TES_event))+
  geom_point()+
  labs(x="Ratio of isoform abundance(log10),longer",y="Ratio of isoform abundance(log10),shorter") +
  scale_color_manual(values = c("#BEBEBE","#00CDCD","#FF4500"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1,linetype = 1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5)) +
  geom_abline(slope = 1, intercept=0, na.rm = FALSE, show.legend = NA,linetype="dashed",size=1)+
  theme(legend.position=c(10,10))+   
  xlim(c(-3.5,3.5))+
  ylim(c(-3.5,3.5))








a <- matrix(c(142,1612,343,818),nrow=2,byrow = T)
a <- log2(a+1)
pheatmap(a,
         col = colorRampPalette(colors = c("white","#483D8B"))(300),
         cluster_rows = F,
         cluster_col = F,
         display_numbers = F,
         show_rownames = T,
         show_colnames = T,
         legend = F,
         angle_col = "45",
         fontsize_row = 12,
         fontsize_col = 12,
         scale = "none")











