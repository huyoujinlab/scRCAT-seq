# fig3 supplementaryfig5
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



load("/home/zjw/zjw/nc/figuregithub/supplementaryfig5.RData")


### fig3c supplementaryfig5a


#getCTSS(myCAGEsettss,removeFirstG = FALSE,correctSystematicG = FALSE)
#getCTSS(myCAGEsettes,removeFirstG = FALSE,correctSystematicG = FALSE)

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
          RowSideColors = c(rep("#0000CD",18),rep("#A52A2A",8)),
          ColSideColors = c(rep("#0000CD",18),rep("#A52A2A",8)),
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
          RowSideColors = c(rep("#0000CD",18),rep("#A52A2A",8)),
          ColSideColors = c(rep("#0000CD",18),rep("#A52A2A",8)),
          labCol = "",
          labRow = "",
          density.info="none") 

dev.off() 




























#supplementaryfig5d

rownames(smart_seq2_readcount) <- smart_seq2_readcount[,1]
smart_seq2_readcount <- smart_seq2_readcount[,-1]
smart_seq2_readcount <- smart_seq2_readcount[,1:23]

smart_seq2_readcount <- smart_seq2_readcount[!rowSums(smart_seq2_readcount)==0,]

condition <- factor(c(rep("D",19),rep("O",4)))

dds <- DESeqDataSetFromMatrix(smart_seq2_readcount[1:28351,], DataFrame(condition), design= ~ condition)

dds <- DESeq(dds)

res <- results(dds)

mcols(res)

summary(res)

a <- as.data.frame(res)

no_sign <- subset(res, padj >= 0.05)

no_sign <- as.data.frame(no_sign)

smart_seq2_RPM <- smart_seq2_readcount

for (i in 1:ncol(smart_seq2_RPM)) {
  smart_seq2_RPM[,i] <- smart_seq2_RPM[,i]/sum(smart_seq2_RPM[,i])*1000000
}

boxplot <- melt(smart_seq2_RPM[rownames(smart_seq2_RPM) %in% "Arpc3",])
boxplot[,1] <- as.character(boxplot[,1])
boxplot[1:19,1] <- "DRG"
boxplot[20:23,1] <- "oocyte D4"
boxplot[,2] <- log10(boxplot[,2]+1)

my_comparison <- list(c("DRG","oocyte D4"))
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




###fig3d supplementaryfig5e supplementaryfig5f

hm_df_5cap[is.na(hm_df_5cap)] <- 0

rownames(hm_df_5cap) <- hm_df_5cap[,1]
hm_df_5cap <- hm_df_5cap[-1]
hm_df_5cap <- as.matrix(hm_df_5cap)


hm_df_5cap <- log10(hm_df_5cap+1)


################选择没有差异的基因
rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))

hm_df_5cap <- hm_df_5cap[rowname %in% rownames(no_sign),]


#### remove  1
double <- as.data.frame(table(rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))))
double <- as.character(double[double[,2]==2,1])
hm_df_5cap <- hm_df_5cap[unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1])) %in% double,]




boxplot3 <- melt(hm_df_5cap[rownames(hm_df_5cap) %in% "Arpc3 _D",])
boxplot4 <- melt(hm_df_5cap[rownames(hm_df_5cap) %in% "Arpc3 _O",])

boxplot3 <- data.frame(variable=rownames(boxplot3),value=boxplot3[,1])
boxplot4 <- data.frame(variable=rownames(boxplot4),value=boxplot4[,1])
boxplot3[,1] <- as.character(boxplot3[,1])
boxplot3[1:18,1] <- "DRG"
boxplot3[19:26,1] <- "oocyte D4"

boxplot4[,1] <- as.character(boxplot4[,1])
boxplot4[1:18,1] <- "DRG"
boxplot4[19:26,1] <- "oocyte D4"


my_comparison <- list(c("DRG","oocyte D4"))
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
#a <- pheatmap(hm_df_5cap)

summary(a)



order_row= c(seq(1,nrow(hm_df_5cap),2),seq(2,nrow(hm_df_5cap),2))
hm_df_5cap <- hm_df_5cap[order_row,]

annotation_col = data.frame(celltype=factor(rep(c("Oocyte_D4", "DRG"), c(8,18))))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]



temp <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = "_"), function(x) x[2]))  ###根据isoform分开

temp[temp=="D"] <- "isoform_1"  ###

temp[temp=="O"] <- "isoform_2"  ###




annotation_row = data.frame(isoform=factor(temp))

ann_colors <- list(celltype = c(Oocyte_D4 = "#A52A2A", DRG = "#0000CD"),
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





#supplementaryfig5b



hm_df_5cap_normal[is.na(hm_df_5cap_normal)] <- 0




rownames(hm_df_5cap_normal) <- hm_df_5cap_normal[,1]
hm_df_5cap_normal <- hm_df_5cap_normal[,-1]
hm_df_5cap_normal <- as.matrix(hm_df_5cap_normal)


hm_df_5cap_normal <- as.matrix(hm_df_5cap_normal)
hm_df_5cap_normal <- log10(hm_df_5cap_normal+1)


hm_df_5cap_normal <- hm_df_5cap_normal[rownames(hm_df_5cap_normal) %in% as.character(as.data.frame(table(unlist(strsplit(rownames(hm_df_5cap),split = " _"))))[,1]),]

##############为了跟usage heatmap一样的gene后

a <- pheatmap(hm_df_5cap_normal,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
#a <- pheatmap(hm_df_5cap_normal)
summary(a)
order_row = a$tree_row$order
hm_df_5cap_normal <- hm_df_5cap_normal[order_row,]


annotation_col = data.frame(
  celltype=factor(rep(c("DRG","Oocyte_D4", "DRG","Oocyte_D4"), c(2,4,16,4))))

ann_colors <- list(celltype = c(Oocyte_D4 = "#A52A2A", DRG = "#0000CD"))

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










### supplementaryfig5c

hm_df_3tail[is.na(hm_df_3tail)] <- 0

rownames(hm_df_3tail) <- hm_df_3tail[,1]
hm_df_3tail <- hm_df_3tail[-1]
hm_df_3tail <- as.matrix(hm_df_3tail)


hm_df_3tail <- log10(hm_df_3tail+1)


################选择没有差异的基因
rowname <- unlist(lapply(strsplit(rownames(hm_df_3tail),split = " _"), function(x) x[1]))

hm_df_3tail <- hm_df_3tail[rowname %in% rownames(no_sign),]


#### remove  1
double <- as.data.frame(table(rowname <- unlist(lapply(strsplit(rownames(hm_df_3tail),split = " _"), function(x) x[1]))))
double <- as.character(double[double[,2]==2,1])
hm_df_3tail <- hm_df_3tail[unlist(lapply(strsplit(rownames(hm_df_3tail),split = " _"), function(x) x[1])) %in% double,]











a <- pheatmap(hm_df_3tail,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
#a <- pheatmap(hm_df_3tail)

summary(a)



order_row= c(seq(1,nrow(hm_df_3tail),2),seq(2,nrow(hm_df_3tail),2))
hm_df_3tail <- hm_df_3tail[order_row,]

annotation_col = data.frame(celltype=factor(rep(c("Oocyte_D4", "DRG"), c(8,18))))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]



temp <- unlist(lapply(strsplit(rownames(hm_df_3tail),split = "_"), function(x) x[2]))  ###根据isoform分开

temp[temp=="D"] <- "isoform_1"  ###

temp[temp=="O"] <- "isoform_2"  ###




annotation_row = data.frame(isoform=factor(temp))

ann_colors <- list(celltype = c(Oocyte_D4 = "#A52A2A", DRG = "#0000CD"),
                   isoform = c(isoform_1 = "#1B9E77", isoform_2 = "#D95F02"))



rownames(annotation_row) <- rownames(hm_df_3tail)   ###



pheatmap(hm_df_3tail,
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









hm_df_3tail_normal[is.na(hm_df_3tail_normal)] <- 0




rownames(hm_df_3tail_normal) <- hm_df_3tail_normal[,1]
hm_df_3tail_normal <- hm_df_3tail_normal[,-1]
hm_df_3tail_normal <- as.matrix(hm_df_3tail_normal)


hm_df_3tail_normal <- as.matrix(hm_df_3tail_normal)
hm_df_3tail_normal <- log10(hm_df_3tail_normal+1)


hm_df_3tail_normal <- hm_df_3tail_normal[rownames(hm_df_3tail_normal) %in% as.character(as.data.frame(table(unlist(strsplit(rownames(hm_df_3tail),split = " _"))))[,1]),]

##############为了跟usage heatmap一样的gene后

a <- pheatmap(hm_df_3tail_normal,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
#a <- pheatmap(hm_df_3tail_normal)
summary(a)
order_row = a$tree_row$order
hm_df_3tail_normal <- hm_df_3tail_normal[order_row,]


annotation_col = data.frame(
  celltype=factor(rep(c("DRG","Oocyte_D4", "DRG"), c(6,8,12))))

ann_colors <- list(celltype = c(Oocyte_D4 = "#A52A2A", DRG = "#0000CD"))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]


pheatmap(hm_df_3tail_normal,
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








###supplementaryfig8f






## Ermp1 Nsf

ind <- data.frame()
for (i in grep("in_D",ls(),value = T)) {
  print(i)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  #print(a)
  print(temp[temp[,4]=="Nsf",5])   ### 
  if (length(temp[temp[,4]=="Nsf",5])==0) {
    value=0
  } else {
    value=temp[temp[,4]=="Nsf",5]
  }
  ind <- rbind(ind,
               data.frame(value=value,group=i))
}


ino <- data.frame()
for (i in grep("in_O",ls(),value = T)) {
  print(i)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  #print(a)
  print(temp[temp[,4]=="Nsf",5])
  if (length(temp[temp[,4]=="Nsf",5])==0) {
    value=0
  } else {
    value=temp[temp[,4]=="Nsf",5]
  }
  ino <- rbind(ino,
               data.frame(value=value,group=i))
}






boxplot <- rbind(ind,ino)


boxplot[,2] <- c(rep("DRG",18),rep("Oocyte D4",8),rep("DRG",18),rep("Oocyte D4",8))
boxplot[,3] <- c(rep("Isoform 2",26),rep("Isoform 1",26))

boxplot[,1] <- log10(boxplot[,1]+1)

boxplot[,3] <- factor(boxplot[,3],levels = c("Isoform 1","Isoform 2"))

my_comparison <- list(c("DRG","Oocyte D4"))
ggboxplot(boxplot, x="V3", y="value", fill = "group",palette = "jco",xlab="",legend = "right") + 
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  labs(y="log10 (RPM+1)")+
  #labs(y="log10 (absolute gene expression +1)")+
  stat_compare_means(aes(group=group),label="..p.format..", method = "wilcox.test") +
  scale_y_continuous(limits = c(0,3))






boxplot <- data.frame(value=c(3,2,1,0,0,0,0,0,0,0,0,0,0,0,0,1),
                      group=c(rep("Oocyte D3",6),rep("DRG",2),rep("Oocyte D3",6),rep("DRG",2)),
                      V3=c(rep("Isoform 1",8),rep("Isoform 2",8)))



boxplot[,3] <- factor(boxplot[,3],levels = c("Isoform 1","Isoform 2"))

my_comparison <- list(c("DRG","Oocyte D4"))
ggboxplot(boxplot, x="V3", y="value", fill = "group",palette = "jco",xlab="",legend = "right") + 
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  labs(y="log10 (RPM+1)")+
  #labs(y="log10 (absolute gene expression +1)")+
  stat_compare_means(aes(group=group),label="..p.format..", method = "wilcox.test") +
  scale_y_continuous(limits = c(0,3))
