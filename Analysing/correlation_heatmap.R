options(stringsAsFactors = FALSE)
options(scipen = 100)
setwd("G:/CAGEr/CAGEr20190318isoform_heatmap/")
library(CAGEr)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(heatmap.plus)
#load("20190316Dpool,Opool,Ppool,single-cell.RData")

#rm(myCAGEsetDtss)
#rm(myCAGEsetOtss)
#rm(myCAGEsetPtss)
#rm(myCAGEsetDtes)
#rm(myCAGEsetOtes)
#rm(myCAGEsetPtes)
#rm(myCAGEsetD47_71_3tail)
#rm(myCAGEsetD47_71_5cap)


####5' 单细胞callpeak
#for(i in objects()[grep("5cap",objects())]) {
#  a <- paste('normalizeTagCount(',i,', method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)',sep = "\t")
#  eval(parse(text=a))
#  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 10, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
#  eval(parse(text=a))
#}


####3' 单细胞callpeak
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

##输出dominant plus
#for(i in grep("plus",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
#  eval(parse(text=a))
#}

##输出dominant minus
#for(i in grep("minus",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
#  eval(parse(text=a))
#}












####bedtools intersect

rm(list = ls())


for(i in grep("bed",grep("stream",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split = "[.]")[[1]][1],' <- read.table("',i,'",header = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)
}

save.image("correlation_heatmap.RData")

for(i in grep("stream",objects(),value = T)) {
  a <- paste(i,' <- ',i,'[,c(1,3,6,4)]',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste(i,' <- unique(',i,')',sep = "")
  eval(parse(text=a))
  print(a)
}






###正负合并
for(i in grep("plus",grep("stream",objects(),value = T),value = T)) {
  a <- paste(strsplit(i,split = "plus")[[1]][1],strsplit(i,split = "plus_")[[1]][2],' <- rbind(',i,',',strsplit(i,split = "plus")[[1]][1],'minus',strsplit(i,split = "plus")[[1]][2],')',sep = "")
  eval(parse(text=a))
  print(a)
}

#找回原来的值


for(i in grep("cap",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste('b <- CTSStagCount(',i,')',sep = "")
  eval(parse(text=a))
  a <- 'b <- sum(b[,4])'
  eval(parse(text=a))
  a <- paste(strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_upstream2k_and_genebody[,4] <- round(',strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_upstream2k_and_genebody[,4]/1000000*',b,')',sep = "")
  eval(parse(text=a))
}


for(i in grep("tail",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste('b <- CTSStagCount(',i,')',sep = "")
  eval(parse(text=a))
  a <- 'b <- sum(b[,4])'
  eval(parse(text=a))
  a <- paste(strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_gene_genebody_and_downstream20k[,4] <- round(',strsplit(i,split = "myCAGEset")[[1]][2],'_dominant_tss_gene_genebody_and_downstream20k[,4]/1000000*',b,')',sep = "")
  eval(parse(text=a))
}






##输出
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


pdf(file=paste("heat_cor_5cap.pdf", sep=""), height = 10, width = 10)

heatmap.2(cor_tss_df, ##对象
          Rowv=TRUE, Colv="Rowv", dendrogram='row', ##聚类对象，要保持一致
          trace='none', # trace可以给每个色块中添加一条线，与行平行或者与列平行。其与色块中心的距离代表了这个值被显示的比例
          scale = "none", # 标准化
          # 调整热图大小比例
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
legend("topright",      # location of the legend on the heatmap plot
       legend = c("DRG", "Oocyte D3", "Oocyte D4"), # category labels
       #col = c("#0000CD", "#1C1C1C", "#A52A2A"),  # color key
       fill = c("#0000CD", "#1C1C1C", "#A52A2A"),  # color key
       border = c("#0000CD", "#1C1C1C", "#A52A2A"))
dev.off() ##这个在图形出不来时有用

pdf(file=paste("heat_cor_3tail.pdf", sep=""), height = 10, width = 10)
heatmap.2(cor_tes_df, ##对象
          Rowv=TRUE, Colv="Rowv", dendrogram='row', ##聚类对象，要保持一致
          trace='none', # trace可以给每个色块中添加一条线，与行平行或者与列平行。其与色块中心的距离代表了这个值被显示的比例
          scale = "none", # 标准化
          # 调整热图大小比例
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

dev.off() ##这个在图形出不来时有用

























