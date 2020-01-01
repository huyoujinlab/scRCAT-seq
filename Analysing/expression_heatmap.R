options(stringsAsFactors = FALSE)
options(scipen = 100)
library(reshape)

setwd("G:/CAGEr/CAGEr20190616/")

load("G:/CAGEr/CAGEr20190506usageheatmap/20190316Dpool,Opool,Ppool,single-cell.RData")


ENSM2symbol <- read.table("G:/annotation/gencode.vM18.annotation.sorted.all.gene.gtf",sep = "\t")

gene_id <- strsplit(ENSM2symbol[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

gene_name <- strsplit(ENSM2symbol[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

ENSM2symbol <- cbind(gene_id,gene_name)

normalizeTagCount(myCAGEsetDtss, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
normalizeTagCount(myCAGEsetOtss, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
normalizeTagCount(myCAGEsetPtss, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
normalizeTagCount(myCAGEsetDtes, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
normalizeTagCount(myCAGEsetOtes, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)
normalizeTagCount(myCAGEsetPtes, method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)

i <- 10
clusterCTSS(object = myCAGEsetDtss, method = "distclu", threshold = i, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)
clusterCTSS(object = myCAGEsetOtss, method = "distclu", threshold = i, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)
clusterCTSS(object = myCAGEsetPtss, method = "distclu", threshold = i, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)
clusterCTSS(object = myCAGEsetDtes, method = "distclu", threshold = i, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)
clusterCTSS(object = myCAGEsetOtes, method = "distclu", threshold = i, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)
clusterCTSS(object = myCAGEsetPtes, method = "distclu", threshold = i, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)



tc_Dtss <- tagClusters(myCAGEsetDtss)[[1]]
tc_Otss <- tagClusters(myCAGEsetOtss)[[1]]
tc_Ptss <- tagClusters(myCAGEsetPtss)[[1]]
tc_Dtes <- tagClusters(myCAGEsetDtes)[[1]]
tc_Otes <- tagClusters(myCAGEsetOtes)[[1]]
tc_Ptes <- tagClusters(myCAGEsetPtes)[[1]]


tc_Dtss <- tc_Dtss[,c(2,3,4,9,6,5)]
tc_Otss <- tc_Otss[,c(2,3,4,9,6,5)]
tc_Ptss <- tc_Ptss[,c(2,3,4,9,6,5)]
tc_Dtes <- tc_Dtes[,c(2,3,4,9,6,5)]
tc_Otes <- tc_Otes[,c(2,3,4,9,6,5)]
tc_Ptes <- tc_Ptes[,c(2,3,4,9,6,5)]


tc_Dtss_plus <- tc_Dtss[tc_Dtss[,6]=="+",]
tc_Otss_plus <- tc_Otss[tc_Otss[,6]=="+",]
tc_Ptss_plus <- tc_Ptss[tc_Ptss[,6]=="+",]
tc_Dtss_minus <- tc_Dtss[tc_Dtss[,6]=="-",]
tc_Otss_minus <- tc_Otss[tc_Otss[,6]=="-",]
tc_Ptss_minus <- tc_Ptss[tc_Ptss[,6]=="-",]

tc_Dtes_plus <- tc_Dtes[tc_Dtes[,6]=="+",]
tc_Otes_plus <- tc_Otes[tc_Otes[,6]=="+",]
tc_Ptes_plus <- tc_Ptes[tc_Ptes[,6]=="+",]
tc_Dtes_minus <- tc_Dtes[tc_Dtes[,6]=="-",]
tc_Otes_minus <- tc_Otes[tc_Otes[,6]=="-",]
tc_Ptes_minus <- tc_Ptes[tc_Ptes[,6]=="-",]


write.table(tc_Dtss_plus,file = "tc_Dtss_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Otss_plus,file = "tc_Otss_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptss_plus,file = "tc_Ptss_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Dtss_minus,file = "tc_Dtss_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Otss_minus,file = "tc_Otss_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptss_minus,file = "tc_Ptss_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")


write.table(tc_Dtes_plus,file = "tc_Dtes_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Otes_plus,file = "tc_Otes_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptes_plus,file = "tc_Ptes_plus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Dtes_minus,file = "tc_Dtes_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Otes_minus,file = "tc_Otes_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptes_minus,file = "tc_Ptes_minus.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")



##########bedtools intersect

#bedtools intersect -a tc_Dtss_plus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb > tc_Dtss_plus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Otss_plus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb > tc_Otss_plus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Ptss_plus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb > tc_Ptss_plus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Dtss_minus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb > tc_Dtss_minus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Otss_minus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb > tc_Otss_minus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Ptss_minus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb > tc_Ptss_minus_upstream_downstream_2k.bed

#bedtools intersect -a tc_Dtes_plus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb > tc_Dtes_plus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Otes_plus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb > tc_Otes_plus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Ptes_plus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb > tc_Ptes_plus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Dtes_minus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb > tc_Dtes_minus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Otes_minus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb > tc_Otes_minus_upstream_downstream_2k.bed
#bedtools intersect -a tc_Ptes_minus.bed -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb > tc_Ptes_minus_upstream_downstream_2k.bed



tc_Dtss_plus_upstream_downstream_2k <- read.table("tc_Dtss_plus_upstream_downstream_2k.bed")
tc_Otss_plus_upstream_downstream_2k <- read.table("tc_Otss_plus_upstream_downstream_2k.bed")
tc_Ptss_plus_upstream_downstream_2k <- read.table("tc_Ptss_plus_upstream_downstream_2k.bed")
tc_Dtss_minus_upstream_downstream_2k <- read.table("tc_Dtss_minus_upstream_downstream_2k.bed")
tc_Otss_minus_upstream_downstream_2k <- read.table("tc_Otss_minus_upstream_downstream_2k.bed")
tc_Ptss_minus_upstream_downstream_2k <- read.table("tc_Ptss_minus_upstream_downstream_2k.bed")

tc_Dtes_plus_upstream_downstream_2k <- read.table("tc_Dtes_plus_upstream_downstream_2k.bed")
tc_Otes_plus_upstream_downstream_2k <- read.table("tc_Otes_plus_upstream_downstream_2k.bed")
tc_Ptes_plus_upstream_downstream_2k <- read.table("tc_Ptes_plus_upstream_downstream_2k.bed")
tc_Dtes_minus_upstream_downstream_2k <- read.table("tc_Dtes_minus_upstream_downstream_2k.bed")
tc_Otes_minus_upstream_downstream_2k <- read.table("tc_Otes_minus_upstream_downstream_2k.bed")
tc_Ptes_minus_upstream_downstream_2k <- read.table("tc_Ptes_minus_upstream_downstream_2k.bed")



tc_Dtss_plus_upstream_downstream_2k <- tc_Dtss_plus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Otss_plus_upstream_downstream_2k <- tc_Otss_plus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Ptss_plus_upstream_downstream_2k <- tc_Ptss_plus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Dtss_minus_upstream_downstream_2k <- tc_Dtss_minus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Otss_minus_upstream_downstream_2k <- tc_Otss_minus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Ptss_minus_upstream_downstream_2k <- tc_Ptss_minus_upstream_downstream_2k[,c(1,2,3,10,5,6)]

tc_Dtes_plus_upstream_downstream_2k <- tc_Dtes_plus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Otes_plus_upstream_downstream_2k <- tc_Otes_plus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Ptes_plus_upstream_downstream_2k <- tc_Ptes_plus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Dtes_minus_upstream_downstream_2k <- tc_Dtes_minus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Otes_minus_upstream_downstream_2k <- tc_Otes_minus_upstream_downstream_2k[,c(1,2,3,10,5,6)]
tc_Ptes_minus_upstream_downstream_2k <- tc_Ptes_minus_upstream_downstream_2k[,c(1,2,3,10,5,6)]



######## select hihgest peak
for(i in objects()[grep("upstream_downstream_2k",objects())]) {
  b <- data.frame()
  d <- paste('for(j in unique(',i,'[,4])) {',
             'a <- ',i,'[',i,'[,4]==j,];',
             'a <- a[order(a[,5],decreasing = TRUE),];',
             'b <- rbind(b,a[1,])}',sep = "")
  eval(parse(text=d))
  d <- paste(i,' <- b',sep = "")
  eval(parse(text=d))
}


########## DP

tc_Dtss_plus_upstream_downstream_2k_intersect <- tc_Dtss_plus_upstream_downstream_2k[tc_Dtss_plus_upstream_downstream_2k[,4] %in% tc_Ptss_plus_upstream_downstream_2k[,4],]
tc_Ptss_plus_upstream_downstream_2k_intersect <- tc_Ptss_plus_upstream_downstream_2k[tc_Ptss_plus_upstream_downstream_2k[,4] %in% tc_Dtss_plus_upstream_downstream_2k[,4],]
tc_Dtss_minus_upstream_downstream_2k_intersect <- tc_Dtss_minus_upstream_downstream_2k[tc_Dtss_minus_upstream_downstream_2k[,4] %in% tc_Ptss_minus_upstream_downstream_2k[,4],]
tc_Ptss_minus_upstream_downstream_2k_intersect <- tc_Ptss_minus_upstream_downstream_2k[tc_Ptss_minus_upstream_downstream_2k[,4] %in% tc_Dtss_minus_upstream_downstream_2k[,4],]


tc_Dtes_plus_upstream_downstream_2k_intersect <- tc_Dtes_plus_upstream_downstream_2k[tc_Dtes_plus_upstream_downstream_2k[,4] %in% tc_Ptes_plus_upstream_downstream_2k[,4],]
tc_Ptes_plus_upstream_downstream_2k_intersect <- tc_Ptes_plus_upstream_downstream_2k[tc_Ptes_plus_upstream_downstream_2k[,4] %in% tc_Dtes_plus_upstream_downstream_2k[,4],]
tc_Dtes_minus_upstream_downstream_2k_intersect <- tc_Dtes_minus_upstream_downstream_2k[tc_Dtes_minus_upstream_downstream_2k[,4] %in% tc_Ptes_minus_upstream_downstream_2k[,4],]
tc_Ptes_minus_upstream_downstream_2k_intersect <- tc_Ptes_minus_upstream_downstream_2k[tc_Ptes_minus_upstream_downstream_2k[,4] %in% tc_Dtes_minus_upstream_downstream_2k[,4],]




write.table(tc_Dtss_plus_upstream_downstream_2k_intersect,"tc_Dtss_plus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptss_plus_upstream_downstream_2k_intersect,"tc_Ptss_plus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Dtss_minus_upstream_downstream_2k_intersect,"tc_Dtss_minus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptss_minus_upstream_downstream_2k_intersect,"tc_Ptss_minus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

write.table(tc_Dtes_plus_upstream_downstream_2k_intersect,"tc_Dtes_plus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptes_plus_upstream_downstream_2k_intersect,"tc_Ptes_plus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Dtes_minus_upstream_downstream_2k_intersect,"tc_Dtes_minus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(tc_Ptes_minus_upstream_downstream_2k_intersect,"tc_Ptes_minus_upstream_downstream_2k_intersect.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")



#bedtools subtract -a tc_Dtss_plus_upstream_downstream_2k_intersect.bed -A -b tc_Ptss_plus_upstream_downstream_2k_intersect.bed >tc_Dtss_plus_upstream_downstream_2k_intersect_unique.bed
#bedtools subtract -a tc_Ptss_plus_upstream_downstream_2k_intersect.bed -A -b tc_Dtss_plus_upstream_downstream_2k_intersect.bed >tc_Ptss_plus_upstream_downstream_2k_intersect_unique.bed
#bedtools subtract -a tc_Dtss_minus_upstream_downstream_2k_intersect.bed -A -b tc_Ptss_minus_upstream_downstream_2k_intersect.bed >tc_Dtss_minus_upstream_downstream_2k_intersect_unique.bed
#bedtools subtract -a tc_Ptss_minus_upstream_downstream_2k_intersect.bed -A -b tc_Dtss_minus_upstream_downstream_2k_intersect.bed >tc_Ptss_minus_upstream_downstream_2k_intersect_unique.bed

#bedtools subtract -a tc_Dtes_plus_upstream_downstream_2k_intersect.bed -A -b tc_Ptes_plus_upstream_downstream_2k_intersect.bed >tc_Dtes_plus_upstream_downstream_2k_intersect_unique.bed
#bedtools subtract -a tc_Ptes_plus_upstream_downstream_2k_intersect.bed -A -b tc_Dtes_plus_upstream_downstream_2k_intersect.bed >tc_Ptes_plus_upstream_downstream_2k_intersect_unique.bed
#bedtools subtract -a tc_Dtes_minus_upstream_downstream_2k_intersect.bed -A -b tc_Ptes_minus_upstream_downstream_2k_intersect.bed >tc_Dtes_minus_upstream_downstream_2k_intersect_unique.bed
#bedtools subtract -a tc_Ptes_minus_upstream_downstream_2k_intersect.bed -A -b tc_Dtes_minus_upstream_downstream_2k_intersect.bed >tc_Ptes_minus_upstream_downstream_2k_intersect_unique.bed


tc_Dtss_plus_unique <- read.table("tc_Dtss_plus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Ptss_plus_unique <- read.table("tc_Ptss_plus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Dtss_minus_unique <- read.table("tc_Dtss_minus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Ptss_minus_unique <- read.table("tc_Ptss_minus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Dtes_plus_unique <- read.table("tc_Dtes_plus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Ptes_plus_unique <- read.table("tc_Ptes_plus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Dtes_minus_unique <- read.table("tc_Dtes_minus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")
tc_Ptes_minus_unique <- read.table("tc_Ptes_minus_upstream_downstream_2k_intersect_unique.bed",sep = "\t")


tc_Dtss_unique <- rbind(tc_Dtss_plus_unique,tc_Dtss_minus_unique)
rm(tc_Dtss_plus_unique)
rm(tc_Dtss_minus_unique)
tc_Ptss_unique <- rbind(tc_Ptss_plus_unique,tc_Ptss_minus_unique)
rm(tc_Ptss_plus_unique)
rm(tc_Ptss_minus_unique)
tc_Dtes_unique <- rbind(tc_Dtes_plus_unique,tc_Dtes_minus_unique)
rm(tc_Dtes_plus_unique)
rm(tc_Dtes_minus_unique)
tc_Ptes_unique <- rbind(tc_Ptes_plus_unique,tc_Ptes_minus_unique)
rm(tc_Ptes_plus_unique)
rm(tc_Ptes_minus_unique)

a <- as.data.frame(table(c(tc_Dtss_unique[,4],tc_Ptss_unique[,4],tc_Dtes_unique[,4],tc_Ptes_unique[,4])))
a[,1] <- as.character(a[,1])
a <- a[a[,2] %in% c(4),]

b <- rbind(tc_Dtss_unique,tc_Ptss_unique,tc_Dtes_unique,tc_Ptes_unique)

intersect(rownames(no_sign),a[,1])





#--------------------------------------------single cell smart-seq2 DEseq2


#########single cell read count
for(i in grep("TKD",grep("count",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split = "_TKD")[[1]][1],'_count <- read.table("',i,'")',sep = "")
  eval(parse(text=a))
  print(a)
}

df <- data.frame(gene=D41_110_count[,1])
for(i in grep("_count",objects(),value = T)) {
  a <- paste('df <- cbind(df,data.frame(',strsplit(i,split = "_count"),'=',i,'[,2]))',sep = "")
  eval(parse(text=a))
  print(a)
}

df <- df[,c(1:20,25:36)]    ############ DP



for(i in 1:nrow(df)) {
  if(df[i,1] %in% ENSM2symbol[,1]) df[i,1] <- ENSM2symbol[ENSM2symbol[,1]==df[i,1],2]
}
df <- df[!df[,1] %in% df[,1][duplicated(df[,1])],]
rownames(df) <- df[,1]
df <- df[,-1]
df <- df[-(54029:54033),]



condition <- factor(c(rep("D",19),rep("P",12)))

dds <- DESeqDataSetFromMatrix(df, DataFrame(condition), design= ~ condition )

dds <- DESeq(dds)

res <- results(dds)

mcols(res)

summary(res)

a <- as.data.frame(res)

no_sign <- subset(res, padj >= 0.05)

no_sign <- as.data.frame(no_sign)




#--------------------------------------------






#########单细胞 read count  这里可以计算nlrp5、tet3、dnmt1
for(i in grep("TKD",grep("count",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split = "_TKD")[[1]][1],'_count <- read.table("',i,'")',sep = "")
  eval(parse(text=a))
  print(a)
}

df <- data.frame(gene=D41_110_count[,1])
for(i in grep("_count",objects(),value = T)) {
  a <- paste('df <- cbind(df,data.frame(',strsplit(i,split = "_count"),'=',i,'[,2]))',sep = "")
  eval(parse(text=a))
  print(a)
}

for(i in 1:nrow(df)) {
  if(df[i,1] %in% ENSM2symbol[,1]) df[i,1] <- ENSM2symbol[ENSM2symbol[,1]==df[i,1],2]
}
df <- df[!df[,1] %in% df[,1][duplicated(df[,1])],]
rownames(df) <- df[,1]
df <- df[,-1]
##write.table(df,"smart_seq2_single_cell_martix.txt",quote = FALSE,row.names = T,col.names = T,sep = "\t")

df <- read.table("smart_seq2_single_cell_martix.txt",header = T)
for(i in 1:ncol(df)) {
  df[,i] <- df[,i]/sum(df[,i])*1000000
}
df <- df[,20:35]












df <- read.table("smart_seq2_single_cell_martix.txt",header = T)  ###DP
for(i in 1:ncol(df)) {
  df[,i] <- df[,i]/sum(df[,i])*1000000
}
df <- df[,c(-20,-21,-22,-23)]  ###DP
boxplot <- melt(df[rownames(df) %in% "Grpel1",])
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

#####################smart-seq2 end




load("G:/CAGEr/CAGEr20190506usageheatmap/20190316Dpool,Opool,Ppool,single-cell.RData")





for(i in grep("my",objects(),value = T)) {
  a <- paste('normalizeTagCount(',i,', method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)',sep = "")
  print(a)
  eval(parse(text=a))
}

for(i in grep("5cap",grep("myCAGEsetD",objects(),value = T),value = T)) {
  a <- paste('CTSS_',strsplit(i,split="myCAGEset")[[1]][2],' <- CTSSnormalizedTpm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}


for(i in grep("5cap",grep("myCAGEsetP",objects(),value = T),value = T)) {
  a <- paste('CTSS_',strsplit(i,split="myCAGEset")[[1]][2],' <- CTSSnormalizedTpm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}

for(i in grep("^CTSS",objects(),value = T)) {
  a <- paste(i,' <- ',i,'[',i,'[,4]>10,]',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(i,' <- data.frame(',i,'[,1],',i,'[,2]-1,',i,'[,c(2,4,4,3)])',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(i,'_plus <- ',i,'[',i,'[,6]=="+",]',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(i,'_minus <- ',i,'[',i,'[,6]=="-",]',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('write.table(',i,'_plus,"',i,'_plus.bed"',',quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('write.table(',i,'_minus,"',i,'_minus.bed"',',quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  print(a)
  eval(parse(text=a))
}


for i in `ls|grep "^CTSS_D"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtss_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTSS_D"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptss_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done
for i in `ls|grep "^CTSS_P"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtss_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTSS_P"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptss_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done
for i in `ls|grep "^CTSS_D"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtss_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTSS_D"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptss_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done
for i in `ls|grep "^CTSS_P"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtss_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTSS_P"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptss_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done

for i in `ls|grep "^CTSS_"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb >${i%%.*}_upstream_downstream_2k.bed; done
for i in `ls|grep "^CTSS_"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb >${i%%.*}_upstream_downstream_2k.bed; done






#-------------------------------------------------------------------usage heatmap前----------------------------------------------------#

for(i in grep("CTSS",objects(),value = T)) {
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}


for(i in grep("_in_",grep("CTSS",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- read.table("',i,'")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- ',strsplit(i,split="[.]")[[1]][1],'[,c(1,2,3,10,5,6)]',sep = "")
  print(a)
  eval(parse(text=a))
}




###########正负合并及去掉无用对象

for(i in grep("plus",objects(),value = T)) {
  a <- paste(strsplit(i,split="_plus")[[1]][1],strsplit(i,split="_plus")[[1]][2],' <- rbind(',i,',',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
}


###########取合并峰
for(i in grep("CTSS",objects(),value = T)) {
  b <- data.frame()
  d <- paste('for(j in unique(',i,'[,4])) {',
             'a <- ',i,'[',i,'[,4]==j,];',
             'a <- a[order(a[,5],decreasing = TRUE),];',
             'a[1,5] <- sum(a[,5]);',
             'b <- rbind(b,a[1,])}',sep = "")
  eval(parse(text=d))
  print(d)
  d <- paste(i,' <- b',sep = "")
  eval(parse(text=d))
  print(d)
}


############改基因名字
for(i in grep("_D$",objects(),value = T)) {
  a <- paste(i,'[,4] <- paste(',i,'[,4],"_D")',sep = "")
  print(a)
  eval(parse(text=a))
}


for(i in grep("_P$",objects(),value = T)) {
  a <- paste(i,'[,4] <- paste(',i,'[,4],"_P")',sep = "")
  print(a)
  eval(parse(text=a))
}


####in_D和in_P合并
for(i in grep("_D$",objects(),value = T)) {
  a <- paste(strsplit(i,split="_in_")[[1]][1],' <- rbind(',i,',',strsplit(i,split="_in_")[[1]][1],'_in_P)',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',strsplit(i,split="_in_")[[1]][1],'_in_P)',sep = "")
  print(a)
  eval(parse(text=a))
}


########改列名
for(i in grep("CTSS",objects(),value = T)) {
  a <- paste('colnames(',i,')[5] <- "',strsplit(strsplit(i,split="CTSS_")[[1]][2],split="_5cap")[[1]],'"',sep = "")
  print(a)
  eval(parse(text=a))
}


########heatmap matrix
hm_df_5cap <- CTSS_D41_52_5cap[,c(4,5)]
for(i in 2:length(objects()[grep("CTSS",objects())])) {
  b <- objects()[grep("CTSS",objects())][i]
  a <- paste()
  a <- paste('hm_df_5cap <- merge(hm_df_5cap,',b,'[,c(4,5)],by="V10",all=TRUE)',sep = "")
  eval(parse(text=a))
}

hm_df_5cap[is.na(hm_df_5cap)] <- 0

rownames(hm_df_5cap) <- hm_df_5cap[,1]
hm_df_5cap <- hm_df_5cap[-1]
hm_df_5cap <- as.matrix(hm_df_5cap)

write.table(hm_df_5cap,"hm_df_5cap.bed",quote = FALSE,row.names = T,col.names = T,sep = "\t")

##########上面形成heatmap的过程可以不运行，直接就read.table
hm_df_5cap <- read.table("hm_df_5cap.bed",sep = "\t")

hm_df_5cap <- log10(hm_df_5cap+1)




################选择没有差异的基因
rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))

hm_df_5cap <- hm_df_5cap[rowname %in% rownames(no_sign),]


################选择没有差异的基因hou






##############去除只有一处的基因

double <- as.data.frame(table(rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))))
double <- as.character(double[double[,2]==2,1])
hm_df_5cap <- hm_df_5cap[unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1])) %in% double,]

##############去除只有一处的基因 hou







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

temp[temp=="D"] <- "isoform_1"  ###根据isoform分开

temp[temp=="P"] <- "isoform_2"  ###根据isoform分开




annotation_row = data.frame(isoform=factor(temp))

ann_colors <- list(celltype = c(Oocyte_D3 = "#1C1C1C", DRG = "#0000CD"),
                   isoform = c(isoform_1 = "#1B9E77", isoform_2 = "#D95F02"))



rownames(annotation_row) <- rownames(hm_df_5cap)   ###根据iso分开



pheatmap(hm_df_5cap,
         scale="row",
         #color = colorRampPalette(colors = c("#6E46AF","#FAF7C4","red"))(100),
         #color = colorRampPalette(colors = c("#6D49AE","#FAF7C4","#A1173D"))(250),
         col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(200)),
         #col=as.character(jdb_palette("brewer_green",type="continuous")),
         cluster_rows = F,
         fontsize_row = 5,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_rownames = F,
         show_colnames = F,
         border_color = F)


#-------------------------------------------------------------------usage heatmap后----------------------------------------------------#


################normal heatmap################前

for(i in grep("CTSS",objects(),value = T)) {
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}
for(i in grep("upstream",grep("CTSS",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- read.table("',i,'")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- ',strsplit(i,split="[.]")[[1]][1],'[,c(1,2,3,10,5,6)]',sep = "")
  print(a)
  eval(parse(text=a))
}


###########正负合并及去掉无用对象

for(i in grep("plus",objects(),value = T)) {
  a <- paste(strsplit(i,split="_plus")[[1]][1],strsplit(i,split="_plus")[[1]][2],' <- rbind(',i,',',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
}

###########取合并峰
for(i in grep("CTSS",objects(),value = T)) {
  b <- data.frame()
  d <- paste('for(j in unique(',i,'[,4])) {',
             'a <- ',i,'[',i,'[,4]==j,];',
             'a <- a[order(a[,5],decreasing = TRUE),];',
             'a[1,5] <- sum(a[,5]);',
             'b <- rbind(b,a[1,])}',sep = "")
  eval(parse(text=d))
  print(d)
  d <- paste(i,' <- b',sep = "")
  eval(parse(text=d))
  print(d)
}


########改列名
for(i in grep("CTSS",objects(),value = T)) {
  a <- paste('colnames(',i,')[5] <- "',strsplit(strsplit(i,split="CTSS_")[[1]][2],split="_5cap")[[1]][1],'"',sep = "")
  print(a)
  eval(parse(text=a))
}

########heatmap matrix
hm_df_5cap_normal <- CTSS_D41_52_5cap_upstream_downstream_2k[,c(4,5)]
for(i in 2:length(objects()[grep("CTSS",objects())])) {
  b <- objects()[grep("CTSS",objects())][i]
  a <- paste()
  a <- paste('hm_df_5cap_normal <- merge(hm_df_5cap_normal,',b,'[,c(4,5)],by="V10",all=TRUE)',sep = "")
  eval(parse(text=a))
}

hm_df_5cap_normal[is.na(hm_df_5cap_normal)] <- 0




rownames(hm_df_5cap_normal) <- hm_df_5cap_normal[,1]
hm_df_5cap_normal <- hm_df_5cap_normal[,-1]
hm_df_5cap_normal <- as.matrix(hm_df_5cap_normal)


#write.table(hm_df_5cap_normal,"hm_df_5cap_normal.bed",quote = FALSE,row.names = T,col.names = T,sep = "\t")


############可以直接运行read.table()
hm_df_5cap_normal <- read.table("hm_df_5cap_normal.bed",sep = "\t")
hm_df_5cap_normal <- as.matrix(hm_df_5cap_normal)
hm_df_5cap_normal <- log10(hm_df_5cap_normal+1)


##############为了跟usage heatmap一样的gene

hm_df_5cap_normal <- hm_df_5cap_normal[rownames(hm_df_5cap_normal) %in% as.character(as.data.frame(table(unlist(strsplit(rownames(hm_df_5cap),split = " _"))))[,1]),]

##############为了跟usage heatmap一样的gene后

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




################normal heatmap################后





































###3'


load("G:/CAGEr/CAGEr20190506usageheatmap/20190316Dpool,Opool,Ppool,single-cell.RData")





for(i in grep("my",objects(),value = T)) {
  a <- paste('normalizeTagCount(',i,', method = "simpleTpm", fitInRange = c(10,1000),alpha = 1.19, T = 10^6)',sep = "")
  print(a)
  eval(parse(text=a))
}

for(i in grep("3tail",grep("myCAGEsetD",objects(),value = T),value = T)) {
  a <- paste('CTES_',strsplit(i,split="myCAGEset")[[1]][2],' <- CTSSnormalizedTpm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}


for(i in grep("3tail",grep("myCAGEsetP",objects(),value = T),value = T)) {
  a <- paste('CTES_',strsplit(i,split="myCAGEset")[[1]][2],' <- CTSSnormalizedTpm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}

for(i in grep("^CTES",objects(),value = T)) {
  a <- paste(i,' <- ',i,'[',i,'[,4]>10,]',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(i,' <- data.frame(',i,'[,1],',i,'[,2]-1,',i,'[,c(2,4,4,3)])',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(i,'_plus <- ',i,'[',i,'[,6]=="+",]',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(i,'_minus <- ',i,'[',i,'[,6]=="-",]',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('write.table(',i,'_plus,"',i,'_plus.bed"',',quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('write.table(',i,'_minus,"',i,'_minus.bed"',',quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  print(a)
  eval(parse(text=a))
}


for i in `ls|grep "^CTES_D"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtes_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTES_D"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptes_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done
for i in `ls|grep "^CTES_P"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtes_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTES_P"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptes_plus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done
for i in `ls|grep "^CTES_D"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtes_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTES_D"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptes_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done
for i in `ls|grep "^CTES_P"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Dtes_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_D.bed; done
for i in `ls|grep "^CTES_P"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b tc_Ptes_minus_upstream_downstream_2k_intersect_unique.bed -wa -wb >${i%%.*}_in_P.bed; done

for i in `ls|grep "^CTES_"|grep "plus.bed$"`; do bedtools intersect -a ${i} -b gencode_mm10_all_gene_upstream_downstream_2k_plus.bed -wa -wb >${i%%.*}_upstream_downstream_2k.bed; done
for i in `ls|grep "^CTES_"|grep "minus.bed$"`; do bedtools intersect -a ${i} -b gencode_mm10_all_gene_upstream_downstream_2k_minus.bed -wa -wb >${i%%.*}_upstream_downstream_2k.bed; done






#-------------------------------------------------------------------usage heatmap前----------------------------------------------------#
for(i in grep("CTES",objects(),value = T)) {
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}

for(i in grep("_in_",grep("CTES",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- read.table("',i,'")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- ',strsplit(i,split="[.]")[[1]][1],'[,c(1,2,3,10,5,6)]',sep = "")
  print(a)
  eval(parse(text=a))
}




###########正负合并及去掉无用对象

for(i in grep("plus",objects(),value = T)) {
  a <- paste(strsplit(i,split="_plus")[[1]][1],strsplit(i,split="_plus")[[1]][2],' <- rbind(',i,',',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
}


###########取合并峰
for(i in grep("CTES",objects(),value = T)) {
  b <- data.frame()
  d <- paste('for(j in unique(',i,'[,4])) {',
             'a <- ',i,'[',i,'[,4]==j,];',
             'a <- a[order(a[,5],decreasing = TRUE),];',
             'a[1,5] <- sum(a[,5]);',
             'b <- rbind(b,a[1,])}',sep = "")
  eval(parse(text=d))
  print(d)
  d <- paste(i,' <- b',sep = "")
  eval(parse(text=d))
  print(d)
}


############改基因名字
for(i in grep("_D$",objects(),value = T)) {
  a <- paste(i,'[,4] <- paste(',i,'[,4],"_D")',sep = "")
  print(a)
  eval(parse(text=a))
}


for(i in grep("_P$",objects(),value = T)) {
  a <- paste(i,'[,4] <- paste(',i,'[,4],"_P")',sep = "")
  print(a)
  eval(parse(text=a))
}


####in_D和in_P合并
for(i in grep("_D$",objects(),value = T)) {
  a <- paste(strsplit(i,split="_in_")[[1]][1],' <- rbind(',i,',',strsplit(i,split="_in_")[[1]][1],'_in_P)',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',strsplit(i,split="_in_")[[1]][1],'_in_P)',sep = "")
  print(a)
  eval(parse(text=a))
}


########改列名
for(i in grep("CTES",objects(),value = T)) {
  a <- paste('colnames(',i,')[5] <- "',strsplit(strsplit(i,split="CTES_")[[1]][2],split="_5cap")[[1]],'"',sep = "")
  print(a)
  eval(parse(text=a))
}


########heatmap matrix
hm_df_5cap <- CTES_D41_52_3tail[,c(4,5)]
for(i in 2:length(objects()[grep("CTES",objects())])) {
  b <- objects()[grep("CTES",objects())][i]
  a <- paste()
  a <- paste('hm_df_5cap <- merge(hm_df_5cap,',b,'[,c(4,5)],by="V10",all=TRUE)',sep = "")
  eval(parse(text=a))
}

hm_df_5cap[is.na(hm_df_5cap)] <- 0

rownames(hm_df_5cap) <- hm_df_5cap[,1]
hm_df_5cap <- hm_df_5cap[-1]
hm_df_5cap <- as.matrix(hm_df_5cap)
hm_df_5cap <- log10(hm_df_5cap+1)




################选择没有差异的基因
rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))

hm_df_5cap <- hm_df_5cap[rowname %in% rownames(no_sign),]


################选择没有差异的基因hou

##############去除只有一处的基因

double <- as.data.frame(table(rowname <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1]))))
double <- as.character(double[double[,2]==2,1])
hm_df_5cap <- hm_df_5cap[unlist(lapply(strsplit(rownames(hm_df_5cap),split = " _"), function(x) x[1])) %in% double,]

##############去除只有一处的基因 hou





a <- pheatmap(hm_df_5cap,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)

summary(a)



order_row= c(seq(1,nrow(hm_df_5cap),2),seq(2,nrow(hm_df_5cap),2))
hm_df_5cap <- hm_df_5cap[order_row,]



annotation_col = data.frame(celltype=factor(rep(c("DRG", "Oocyte_D3"), c(25,22))))

rownames(annotation_col) <- a$tree_col$labels[a$tree_col$order]


temp <- unlist(lapply(strsplit(rownames(hm_df_5cap),split = "_"), function(x) x[2]))

temp[temp=="D"] <- "isoform_1"

temp[temp=="P"] <- "isoform_2"



annotation_row = data.frame(isoform=factor(temp))



################下面2条命令，2选1！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
rownames(annotation_row) <- rownames(hm_df_5cap)


ann_colors <- list(celltype = c(Oocyte_D3 = "#1C1C1C", DRG = "#0000CD"),
                   isoform = c(isoform_1 = "#1B9E77", isoform_2 = "#D95F02"))



pheatmap(hm_df_5cap,
         scale="row",
         #color = colorRampPalette(colors = c("blue","white","red"))(100),
         col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),
         cluster_rows = F,
         fontsize_row = 5,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_rownames = F,
         show_colnames = F,
         border_color = F)


#-------------------------------------------------------------------usage heatmap后----------------------------------------------------#



################normal heatmap################前

for(i in grep("CTES",objects(),value = T)) {
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
}
for(i in grep("upstream",grep("CTES",list.files(),value = T),value = T)) {
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- read.table("',i,'")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste(strsplit(i,split="[.]")[[1]][1],' <- ',strsplit(i,split="[.]")[[1]][1],'[,c(1,2,3,10,5,6)]',sep = "")
  print(a)
  eval(parse(text=a))
}


###########正负合并及去掉无用对象

for(i in grep("plus",objects(),value = T)) {
  a <- paste(strsplit(i,split="_plus")[[1]][1],strsplit(i,split="_plus")[[1]][2],' <- rbind(',i,',',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',i,')',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('rm(',strsplit(i,split="_plus")[[1]][1],'_minus',strsplit(i,split="_plus")[[1]][2],')',sep = "")
  print(a)
  eval(parse(text=a))
}

###########取合并峰
for(i in grep("CTES",objects(),value = T)) {
  b <- data.frame()
  d <- paste('for(j in unique(',i,'[,4])) {',
             'a <- ',i,'[',i,'[,4]==j,];',
             'a <- a[order(a[,5],decreasing = TRUE),];',
             'a[1,5] <- sum(a[,5]);',
             'b <- rbind(b,a[1,])}',sep = "")
  eval(parse(text=d))
  print(d)
  d <- paste(i,' <- b',sep = "")
  eval(parse(text=d))
  print(d)
}


########改列名
for(i in grep("CTES",objects(),value = T)) {
  a <- paste('colnames(',i,')[5] <- "',strsplit(strsplit(i,split="CTES_")[[1]][2],split="_3tail")[[1]][1],'"',sep = "")
  print(a)
  eval(parse(text=a))
}

########heatmap matrix
hm_df_5cap_normal <- CTES_D41_52_3tail_upstream_downstream_2k[,c(4,5)]
for(i in 2:length(objects()[grep("CTES",objects())])) {
  b <- objects()[grep("CTES",objects())][i]
  a <- paste('hm_df_5cap_normal <- merge(hm_df_5cap_normal,',b,'[,c(4,5)],by="V10",all=TRUE)',sep = "")
  eval(parse(text=a))
}

hm_df_5cap_normal[is.na(hm_df_5cap_normal)] <- 0


##############为了跟usage heatmap一样的gene

hm_df_5cap_normal <- hm_df_5cap_normal[hm_df_5cap_normal[,1] %in% as.character(as.data.frame(table(unlist(strsplit(rownames(hm_df_5cap),split = " _"))))[,1]),]

##############为了跟usage heatmap一样的gene后

rownames(hm_df_5cap_normal) <- hm_df_5cap_normal[,1]
hm_df_5cap_normal <- hm_df_5cap_normal[,-1]
hm_df_5cap_normal <- as.matrix(hm_df_5cap_normal)
hm_df_5cap_normal <- log10(hm_df_5cap_normal+1)

a <- pheatmap(hm_df_5cap_normal,scale="row",col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(250)),border_color = F)
summary(a)
order_row = a$tree_row$order
hm_df_5cap_normal <- hm_df_5cap_normal[order_row,]


annotation_col = data.frame(celltype=factor(rep(c("Oocyte_D3","DRG"), c(22,25))))

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




################normal heatmap################后
