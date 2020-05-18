

## supplementaryfig6a
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


load("supplementaryfig6.RData")


tss_df <- CTSStagCount(myCAGEsettss)

cor_tss_df <- cor(tss_df[,c(-1,-2,-3)],method = "pearson")
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

heatmap.2(cor_tss_df)

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
          #RowSideColors = c(rep("#0000CD",25),rep("#A52A2A",8),rep("#1C1C1C",22)),
          #ColSideColors = c(rep("#0000CD",25),rep("#A52A2A",8),rep("#1C1C1C",22)),
          #labCol = "",
          #labRow = "",
          density.info="none") 





























## supplementaryfig6b supplementaryfig6c









############################################## RC


####HEK293T  RC


#HEK293TtssRC <- read.csv("~/zjw/nc/novel20200428/result_read_counts/read_count_tss/tc_HEK293TRC5sccat_5cap_final_gene_rf_prediction.csv")
HEK293TtssRC <- HEK293TtssRC[HEK293TtssRC$model.prediction==1,]
HEK293TtssRC <- HEK293TtssRC[,c(1,2,3,4,6,5)]


#HEK293TtesRC <- read.csv("~/zjw/nc/novel20200428/result_read_counts/read_count_tes/tc_HEK293TRC3sccat_3tail_final_gene_rf_prediction.csv")
HEK293TtesRC <- HEK293TtesRC[HEK293TtesRC$model.prediction==1,]
HEK293TtesRC <- HEK293TtesRC[,c(1,2,3,4,6,5)]

HEK293TtssRC[,5] <- HEK293TtssRC[,5]/3982038*1000000
HEK293TtesRC[,5] <- HEK293TtesRC[,5]/28504810*1000000

####hESC RC

#hESCtssRC <- read.csv("~/zjw/nc/novel20200428/result_read_counts/read_count_tss/tc_hESCRC5sccat_5cap_final_gene_rf_prediction.csv")
hESCtssRC <- hESCtssRC[hESCtssRC$model.prediction==1,]
hESCtssRC <- hESCtssRC[,c(1,2,3,4,6,5)]


#hESCtesRC <- read.csv("~/zjw/nc/novel20200428/result_read_counts/read_count_tes/tc_hESCRC3sccat_3tail_final_gene_rf_prediction.csv")
hESCtesRC <- hESCtesRC[hESCtesRC$model.prediction==1,]
hESCtesRC <- hESCtesRC[,c(1,2,3,4,6,5)]


hESCtssRC[,5] <- hESCtssRC[,5]/6696786*1000000
hESCtesRC[,5] <- hESCtesRC[,5]/51061170*1000000












O_dominant_tss_in_gene_1 <- HEK293TtssRC

O_dominant_tes_in_gene_1 <- HEK293TtesRC

P_dominant_tss_in_gene_1 <- hESCtssRC

P_dominant_tes_in_gene_1 <- hESCtesRC


O_merge_1 <- merge(O_dominant_tss_in_gene_1,O_dominant_tes_in_gene_1,by="gene")
P_merge_1 <- merge(P_dominant_tss_in_gene_1,P_dominant_tes_in_gene_1,by="gene")



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

a_temp <- data.frame()
for (i in unique(union_2_same[,1])) {
  temp <- union_2_same[union_2_same[,1]==i,]
  tssvalue=log10(temp[temp[,8]=="O",4]/temp[temp[,8]=="P",4])
  tesvalue=log10(temp[temp[,8]=="O",6]/temp[temp[,8]=="P",6])
  temp <- data.frame(i,"chrX",1,1,1,1,"+","O",1,"short","same_tss","same","short",1,"same",tssvalue,tssvalue,"no significance",1,"same",tesvalue,tesvalue,"no significance")
  colnames(temp) <- colnames(a)
  a_temp <- rbind(a_temp,temp)
  
}

a_temp <- rbind(a,a_temp)

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



ggplot(a_temp[a_temp[,23]=="no significance",],aes(x=TSS_long_iso_O_devide_P,y=TSS_short_iso_O_devide_P,color=TSS_event))+
  geom_point()+
  labs(x="Log10 fold change of abundance of isoforms with distal TSSs",y="Log10 fold change of abundance of isoforms with proximal TSSs") +
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


ggplot(a_temp[a_temp[,18]=="no significance",],aes(x=TES_long_iso_O_devide_P,y=TES_short_iso_O_devide_P,color=TES_event))+
  geom_point()+
  labs(x="Log10 fold change of abundance of isoforms with distal TESs",y="Log10 fold change of abundance of isoforms with proximal TESs") +
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

















##supplementaryfig6d supplementaryfig6e


############################################## UMI
####HEK293T  UMI


#HEK293TtssUMI <- read.csv("~/zjw/nc/novel20200428/results/hESC_as_model/model_3/result/HEK293T_tss_threshold3/tc_HEK293TnofiltersccatUMI5threshold3_5cap_final_gene_rf_prediction.csv")
HEK293TtssUMI <- HEK293TtssUMI[HEK293TtssUMI$model.prediction==1,]
HEK293TtssUMI <- HEK293TtssUMI[,c(1,2,3,4,6,5)]


#HEK293TtesUMI <- read.csv("~/zjw/nc/novel20200428/results/hESC_as_model/model_3/result/HEK293T_tes_threshold3/tc_HEK293TnofiltersccatUMI3threshold3_3tail_final_gene_rf_prediction.csv")
HEK293TtesUMI <- HEK293TtesUMI[HEK293TtesUMI$model.prediction==1,]
HEK293TtesUMI <- HEK293TtesUMI[,c(1,2,3,4,6,5)]

HEK293TtssUMI[,5] <- HEK293TtssUMI[,5]/2378738*1000000
HEK293TtesUMI[,5] <- HEK293TtesUMI[,5]/22759497*1000000

####hESC UMI

#hESCtssUMI <- read.csv("~/zjw/nc/novel20200428/results/hESC_as_model/model_3/result/hESC_tss_threshold3/tc_hESCnofiltersccatUMI5threshold3_5cap_final_gene_rf_prediction.csv")
hESCtssUMI <- hESCtssUMI[hESCtssUMI$model.prediction==1,]
hESCtssUMI <- hESCtssUMI[,c(1,2,3,4,6,5)]


#hESCtesUMI <- read.csv("~/zjw/nc/novel20200428/results/hESC_as_model/model_3/result/hESC_tes_threshold3/tc_hESCnofiltersccatUMI3threshold3_3tail_final_gene_rf_prediction.csv")
hESCtesUMI <- hESCtesUMI[hESCtesUMI$model.prediction==1,]
hESCtesUMI <- hESCtesUMI[,c(1,2,3,4,6,5)]


hESCtssUMI[,5] <- hESCtssUMI[,5]/4235735*1000000
hESCtesUMI[,5] <- hESCtesUMI[,5]/40595523*1000000














O_dominant_tss_in_gene_1 <- HEK293TtssUMI

O_dominant_tes_in_gene_1 <- HEK293TtesUMI

P_dominant_tss_in_gene_1 <- hESCtssUMI

P_dominant_tes_in_gene_1 <- hESCtesUMI


O_merge_1 <- merge(O_dominant_tss_in_gene_1,O_dominant_tes_in_gene_1,by="gene")
P_merge_1 <- merge(P_dominant_tss_in_gene_1,P_dominant_tes_in_gene_1,by="gene")



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

a_temp <- data.frame()
for (i in unique(union_2_same[,1])) {
  temp <- union_2_same[union_2_same[,1]==i,]
  tssvalue=log10(temp[temp[,8]=="O",4]/temp[temp[,8]=="P",4])
  tesvalue=log10(temp[temp[,8]=="O",6]/temp[temp[,8]=="P",6])
  temp <- data.frame(i,"chrX",1,1,1,1,"+","O",1,"short","same_tss","same","short",1,"same",tssvalue,tssvalue,"no significance",1,"same",tesvalue,tesvalue,"no significance")
  colnames(temp) <- colnames(a)
  a_temp <- rbind(a_temp,temp)
  
}

a_temp <- rbind(a,a_temp)

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



ggplot(a_temp[a_temp[,23]=="no significance",],aes(x=TSS_long_iso_O_devide_P,y=TSS_short_iso_O_devide_P,color=TSS_event))+
  geom_point()+
  labs(x="Log10 fold change of abundance of isoforms with distal TSSs",y="Log10 fold change of abundance of isoforms with proximal TSSs") +
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


ggplot(a_temp[a_temp[,18]=="no significance",],aes(x=TES_long_iso_O_devide_P,y=TES_short_iso_O_devide_P,color=TES_event))+
  geom_point()+
  labs(x="Log10 fold change of abundance of isoforms with distal TESs",y="Log10 fold change of abundance of isoforms with proximal TESs") +
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
