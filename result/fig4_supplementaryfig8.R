# fig4 supplementaryfig8
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


load("fig4_supplementaryfig8.RData")






D80_5_01 <- D80_5_01[D80_5_01$model.prediction==1,]
D80_5_02 <- D80_5_02[D80_5_02$model.prediction==1,]
D80_5_03 <- D80_5_03[D80_5_03$model.prediction==1,]
D80_5_04 <- D80_5_04[D80_5_04$model.prediction==1,]
D80_5_05 <- D80_5_05[D80_5_05$model.prediction==1,]
D80_5_06 <- D80_5_06[D80_5_06$model.prediction==1,]




D80_5_01 <- data.frame(D80_5_01[,2],D80_5_01[,7]-1,D80_5_01[,c(7,1,6,5)])
D80_5_02 <- data.frame(D80_5_02[,2],D80_5_02[,7]-1,D80_5_02[,c(7,1,6,5)])
D80_5_03 <- data.frame(D80_5_03[,2],D80_5_03[,7]-1,D80_5_03[,c(7,1,6,5)])
D80_5_04 <- data.frame(D80_5_04[,2],D80_5_04[,7]-1,D80_5_04[,c(7,1,6,5)])
D80_5_05 <- data.frame(D80_5_05[,2],D80_5_05[,7]-1,D80_5_05[,c(7,1,6,5)])
D80_5_06 <- data.frame(D80_5_06[,2],D80_5_06[,7]-1,D80_5_06[,c(7,1,6,5)])











D80_3_01 <- D80_3_01[D80_3_01$model.prediction==1,]
D80_3_02 <- D80_3_02[D80_3_02$model.prediction==1,]
D80_3_03 <- D80_3_03[D80_3_03$model.prediction==1,]
D80_3_04 <- D80_3_04[D80_3_04$model.prediction==1,]
D80_3_05 <- D80_3_05[D80_3_05$model.prediction==1,]
D80_3_06 <- D80_3_06[D80_3_06$model.prediction==1,]




D80_3_01 <- data.frame(D80_3_01[,2],D80_3_01[,7]-1,D80_3_01[,c(7,1,6,5)])
D80_3_02 <- data.frame(D80_3_02[,2],D80_3_02[,7]-1,D80_3_02[,c(7,1,6,5)])
D80_3_03 <- data.frame(D80_3_03[,2],D80_3_03[,7]-1,D80_3_03[,c(7,1,6,5)])
D80_3_04 <- data.frame(D80_3_04[,2],D80_3_04[,7]-1,D80_3_04[,c(7,1,6,5)])
D80_3_05 <- data.frame(D80_3_05[,2],D80_3_05[,7]-1,D80_3_05[,c(7,1,6,5)])
D80_3_06 <- data.frame(D80_3_06[,2],D80_3_06[,7]-1,D80_3_06[,c(7,1,6,5)])


####################################################################
D80_5_01_new <- data.frame()
for (i in unique(D80_5_01[,4])) {
  temp <- D80_5_01[D80_5_01[,4]==i,]
  temp <- temp[order(temp[,5],decreasing = T),][1,]
  D80_5_01_new <- rbind(D80_5_01_new,temp)
}



D80_3_01_new <- data.frame()
for (i in unique(D80_3_01[,4])) {
  temp <- D80_3_01[D80_3_01[,4]==i,]
  temp <- temp[order(temp[,5],decreasing = T),][1,]
  D80_3_01_new <- rbind(D80_3_01_new,temp)
}



D80_5_06_new <- data.frame()
for (i in unique(D80_5_06[,4])) {
  temp <- D80_5_06[D80_5_06[,4]==i,]
  temp <- temp[order(temp[,5],decreasing = T),][1,]
  D80_5_06_new <- rbind(D80_5_06_new,temp)
}


D80_3_06_new <- data.frame()
for (i in unique(D80_3_06[,4])) {
  temp <- D80_3_06[D80_3_06[,4]==i,]
  temp <- temp[order(temp[,5],decreasing = T),][1,]
  D80_3_06_new <- rbind(D80_3_06_new,temp)
}

colnames(D80_5_01_new) <- c("V1","V2","V3","V4","V5","V6")
colnames(D80_5_06_new) <- c("V1","V2","V3","V4","V5","V6")
colnames(D80_3_01_new) <- c("V1","V2","V3","V4","V5","V6")
colnames(D80_3_06_new) <- c("V1","V2","V3","V4","V5","V6")


D80_5_all <- rbind(D80_5_01_new,D80_5_06_new)
D80_3_all <- rbind(D80_3_01_new,D80_3_06_new)

D80_5_all[,5] <- 1
D80_3_all[,5] <- 1

D80_5_all <- unique(D80_5_all)
D80_3_all <- unique(D80_3_all)

D80_5_all_new <- data.frame()
for (i in unique(D80_5_all[,4])) {
  temp <- D80_5_all[D80_5_all[,4]==i,]
  if (nrow(temp)>=2) {
    temp <- temp[order(temp[,5],decreasing = T),][1:2,]
    D80_5_all_new <- rbind(D80_5_all_new,
                           temp)
  }
}


D80_3_all_new <- data.frame()
for (i in unique(D80_3_all[,4])) {
  temp <- D80_3_all[D80_3_all[,4]==i,]
  if (nrow(temp)>=2) {
    temp <- temp[order(temp[,5],decreasing = T),][1:2,]
    D80_3_all_new <- rbind(D80_3_all_new,
                           temp)
  }
}

rm(D80_5_01_new)
rm(D80_3_01_new)
rm(D80_5_06_new)
rm(D80_3_06_new)
###############




D80_5_01_new <- D80_5_01
D80_3_01_new <- D80_3_01
D80_5_06_new <- D80_5_06
D80_3_06_new <- D80_3_06

D80_5_01_new[,5] <- D80_5_01_new[,5]/1234930*1000000
D80_3_01_new[,5] <- D80_3_01_new[,5]/1928470*1000000
D80_5_06_new[,5] <- D80_5_06_new[,5]/1308374*1000000
D80_3_06_new[,5] <- D80_3_06_new[,5]/1463590*1000000



###



O_dominant_tss_in_gene_1 <- D80_5_06_new

O_dominant_tes_in_gene_1 <- D80_3_06_new

P_dominant_tss_in_gene_1 <- D80_5_01_new

P_dominant_tes_in_gene_1 <- D80_3_01_new









O_merge_1 <- merge(O_dominant_tss_in_gene_1,O_dominant_tes_in_gene_1,by="gene")
P_merge_1 <- merge(P_dominant_tss_in_gene_1,P_dominant_tes_in_gene_1,by="gene")



colnames(O_merge_1) <- c("gene_ID","chr","V3","tss_cor","tss_tpm","strand","V7","V8","tes_cor","tes_tpm","V11")
colnames(P_merge_1) <- c("gene_ID","chr","V3","tss_cor","tss_tpm","strand","V7","V8","tes_cor","tes_tpm","V11")



O_merge_1 <- O_merge_1[,c(1,2,4,5,9,10,6)]
P_merge_1 <- P_merge_1[,c(1,2,4,5,9,10,6)]






#for(i in row(O_merge_1)) {
#  if(O_merge_1[i,7]=="+") {if(O_merge_1[i,3]>O_merge_1[i,5]) O_merge_1[i,7] <- "."}
#  else {if(O_merge_1[i,3]<O_merge_1[i,5]) O_merge_1[i,7] <- "."}
#}
#for(i in row(P_merge_1)) {
#  if(P_merge_1[i,7]=="+") {if(P_merge_1[i,3]>P_merge_1[i,5]) P_merge_1[i,7] <- "."}
#  else {if(P_merge_1[i,3]<P_merge_1[i,5]) P_merge_1[i,7] <- "."}
#}

for(i in row(O_merge_1)) {
  if(O_merge_1[i,7]=="+") {
    if(O_merge_1[i,3]>O_merge_1[i,5]) {
      O_merge_1[i,7] <- "."
    }
  } else {
    if(O_merge_1[i,3]<O_merge_1[i,5]) {
      O_merge_1[i,7] <- "."
    }
  }
}


for(i in row(P_merge_1)) {
  if(P_merge_1[i,7]=="+") {
    if(P_merge_1[i,3]>P_merge_1[i,5]) {
      P_merge_1[i,7] <- "."
    }
  } else {
    if(P_merge_1[i,3]<P_merge_1[i,5]) {
      P_merge_1[i,7] <- "."
    }
  }
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
    union_2_diff_true[4*i-2,4] <- 0.5
  }
  if(is.na(union_2_diff_true[4*i-2,5])) {
    union_2_diff_true[4*i-2,5] <- union_2_diff_true[4*i-3,5]
    union_2_diff_true[4*i-2,6] <- 0.5
  }
  if(is.na(union_2_diff_true[4*i,3])) {
    union_2_diff_true[4*i,3] <- union_2_diff_true[4*i-1,3]
    union_2_diff_true[4*i,4] <- 0.5
  }
  if(is.na(union_2_diff_true[4*i,5])) {
    union_2_diff_true[4*i,5] <- union_2_diff_true[4*i-1,5]
    union_2_diff_true[4*i,6] <- 0.5
  }
}






for(i in 1:(nrow(union_2_diff_true)/4)) {
  a <- abs(union_2_diff_true[4*i-3,3]-union_2_diff_true[4*i-3,5])
  b <- abs(union_2_diff_true[4*i-1,3]-union_2_diff_true[4*i-1,5])
  #print(a)
  #print(b)
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
  if((abs(a[1,3]-a[3,3]))<20) union_2_diff_true_same_tss <- rbind(union_2_diff_true_same_tss,a)
  if((abs(a[1,5]-a[3,5]))<20) union_2_diff_true_same_tes <- rbind(union_2_diff_true_same_tes,a)
  if((abs(a[1,5]-a[3,5]))>19 & ((abs(a[1,3]-a[3,3]))>19)) union_2_diff_true_total_diff <- rbind(union_2_diff_true_total_diff,a)
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

union_2_diff_0515true <- union_2_diff_0515true[-grep("^AC[0-9]",union_2_diff_0515true[,1]),]

union_2_diff_0515true <- union_2_diff_0515true[-grep("^AP[0-9]",union_2_diff_0515true[,1]),]

union_2_diff_0515true <- union_2_diff_0515true[-grep("^AL[0-9]",union_2_diff_0515true[,1]),]

union_2_diff_0515true <- union_2_diff_0515true[-grep("-",union_2_diff_0515true[,1]),]



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
    } else {
      union_2_diff_0515true[4*i-3,18] <- "significant shorter"
    }
    
  } else {
    union_2_diff_0515true[4*i-3,18] <- "no significance"
  }
  
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
    } else {
      union_2_diff_0515true[4*i-3,23] <- "significant shorter"
    }
    
  } else {
    union_2_diff_0515true[4*i-3,23] <- "no significance"
  }
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


#write.csv(a[,c(1,14,18,19,23)],'1520vs1015.csv')

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

a_temp <- rbind(a,a_temp)

ggplot(a_temp,aes(x=TSS_long_iso_O_devide_P,y=TSS_short_iso_O_devide_P,color=TSS_event))+
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


ggplot(a_temp,aes(x=TES_long_iso_O_devide_P,y=TES_short_iso_O_devide_P,color=TES_event))+
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



sixone <- a[,c(1,14,18,19,23)]





rm(D80_5_01_new)
rm(D80_3_01_new)
rm(D80_5_06_new)
rm(D80_3_06_new)










##distal proximal 

#####################################################################################################
cor_df <- data.frame()
hm_5cap <- data.frame()
hm_5cap_percentage <- data.frame()
for (gene in unique(D80_5_all_new[,4])) {
  
  
  D80_5_all_new_new <- D80_5_all_new[D80_5_all_new[,4]==gene,]
  
  if (D80_5_all_new_new[1,6]=="+") {
    D80_5_all_new_new[D80_5_all_new_new[,3]==min(D80_5_all_new_new[,3]),4]="distal"
    D80_5_all_new_new[D80_5_all_new_new[,3]==max(D80_5_all_new_new[,3]),4]="proximal"
  } else {
    D80_5_all_new_new[D80_5_all_new_new[,3]==min(D80_5_all_new_new[,3]),4]="proximal"
    D80_5_all_new_new[D80_5_all_new_new[,3]==max(D80_5_all_new_new[,3]),4]="distal"
  }
  
  distance <- D80_5_all_new_new[1,3]-D80_5_all_new_new[2,3]
  
  if (distance<=20) {
    window=20
  }
  
  if (distance>20 & distance<200) {
    window=20
  }
  
  if (distance>200) {
    window=100
  }
  
  
  temp_df_df <- data.frame(V1=c("distal","proximal"))
  for (i in grep("D80_5_[0-9]",ls(),value = T)) {
    a <- paste('temp <- ',i,sep = "")
    eval(parse(text=a))
    if (length(temp[temp[,4]==gene & abs(temp[,3]-D80_5_all_new_new[D80_5_all_new_new[,4]=="distal",3])<window,5])==0) {
      distalvalue <- 1
    } else {
      distalvalue <- sum(temp[temp[,4]==gene & abs(temp[,3]-D80_5_all_new_new[D80_5_all_new_new[,4]=="distal",3])<window,5])
    }
    if (length(temp[temp[,4]==gene & abs(temp[,3]-D80_5_all_new_new[D80_5_all_new_new[,4]=="proximal",3])<window,5])==0) {
      proximalvalue <- 1
    } else {
      proximalvalue <- sum(temp[temp[,4]==gene & abs(temp[,3]-D80_5_all_new_new[D80_5_all_new_new[,4]=="proximal",3])<window,5])
    }
    #print(distalvalue)
    #print(proximalvalue)
    temp_df_df <- cbind(temp_df_df,
                        data.frame(vk=c(distalvalue,proximalvalue)))
    
    
  }
  if (!1 %in% as.numeric(temp_df_df[1,-1]*temp_df_df[2,-1])) {
    temp_df_df <- t(temp_df_df)[-1,]
    print(temp_df_df)
    temp_df_df <- data.frame(temp_df_df)
    temp_df_df[,1] <- as.numeric(temp_df_df[,1])
    temp_df_df[,2] <- as.numeric(temp_df_df[,2])
    temp_df_df[,3] <- log10(temp_df_df[,1]/temp_df_df[,2])
    temp_df_df[,4] <- 1:6
    
    colnames(temp_df_df) <- c("distal","proximal","distaldevideproximal","value")
    
    
    
    
    
    cor_df <- rbind(cor_df,
                    data.frame(gene,
                               cor=cor(temp_df_df[,4],temp_df_df[,3]),
                               nozeronum=length(as.numeric(c(temp_df_df[,1],temp_df_df[,2]))[!as.numeric(c(temp_df_df[,1],temp_df_df[,2]))==1])))
    
    hm_5cap_temp <- t(data.frame(c(temp_df_df[,3])))
    colnames(hm_5cap_temp) <- c("1","2","3","4","5","6")
    rownames(hm_5cap_temp) <- gene
    hm_5cap <- rbind(hm_5cap,hm_5cap_temp)
    
    
    
    temp_df_df[temp_df_df[,2]==1,2]=0
    temp_df_df[,3] <- temp_df_df[,1]/(temp_df_df[,1]+temp_df_df[,2])
    hm_5cap_percentage_temp <- t(data.frame(c(temp_df_df[,3])))
    colnames(hm_5cap_percentage_temp) <- c("1","2","3","4","5","6")
    rownames(hm_5cap_percentage_temp) <- gene
    hm_5cap_percentage <- rbind(hm_5cap_percentage,hm_5cap_percentage_temp)
  }
  
}
#####################################################################################################

cor_cor <- cor_df 

cor_cor[is.na(cor_cor)] <- 0

cor_cor <- cor_cor[order(cor_cor[,2],decreasing = T),]




hm_5cap_new <- hm_5cap[rownames(hm_5cap) %in% c(sixone[!sixone[,3]=="no significance",1]),]
hm_5cap_new <- hm_5cap_new[rownames(hm_5cap_new) %in% cor_cor[,1],]


hm_5cap_new <- hm_5cap_new[!rownames(hm_5cap_new) %in% c("NDUFC1","AC010616.2","EGLN3","ATP6V1G2-DDX39B","AC093512.2","CERKL","IFITM1","IFITM2","AL365205.1","TOMM6","CST3","ATP6V1G1"),]


for (i in 1:nrow(hm_5cap_new)) {
  hm_5cap_new[i,7] <- cor_cor[cor_cor[,1]==rownames(hm_5cap_new)[i],2]
}

hm_5cap_new[,8] <- abs(hm_5cap_new[,7])


hm_5cap_new <- hm_5cap_new[order(hm_5cap_new[,7],decreasing = T),]
top10gene <- rownames(hm_5cap_new[1:5,])
hm_5cap_new <- hm_5cap_new[order(hm_5cap_new[,7],decreasing = F),]
top10gene <- c(top10gene,rownames(hm_5cap_new[1:5,]))


hm_5cap_new <- hm_5cap_new[rownames(hm_5cap_new) %in% top10gene,] #

hm_5cap_new <- hm_5cap_new[order(hm_5cap_new[,7],decreasing = T),1:6]



bk <- seq(-2,2,by=0.04)

pheatmap(hm_5cap_new,
         scale="row",
         #col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(200)),
         cluster_rows = F,
         cluster_cols = F,
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         #fontsize_row = 5,
         show_rownames = T,
         show_colnames = T,
         border_color = F)


1


























#####################################################################################################
cor_df_3 <- data.frame()
hm_3tail <- data.frame()
hm_3tail_percentage <- data.frame()
for (gene in unique(D80_3_all_new[,4])) {
  
  #gene="NDUFC1"
  
  D80_3_all_new_new <- D80_3_all_new[D80_3_all_new[,4]==gene,]
  
  if (D80_3_all_new_new[1,6]=="+") {
    D80_3_all_new_new[D80_3_all_new_new[,3]==min(D80_3_all_new_new[,3]),4]="proximal"
    D80_3_all_new_new[D80_3_all_new_new[,3]==max(D80_3_all_new_new[,3]),4]="distal"
  } else {
    D80_3_all_new_new[D80_3_all_new_new[,3]==min(D80_3_all_new_new[,3]),4]="distal"
    D80_3_all_new_new[D80_3_all_new_new[,3]==max(D80_3_all_new_new[,3]),4]="proximal"
  }
  
  
  distance <- D80_3_all_new_new[1,3]-D80_3_all_new_new[2,3]
  
  if (distance<=20) {
    window=20
  }
  
  if (distance>20 & distance<200) {
    window=20
  }
  
  if (distance>200) {
    window=100
  }
  
  temp_df_df <- data.frame(V1=c("distal","proximal"))
  for (i in grep("D80_3_[0-9]",ls(),value = T)) {
    a <- paste('temp <- ',i,sep = "")
    eval(parse(text=a))
    if (length(temp[temp[,4]==gene & abs(temp[,3]-D80_3_all_new_new[D80_3_all_new_new[,4]=="distal",3])<window,5])==0) {
      distalvalue <- 1
    } else {
      distalvalue <- sum(temp[temp[,4]==gene & abs(temp[,3]-D80_3_all_new_new[D80_3_all_new_new[,4]=="distal",3])<window,5])
    }
    if (length(temp[temp[,4]==gene & abs(temp[,3]-D80_3_all_new_new[D80_3_all_new_new[,4]=="proximal",3])<window,5])==0) {
      proximalvalue <- 1
    } else {
      proximalvalue <- sum(temp[temp[,4]==gene & abs(temp[,3]-D80_3_all_new_new[D80_3_all_new_new[,4]=="proximal",3])<window,5])
    }
    #print(distalvalue)
    #print(proximalvalue)
    temp_df_df <- cbind(temp_df_df,
                        data.frame(vk=c(distalvalue,proximalvalue)))
    
    
  }
  if (!1 %in% as.numeric(temp_df_df[1,-1]*temp_df_df[2,-1])) {
    temp_df_df <- t(temp_df_df)[-1,]
    print(temp_df_df)
    temp_df_df <- data.frame(temp_df_df)
    temp_df_df[,1] <- as.numeric(temp_df_df[,1])
    temp_df_df[,2] <- as.numeric(temp_df_df[,2])
    temp_df_df[,3] <- log10(temp_df_df[,1]/temp_df_df[,2])
    temp_df_df[,4] <- 1:6
    
    colnames(temp_df_df) <- c("distal","proximal","distaldevideproximal","value")
    
    
    
    
    ggplot(temp_df_df, aes(x = value, y = distaldevideproximal )) +  
      geom_point(size=0.8) + 
      geom_smooth(color="black",method = glm) +
      theme_classic()+
      labs(x="pseudotime",y="log10 (distal/proximal)")
    
    
    cor_df_3 <- rbind(cor_df_3,
                      data.frame(gene,
                                 cor=cor(temp_df_df[,4],temp_df_df[,3]),
                                 nozeronum=length(as.numeric(c(temp_df_df[,1],temp_df_df[,2]))[!as.numeric(c(temp_df_df[,1],temp_df_df[,2]))==1])
                      ))
    hm_3tail_temp <- t(data.frame(c(temp_df_df[,3])))
    colnames(hm_3tail_temp) <- c("1","2","3","4","5","6")
    rownames(hm_3tail_temp) <- gene
    
    hm_3tail <- rbind(hm_3tail,hm_3tail_temp)
    
    temp_df_df[temp_df_df[,2]==1,2]=0
    temp_df_df[,3] <- temp_df_df[,1]/(temp_df_df[,1]+temp_df_df[,2])
    hm_3tail_percentage_temp <- t(data.frame(c(temp_df_df[,3])))
    colnames(hm_3tail_percentage_temp) <- c("1","2","3","4","5","6")
    rownames(hm_3tail_percentage_temp) <- gene
    hm_3tail_percentage <- rbind(hm_3tail_percentage,hm_3tail_percentage_temp)
  }
  
}
#####################################################################################################

cor_cor <- merge(cor_df_3,cor_df_3,by="gene")


#cor_cor[,2:4] <- abs(cor_cor[,2:4])


cor_cor <- na.omit(cor_cor)








cor_cor <- cor_cor[order(cor_cor[,4],decreasing = T),]







hm_3tail_new <- hm_3tail[rownames(hm_3tail) %in% c(sixone[!sixone[,5]=="no significance",1]),]

#hm_3tail_new <- hm_3tail_new[-grep('-',rownames(hm_3tail_new)),]

#hm_3tail_new <- hm_3tail_new[-grep('^AC',rownames(hm_3tail_new)),]

#hm_3tail_new <- hm_3tail_new[-grep('^AL',rownames(hm_3tail_new)),]

hm_3tail_new <- hm_3tail_new[!rownames(hm_3tail_new) %in% c("NDUFC1","AC010616.2","EGLN3","ATP6V1G2-DDX39B","AC093512.2","CERKL","CERKL","NEUROD1","AL669918.1","TAP2"),]

hm_3tail_new <- hm_3tail_new[rownames(hm_3tail_new) %in% cor_cor[,1],]
#hm_3tail_new <- hm_3tail_new[rownames(hm_3tail_new) %in% cor_cor[cor_cor[,5]>=6,1],]  #dorpout
for (i in 1:nrow(hm_3tail_new)) {
  hm_3tail_new[i,7] <- cor_cor[cor_cor[,1]==rownames(hm_3tail_new)[i],4]
}



hm_3tail_new[,8] <- abs(hm_3tail_new[,7])


hm_3tail_new <- hm_3tail_new[order(hm_3tail_new[,7],decreasing = T),]
top10gene <- rownames(hm_3tail_new[1:5,])
hm_3tail_new <- hm_3tail_new[order(hm_3tail_new[,7],decreasing = F),]
top10gene <- c(top10gene,rownames(hm_3tail_new[1:5,]))



hm_3tail_new <- hm_3tail_new[rownames(hm_3tail_new) %in% top10gene,] #

hm_3tail_new <- hm_3tail_new[order(hm_3tail_new[,7],decreasing = T),1:6]


bk <- seq(-2,2,by=0.04)


pheatmap(hm_3tail_new,
         scale="row",
         #col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(200)),
         cluster_rows = F,
         cluster_cols = F,
         #fontsize_row = 5,
         legend_breaks=seq(-2,2,1),
         breaks = bk,
         show_rownames = T,
         show_colnames = T,
         border_color = F)









###supplementaryfig8e 



cor_cor <- cor_df 

cor_cor[is.na(cor_cor)] <- 0

cor_cor <- cor_cor[order(cor_cor[,2],decreasing = T),]




hm_5cap_new <- hm_5cap[rownames(hm_5cap) %in% c(sixone[!sixone[,3]=="no significance",1]),]
hm_5cap_new <- hm_5cap_new[rownames(hm_5cap_new) %in% cor_cor[,1],]


hm_5cap_new <- hm_5cap_new[!rownames(hm_5cap_new) %in% c("NDUFC1","AC010616.2","EGLN3","ATP6V1G2-DDX39B","AC093512.2","CERKL","IFITM1","IFITM2","AL365205.1","TOMM6","CST3","ATP6V1G1"),]


for (i in 1:nrow(hm_5cap_new)) {
  hm_5cap_new[i,7] <- cor_cor[cor_cor[,1]==rownames(hm_5cap_new)[i],2]
}

hm_5cap_new[,8] <- abs(hm_5cap_new[,7])


hm_5cap_new <- hm_5cap_new[order(hm_5cap_new[,7],decreasing = T),]
top10gene <- rownames(hm_5cap_new[1:5,])
hm_5cap_new <- hm_5cap_new[order(hm_5cap_new[,7],decreasing = F),]
top10gene <- c(top10gene,rownames(hm_5cap_new[1:5,]))


#hm_5cap_new <- hm_5cap_new[rownames(hm_5cap_new) %in% top10gene,] #

hm_5cap_new <- hm_5cap_new[order(hm_5cap_new[,7],decreasing = T),1:6]




df <- data.frame(V1=c(1,2,3,4,5,6),
                 V2=c(nrow(hm_5cap_new[hm_5cap_new[,1]<0.1 & hm_5cap_new[,1]>-0.12,]),
                      nrow(hm_5cap_new[hm_5cap_new[,2]<0.1 & hm_5cap_new[,2]>-0.12,]),
                      nrow(hm_5cap_new[hm_5cap_new[,3]<0.1 & hm_5cap_new[,3]>-0.12,]),
                      nrow(hm_5cap_new[hm_5cap_new[,4]<0.1 & hm_5cap_new[,4]>-0.12,]),
                      nrow(hm_5cap_new[hm_5cap_new[,5]<0.1 & hm_5cap_new[,5]>-0.12,]),
                      nrow(hm_5cap_new[hm_5cap_new[,6]<0.1 & hm_5cap_new[,6]>-0.12,]))
                 )



ggbarplot(df, x="V1", y="V2",color = "V1", fill = "V1",ylab="Number of genes with two major isoforms equally expressed",xlab="",legend = "right",add.params = list(binwidth = 0.02)) +
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13)) +
  #stat_compare_means(comparisons=my_comparisons,label = "p-value", method = "t.test",paired = T) +
  scale_y_continuous(expand = c(0,0))










###supplementaryfig8f




cor_cor <- merge(cor_df_3,cor_df_3,by="gene")


#cor_cor[,2:4] <- abs(cor_cor[,2:4])


cor_cor <- na.omit(cor_cor)








cor_cor <- cor_cor[order(cor_cor[,4],decreasing = T),]







hm_3tail_new <- hm_3tail[rownames(hm_3tail) %in% c(sixone[!sixone[,5]=="no significance",1]),]


hm_3tail_new <- hm_3tail_new[!rownames(hm_3tail_new) %in% c("NDUFC1","AC010616.2","EGLN3","ATP6V1G2-DDX39B","AC093512.2","CERKL","CERKL","NEUROD1","AL669918.1","TAP2"),]

hm_3tail_new <- hm_3tail_new[rownames(hm_3tail_new) %in% cor_cor[,1],]
for (i in 1:nrow(hm_3tail_new)) {
  hm_3tail_new[i,7] <- cor_cor[cor_cor[,1]==rownames(hm_3tail_new)[i],4]
}



hm_3tail_new[,8] <- abs(hm_3tail_new[,7])


hm_3tail_new <- hm_3tail_new[order(hm_3tail_new[,7],decreasing = T),]
top10gene <- rownames(hm_3tail_new[1:5,])
hm_3tail_new <- hm_3tail_new[order(hm_3tail_new[,7],decreasing = F),]
top10gene <- c(top10gene,rownames(hm_3tail_new[1:5,]))



#hm_3tail_new <- hm_3tail_new[rownames(hm_3tail_new) %in% top10gene,] #

hm_3tail_new <- hm_3tail_new[order(hm_3tail_new[,7],decreasing = T),1:6]








df <- data.frame(V1=c(1,2,3,4,5,6),
                 V2=c(nrow(hm_3tail_new[hm_3tail_new[,1]<0.1 & hm_3tail_new[,1]>-0.12,]),
                      nrow(hm_3tail_new[hm_3tail_new[,2]<0.1 & hm_3tail_new[,2]>-0.12,]),
                      nrow(hm_3tail_new[hm_3tail_new[,3]<0.1 & hm_3tail_new[,3]>-0.12,]),
                      nrow(hm_3tail_new[hm_3tail_new[,4]<0.1 & hm_3tail_new[,4]>-0.12,]),
                      nrow(hm_3tail_new[hm_3tail_new[,5]<0.1 & hm_3tail_new[,5]>-0.12,]),
                      nrow(hm_3tail_new[hm_3tail_new[,6]<0.1 & hm_3tail_new[,6]>-0.12,]))
                 )



ggbarplot(df, x="V1", y="V2",color = "V1", fill = "V1",ylab="Number of genes with two major isoforms equally expressed",xlab="",legend = "right",add.params = list(binwidth = 0.02)) +
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13)) +
  #stat_compare_means(comparisons=my_comparisons,label = "p-value", method = "t.test",paired = T) +
  scale_y_continuous(expand = c(0,0))

