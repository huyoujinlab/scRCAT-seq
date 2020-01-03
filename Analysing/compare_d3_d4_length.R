

options(stringsAsFactors = FALSE)
options(scipen = 100)
library(CAGEr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pheatmap)






load("compare_d3_d4_length.RData")













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








