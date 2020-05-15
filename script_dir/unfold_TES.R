options(stringsAsFactors = FALSE)
options(scipen = 100)


library(rlist)
library(parallel)



args=commandArgs(T)

peak <- read.csv(args[1])

colnames(peak)[11] <- c("dominant_tes_RPM/smart_seq2_RPM")



for (i in 17:ncol(peak)) {
  peak[,i] <- as.character(peak[,i])
}

peak17 <- data.frame(peak[,17])
peak18 <- data.frame(peak[,18])
peak19 <- data.frame(peak[,19])
peak20 <- data.frame(peak[,20])
peak21 <- data.frame(peak[,21])
peak22 <- data.frame(peak[,22])
peak23 <- data.frame(peak[,23])
peak24 <- data.frame(peak[,24])
peak25 <- data.frame(peak[,25])
peak26 <- data.frame(peak[,26])
peak27 <- data.frame(peak[,27])
peak28 <- data.frame(peak[,28])
peak29 <- data.frame(peak[,29])


unfold <- function(peak) {
  if (is.na(peak[1])) {
    data.frame(matrix(0,nrow = 1,ncol = 50))
  } else {
    temp <- data.frame(matrix(0,nrow = 1,ncol = 50))
    temp[1,(51+as.numeric(unlist(strsplit(peak[1],split = ";"))))] <- 1
    temp
  }
}

cl <- makeCluster(30, type = "FORK")

test17 <- parApply(cl,peak17,1,unfold)
test18 <- parApply(cl,peak18,1,unfold)
test19 <- parApply(cl,peak19,1,unfold)
test20 <- parApply(cl,peak20,1,unfold)
test21 <- parApply(cl,peak21,1,unfold)
test22 <- parApply(cl,peak22,1,unfold)
test23 <- parApply(cl,peak23,1,unfold)
test24 <- parApply(cl,peak24,1,unfold)
test25 <- parApply(cl,peak25,1,unfold)
test26 <- parApply(cl,peak26,1,unfold)
test27 <- parApply(cl,peak27,1,unfold)
test28 <- parApply(cl,peak28,1,unfold)
test29 <- parApply(cl,peak29,1,unfold)

test17 <- list.rbind(test17)
test18 <- list.rbind(test18)
test19 <- list.rbind(test19)
test20 <- list.rbind(test20)
test21 <- list.rbind(test21)
test22 <- list.rbind(test22)
test23 <- list.rbind(test23)
test24 <- list.rbind(test24)
test25 <- list.rbind(test25)
test26 <- list.rbind(test26)
test27 <- list.rbind(test27)
test28 <- list.rbind(test28)
test29 <- list.rbind(test29)


colnames(test17) <- paste(colnames(peak)[17],"minus",seq(49,0),sep = "_")
colnames(test18) <- paste(colnames(peak)[18],"minus",seq(49,0),sep = "_")
colnames(test19) <- paste(colnames(peak)[19],"minus",seq(49,0),sep = "_")
colnames(test20) <- paste(colnames(peak)[20],"minus",seq(49,0),sep = "_")
colnames(test21) <- paste(colnames(peak)[21],"minus",seq(49,0),sep = "_")
colnames(test22) <- paste(colnames(peak)[22],"minus",seq(49,0),sep = "_")
colnames(test23) <- paste(colnames(peak)[23],"minus",seq(49,0),sep = "_")
colnames(test24) <- paste(colnames(peak)[24],"minus",seq(49,0),sep = "_")
colnames(test25) <- paste(colnames(peak)[25],"minus",seq(49,0),sep = "_")
colnames(test26) <- paste(colnames(peak)[26],"minus",seq(49,0),sep = "_")
colnames(test27) <- paste(colnames(peak)[27],"minus",seq(49,0),sep = "_")
colnames(test28) <- paste(colnames(peak)[28],"minus",seq(49,0),sep = "_")
colnames(test29) <- paste(colnames(peak)[29],"minus",seq(49,0),sep = "_")

peak <- cbind(peak[,c(1:16)],test17,test18,test19,test20,test21,test22,test23,test24,test25,test26,test27,test28,test29)




write.csv(peak[,c(-11,-12,-13)],paste(strsplit(args[1],split = "_new")[[1]][1],"_final.csv",sep = ""),quote = F,row.names = F)
