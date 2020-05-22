options(stringsAsFactors = FALSE)
options(scipen = 100)


library(rlist)
library(parallel)



args=commandArgs(T)

peak <- read.csv(args[1])

colnames(peak)[c(11,17,18,19)] <- c("dominant_tss_RPM/smart_seq2_RPM","BREu_motif","BREd_motif","TATA_motif")



for (i in 17:ncol(peak)) {
  peak[,i] <- as.character(peak[,i])
}


peak17 <- data.frame(peak[,17])
peak18 <- data.frame(peak[,18])
peak19 <- data.frame(peak[,19])

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

test17 <- list.rbind(test17)
test18 <- list.rbind(test18)
test19 <- list.rbind(test19)

colnames(test17) <- paste(colnames(peak)[17],"minus",seq(49,0),sep = "_")
colnames(test18) <- paste(colnames(peak)[18],"minus",seq(49,0),sep = "_")
colnames(test19) <- paste(colnames(peak)[19],"minus",seq(49,0),sep = "_")

peak <- cbind(peak[,c(1:16)],test17,test18,test19)

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]


write.csv(peak[,c(-11,-12,-13)],paste("outdir/five_prime/peakfile/",strsplit(filename,split = "_new")[[1]][1],"_final.csv",sep = ""),quote = F,row.names = F)

