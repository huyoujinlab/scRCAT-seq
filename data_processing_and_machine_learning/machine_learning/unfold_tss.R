options(stringsAsFactors = FALSE)
options(scipen = 100)





args=commandArgs(T)

peak <- read.csv(args[1])

#peak <- read.csv("tc_D44_52_5cap_peak.csv.csv.csv.csv.csv")

colnames(peak)[c(11,17,18,19)] <- c("dominant_tss_RPM/smart_seq2_RPM","BREu_motif","BREd_motif","TATA_motif")



for (i in 17:ncol(peak)) {
  peak[,i] <- as.character(peak[,i])
}

gene <- read.table("transcript_length.tsv")


temp4 <- data.frame()
for (i in 1:nrow(peak)) {
  #print(i)
  peak[i,9] <- gene[gene[,1]==peak[i,4],2]
  
  
  
  temp3 <- data.frame()
  for (k in 17:ncol(peak)) {
    if (is.na(peak[i,k])) {
      temp3 <- rbind(temp3,
                     data.frame(V1=paste(colnames(peak)[k],"minus",seq(49,0),sep = "_"),
                                V2=0)
      )
    } else {
      temp3 <- rbind(temp3,
                     data.frame(V1=paste(colnames(peak)[k],"minus",seq(49,0),sep = "_"),
                                V2=0)
      )
      temp3[temp3[,1] %in% grep(colnames(peak)[k],grep(paste(-as.numeric(unlist(strsplit(peak[i,k],split = ";"))),collapse = "|"),temp3[,1],value = T),value = T),2] <- 1
    }
    
  }
  temp4 <- rbind(temp4,
                 as.data.frame(t(temp3))[2,])
}

colnames(temp4) <- temp3[,1]

peak <- cbind(peak[,1:16],temp4)



write.csv(peak,paste(args[1],".csv",sep=""),quote = F,row.names = F)
