options(stringsAsFactors = FALSE)
options(scipen = 100)

args=commandArgs(T)


a <- read.csv(args[1])
a <- data.frame(a[,c(1,2,3,4,5,6,8,11)],peak_length=(a[,4]-a[,3]),dominant_tss_RPM.smart2_seq_RPM=a[,8]/a[,10],a[,c(14,15,13)])
a[is.na(a)] <- 0
a[,10][is.infinite(a[,10])] <- 0
colnames(a)[10] <- "dominant_tes_RPM/smart_seq2_RPM"
for(i in 1:nrow(a)) {
  a[i,14] <- a[i,6]/sum(a[a[,1]==a[i,1],6])
}
a[,c(13,14)] <- a[,c(14,13)]
colnames(a)[c(13,14)] <- c("percentage","annotated_peak")
write.csv(a,paste(args[1],'.csv',sep = ""),quote = F,row.names = F)

write.table(a[,c(2,3,4,1,6,5)],"temp_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

