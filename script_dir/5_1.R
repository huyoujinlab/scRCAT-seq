options(stringsAsFactors = FALSE)
options(scipen = 100)

args = commandArgs(T)

a <- read.csv(args[1])

a_bed <- data.frame(a[,1],a[,7]-1,a[,7],a[,c(6,6,5)])

write.table(a_bed,paste(args[1],"_for_gs.bed",sep=""),quote = F,row.names = F,col.names = F,sep = "\t")
