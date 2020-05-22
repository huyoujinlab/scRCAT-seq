######### 

options(stringsAsFactors = FALSE)
options(scipen = 100)

args=commandArgs(T)




peak <- read.csv(args[1])

filename <- strsplit(args[1],split="/")[[1]][length(strsplit(args[1],split="/")[[1]])]

inpolydb <- read.table(paste("outdir/three_prime/peakfile/temp_",strsplit(strsplit(filename,split = "_3tail")[[1]][1],split = "tc_")[[1]][2],"_tes_in_polydb.bed",sep=""))

inpolydb <- unique(data.frame(inpolydb[,c(1,2,3,4,6)],isinPAS=1))

colnames(inpolydb)[c(1,2,3,4,5)] <- c("chr","peak_start","peak_end","gene","strand")

peak <- merge(peak,inpolydb,all = T,by = c("chr","peak_start","peak_end","gene","strand"))

peak[is.na(peak)] <- 0

peak <- data.frame(peak[,1:6],
                   dominant_tes_position=peak[,7],
                   peak[,8:16])
colnames(peak)[11] <- "dominant_tes_RPM/smart_seq2_RPM"

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

write.csv(peak,paste("outdir/three_prime/peakfile/",strsplit(filename,split = "[.]")[[1]][1],'_new.csv',sep = ""),quote = F,row.names = F)
