######### the first args is output file generated by cal_slope_cattss.R script

options(stringsAsFactors = FALSE)
options(scipen = 100)

args=commandArgs(T)




peak <- read.csv(args[1])

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]


infantom <- read.table(paste("outdir/five_prime/peakfile/temp_",strsplit(strsplit(filename,split = "_5cap")[[1]][1],split = "tc_")[[1]][2],"_tss_in_FANTOM.bed",sep=""))


infantom <- unique(data.frame(infantom[,c(1,2,3,4,6)],isindsc=1))

colnames(infantom)[c(1,2,3,4,5)] <- c("chr","peak_start","peak_end","gene","strand")

peak <- merge(peak,infantom,all = T,by = c("chr","peak_start","peak_end","gene","strand"))

peak[is.na(peak)] <- 0

peak <- data.frame(peak[,1:6],
                   dominant_tss_position=peak[,7],
                   peak[,8:16])
colnames(peak)[11] <- "dominant_tss_RPM/smart_seq2_RPM"

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

write.csv(peak,paste("outdir/five_prime/peakfile/",strsplit(filename,split = "[.]")[[1]][1],'_new.csv',sep = ""),quote = F,row.names = F)
