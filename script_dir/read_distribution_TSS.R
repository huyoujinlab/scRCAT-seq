options(stringsAsFactors = FALSE)
options(scipen = 100)

library(CAGEr)
library(rlist)
library(parallel)
library(dplyr)


args=commandArgs(T)

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

load(paste("outdir/five_prime/peakfile/temp_",strsplit(filename,split="_5cap")[[1]][1],".RData",sep=""))

args=commandArgs(T)

args[1]


dominant_tss_in_gene <- read.table(args[1],sep="\t")

temp <- read.table(args[2])

temp1 <- anti_join(temp,dominant_tss_in_gene[,1:6])

dominant_tss_in_gene <- rbind(dominant_tss_in_gene,
                              data.frame(temp1,V7="chr1",V8=10,V9=11,V10="intergenic",V11="intergenic",V12="+"))

dominant_tss_in_gene <- dominant_tss_in_gene[,c(1,2,3,10,6)]

colnames(dominant_tss_in_gene) <- c("chr","start","end","gene","strand")

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

a <- paste('tc <- merge(',paste('tc_',strsplit(filename,split = "_domin")[[1]][1],sep = ''),',dominant_tss_in_gene,by=c("chr","start","end","strand"))',sep = "")
eval(parse(text=a))
print(a)
tc[,6] <- as.numeric(tc[,6])


tc[,c(4,5)] <- tc[,c(5,4)]

colnames(tc)[4:5] <- c("cluster","strand")

tc <- tc[,c(10,1,2,3,5,6,7,8,9)]


smart <- data.frame(gene=tc$gene,smart_seq2_read_count=1,smart_seq2_RPM=1)
smart <- unique(smart)
tc <- merge(tc,smart,by="gene")


transcript_length <- data.frame(gene=tc$gene,gene_length=1)
transcript_length <- unique(transcript_length)
tc <- merge(tc,transcript_length,by="gene")

colnames(tc) <- c("gene","chr","peak_start","peak_end","strand","peak_RPM","number_of_tss_in_peak","dominant_tss_position","dominant_tss_RPM","smart_seq2_read_count","smart_seq2_RPM","gene_length")


tc <- tc[,-7]

tc[,12] <- 1

all_gene_all_transcript_tss <- read.table("reference/gencode_hg38_all_gene_all_transcript_tss.bed",sep = "\t")


label <- function(tc) {
  temp <- all_gene_all_transcript_tss[all_gene_all_transcript_tss[,1]==tc[2] & all_gene_all_transcript_tss[,6]==tc[5] & abs(all_gene_all_transcript_tss[,3]-as.numeric(tc[7]))<=100,]
  if (nrow(temp)==0) {
    data.frame(xx=0)
  } else {
    data.frame(xx=1)
  }
}

cl <- makeCluster(30, type = "FORK")

test <- parApply(cl,tc,1,label)

tc_new <- list.rbind(test)

tc <- cbind(tc,tc_new)


colnames(tc)[c(12,13)] <- c("major_peak","annotated_peak")

tc <- unique(tc)



tc <- cbind(tc,data.frame(slope=1,correlation=1))


#####from caltss3



tc <- data.frame(tc[,c(1,2,3,4,5,6,7,8,11)],peak_width=(tc[,4]-tc[,3]),dominant_tss_RPM.smart2_seq_RPM=tc[,8]/tc[,10],tc[,c(14,15,13)])
tc[is.na(tc)] <- 0
tc[,11][is.infinite(tc[,11])] <- 0
colnames(tc)[11] <- "dominant_tss_RPM/smart_seq2_RPM"

gene_sum <- data.frame(gene=unique(tc[,1]),V2=1)

cal_percentage <- function(gene_sum) {
  temp <- sum(tc[tc[,1]==gene_sum[1],6])
  data.frame(sum=temp)
}

cl <- makeCluster(30, type = "FORK")

test <- parApply(cl,gene_sum,1,cal_percentage)

gene_sum <- data.frame(gene=gene_sum[,1],sum=list.rbind(test))

tc <- merge(tc,gene_sum,by="gene")

tc[,15] <- tc[,6]/tc[,15]

tc[,c(14,15)] <- tc[,c(15,14)]
colnames(tc)[c(14,15)] <- c("percentage","annotated_peak")

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

write.csv(tc,paste('outdir/five_prime/peakfile/tc_',strsplit(filename,split = "_5cap")[[1]][1],"_5cap.csv",sep = ""),quote = F,row.names = F)

write.table(tc[,c(2,3,4,1,6,5)],paste("outdir/five_prime/peakfile/temp_",strsplit(filename,split="_5cap")[[1]][1],"_tss.bed",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

