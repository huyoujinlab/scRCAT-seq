###############the first args is tes_dominant_peak_file
###############the second args is smart-seq2 read count file

options(stringsAsFactors = FALSE)
options(scipen = 100)

library(CAGEr)

load("temp.RData")

args=commandArgs(T)

args[1]

dominant_tes_in_gene <- read.table(args[1],sep="\t")

dominant_tes_in_gene <- dominant_tes_in_gene[,c(1,2,3,10,6)]

colnames(dominant_tes_in_gene) <- c("chr","start","end","gene","strand")


a <- paste('tc <- merge(',paste('tc_',strsplit(args[1],split = "_domin")[[1]][1],sep = ''),',dominant_tes_in_gene,by=c("chr","start","end","strand"))',sep = "")
eval(parse(text=a))
print(a)



tc[,c(4,5)] <- tc[,c(5,4)]

colnames(tc)[4:5] <- c("cluster","strand")

tc <- tc[,c(10,1,2,3,5,6,7,8,9)]

smart <- read.table(args[2])

ENSM2symbol <- read.table("gencode.vM18.annotation.sorted.all.gene.gtf",sep = "\t")

gene_id <- strsplit(ENSM2symbol[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

gene_name <- strsplit(ENSM2symbol[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

ENSM2symbol <- cbind(gene_id,gene_name)


for(i in 1:nrow(smart)) {
  if(smart[i,1] %in% ENSM2symbol[,1]) smart[i,1] <- ENSM2symbol[ENSM2symbol[,1]==smart[i,1],2]
}

smart[,3] <- smart[,2]/sum(smart[,2])*10^6

colnames(smart) <- c("gene","smart_seq2_read_count","smart_seq2_RPM")

tc <- merge(tc,smart,by="gene")

gene <- rbind(read.table("gencode_mm10_all_gene_plus.bed"),
              read.table("gencode_mm10_all_gene_minus.bed"))


gene[,7] <- gene[,3]-gene[,2]

colnames(gene)[c(4,7)] <- c("gene","gene_length")

tc <- merge(tc,gene[,c(4,7)],by="gene")


colnames(tc) <- c("gene","chr","peak_start","peak_end","strand","peak_RPM","number_of_tes_in_peak","dominant_tes_position","dominant_tes_RPM","smart_seq2_read_count","smart_seq2_RPM","gene_length")


tc <- tc[,-7]

tc[,12] <- 1

all_gene_all_transcript_tes <- read.table("gencode_mm10_all_gene_all_transcript_tes.bed",sep = "\t")


for(j in 1:nrow(tc)) {
  temp <- all_gene_all_transcript_tes[all_gene_all_transcript_tes[,1]==tc[j,2] & all_gene_all_transcript_tes[,6]==tc[j,5] & abs(all_gene_all_transcript_tes[,3]-tc[j,7])<=100,]
  ifelse(nrow(temp)==0,tc[j,13] <- 0,tc[j,13] <- 1)
}


colnames(tc)[c(12,13)] <- c("major_peak","annotated_peak")

tc <- unique(tc)

write.csv(tc,paste('tc_',strsplit(args[1],split = '_do')[[1]][1],'_peak.csv',sep = ''),row.names = F,quote = F)

