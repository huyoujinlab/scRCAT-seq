options(stringsAsFactors = FALSE)
options(scipen = 100)

args = commandArgs(T)

a <- read.csv(args[1])

output <- read.table(args[2])

gencode <- read.table(args[3])

a <- cbind(a,output[,7:20])

colnames(a)[664:677] <- c("strandard_5","strandard_10_6","strandard_10_7","strandard_10_8","strandard_15","strandard_20","strandard_30","strandard_50","A_percentage_5","A_percentage_10","A_percentage_15","A_percentage_20","A_percentage_30","A_percentage_50")

gencode <- gencode[,4:5]

gencode <- unique(gencode)

colnames(gencode) <- c("gene","gene_type")

gencode <- rbind(gencode,
                 data.frame(gene="intergenic",gene_type="intergenic"))


a <- merge(a,gencode,by="gene")

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

write.csv(a[,c(-9,-11)],paste("outdir/three_prime/peakfile/",filename,"_final.csv",sep=""),quote = F,row.names = F)
