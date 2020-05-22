options(stringsAsFactors = FALSE)
options(scipen = 100)

args = commandArgs(T)

a <- read.csv(args[1])

output <- read.table(args[2])

gencode <- read.table(args[3])

a <- cbind(a,output[,7:12])

colnames(a)[164:169] <- c("3_G_3","3_G_2","5_G_3","5_G_4","G_percentage_3","G_percentage_5")

gencode <- gencode[,4:5]

gencode <- unique(gencode)

colnames(gencode) <- c("gene","gene_type")

gencode <- rbind(gencode,
                 data.frame(gene="intergenic",gene_type="intergenic"))

a <- merge(a,gencode,by="gene")

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

write.csv(a[,c(-9,-11)],paste("outdir/five_prime/peakfile/",filename,"_final.csv",sep=""),quote = F,row.names = F)
