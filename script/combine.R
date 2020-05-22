options(stringsAsFactors = FALSE)
options(scipen = 100)

args = commandArgs(T)
TSS_file <- read.csv(args[1])
TES_file <- read.csv(args[2])
outdir <- args[3]
if (!dir.exists(outdir)) dir.create(outdir)

TSS_file <- TSS_file[TSS_file$model.prediction==1,]
TES_file <- TES_file[TES_file$model.prediction==1,]
TSS_file <- TSS_file[,c(2,7,7,1,6,5)]
TES_file <- TES_file[,c(2,7,7,1,6,5)]
TSS_file_new <- data.frame(unique(TSS_file[,c(1,4,6)]),Num_TSS=1,TSS_pos=1)
TES_file_new <- data.frame(unique(TES_file[,c(1,4,6)]),Num_TES=1,TES_pos=1)

for (i in 1:nrow(TSS_file_new)) {
  temp <- TSS_file[TSS_file[,4]==TSS_file_new[i,2],]
  temp <- temp[order(temp[,5],decreasing = T),]
  TSS_file_new[i,4] <- nrow(temp)
  TSS_file_new[i,5] <- temp[1,3]
}

for (i in 1:nrow(TES_file_new)) {
  temp <- TES_file[TES_file[,4]==TES_file_new[i,2],]
  temp <- temp[order(temp[,5],decreasing = T),]
  TES_file_new[i,4] <- nrow(temp)
  TES_file_new[i,5] <- temp[1,3]
}

combine <- merge(TSS_file_new,TES_file_new,by=c("gene","chr","strand"),all=F)
combine <- combine[,c(1,2,3,4,6,5,7)]
write.table(combine, file.path(outdir, "combine_result.tsv"), quote = F,row.names = F,col.names = T,sep = "\t")

