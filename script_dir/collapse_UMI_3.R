options(stringsAsFactors = FALSE)
options(scipen = 100)

args=commandArgs(T)

bed <- read.table(args[1])

bed[,c(4,5)]=1

bed_plus <- bed[bed[,6]=="+",]
bed_minus <- bed[bed[,6]=="-",]


bed_plus <- data.frame(V1=bed_plus[,1],V2=bed_plus[,3]-1,bed_plus[,3:8])
bed_minus <- data.frame(bed_minus[,1:2],V3=bed_minus[,2]+1,bed_minus[,4:8])

bed <- rbind(bed_plus,bed_minus)

bed <- unique(bed)

write.table(bed,paste(args[1],'.collapse',sep = ""),quote = F,row.names = F,col.names = F,sep = "\t")
