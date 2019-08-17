options(stringsAsFactors = FALSE)
options(scipen=20)

library(CAGEr)
setwd("~/zjw/20190105/add_header")
tes <- read.table("D_tes.bed", header = FALSE, sep = "\t")
tes_plus <- tes[tes[,6]=="+",]
tes_minus <- tes[tes[,6]=="-",]
tes_plus <- tes_plus[,c(1,3,6)]
tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
tes <- rbind(tes_plus,tes_minus)
tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
write.table(tes,file = "D_tes.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")

tes <- read.table("O_tes.bed", header = FALSE, sep = "\t")
tes_plus <- tes[tes[,6]=="+",]
tes_minus <- tes[tes[,6]=="-",]
tes_plus <- tes_plus[,c(1,3,6)]
tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
tes <- rbind(tes_plus,tes_minus)
tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
write.table(tes,file = "O_tes.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")

tes <- read.table("P_tes.bed", header = FALSE, sep = "\t")
tes_plus <- tes[tes[,6]=="+",]
tes_minus <- tes[tes[,6]=="-",]
tes_plus <- tes_plus[,c(1,3,6)]
tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
tes <- rbind(tes_plus,tes_minus)
tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
write.table(tes,file = "P_tes.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")

myCAGEsetDtes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "D_tes.txt", inputFilesType = "bed",sampleLabels = "D")
myCAGEsetOtes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "O_tes.txt", inputFilesType = "bed",sampleLabels = "O")
myCAGEsetPtes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "P_tes.txt", inputFilesType = "bed",sampleLabels = "P")

getCTSS(myCAGEsetDtes,removeFirstG = FALSE,correctSystematicG = FALSE)
getCTSS(myCAGEsetOtes,removeFirstG = FALSE,correctSystematicG = FALSE)
getCTSS(myCAGEsetPtes,removeFirstG = FALSE,correctSystematicG = FALSE)



for(i in grep("tail",list.files(),value = T)) {
        tes <- read.table(i, header = FALSE, sep = "\t")
        tes_plus <- tes[tes[,6]=="+",]
        tes_minus <- tes[tes[,6]=="-",]
        tes_plus <- tes_plus[,c(1,3,6)]
        tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
        tes <- rbind(tes_plus,tes_minus)
        tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
        write.table(tes,file = "sample.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
        a <- paste('myCAGEset',strsplit(i,split = "[.]")[[1]][1],' <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "sample.txt", inputFilesType = "bed",sampleLabels = "sample")',sep = "")
        eval(parse(text=a))
        a <- paste('getCTSS(','myCAGEset',strsplit(i,split = "[.]")[[1]][1],',removeFirstG = FALSE,correctSystematicG = FALSE)',sep = "")
        eval(parse(text=a))
}

rm(a)
rm(tes)
rm(tes_plus)
rm(tes_minus)


save.image("20190630.RData")
