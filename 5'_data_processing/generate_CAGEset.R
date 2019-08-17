options(stringsAsFactors = FALSE)
options(scipen=20)

library(CAGEr)
setwd("~/zjw/20190109/final_out")


myCAGEsetDtes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "D_tss.bed", inputFilesType = "bed",sampleLabels = "D")
myCAGEsetOtes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "O_tss.bed", inputFilesType = "bed",sampleLabels = "O")
myCAGEsetPtes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "P_tss.bed", inputFilesType = "bed",sampleLabels = "P")

getCTSS(myCAGEsetDtss,removeFirstG = FALSE,correctSystematicG = FALSE)
getCTSS(myCAGEsetOtss,removeFirstG = FALSE,correctSystematicG = FALSE)
getCTSS(myCAGEsetPtss,removeFirstG = FALSE,correctSystematicG = FALSE)



for(i in grep(".bed",list.files(),value = T)) {
        a <- paste('myCAGEset',strsplit(i,split = "[.]")[[1]][1],' <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "sample.txt", inputFilesType = "bed",sampleLabels = "sample")',sep = "")
        eval(parse(text=a))
        a <- paste('getCTSS(','myCAGEset',strsplit(i,split = "[.]")[[1]][1],',removeFirstG = FALSE,correctSystematicG = FALSE)',sep = "")
        eval(parse(text=a))
}




save.image("20190630.RData")
