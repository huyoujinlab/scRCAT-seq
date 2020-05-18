############the first step is to call peaks:
########### on R:

options(stringsAsFactors = FALSE)
options(scipen = 100)

args=commandArgs(T)

args[1]

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

threshold=as.numeric(args[2])

library(CAGEr)


#### generate CAGEset 5'

a <- paste('myCAGEset',strsplit(filename,split = "_TKD")[[1]][1],'_5cap <- new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38",inputFiles = "',args[1],'",inputFilesType = "bed",sampleLabels = "sample")',sep = "")
eval(parse(text=a))
print(a)

a <- paste('getCTSS(','myCAGEset',strsplit(filename,split = "_TKD")[[1]][1],'_5cap,removeFirstG = FALSE,correctSystematicG = FALSE, nrCores = 40)',sep = "")
eval(parse(text=a))
print(a)





####5' single cell call peak
for(i in objects()[grep("5cap",objects())]) {
  a <- paste('librarysize <- as.numeric(librarySizes(',i,'))/1000000*threshold',sep = "")
  eval(parse(text=a))
  a <- paste('normalizeTagCount(',i,', method = "none",  T = 10^6)',sep = "")
  eval(parse(text=a))
  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = librarysize, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = T, nrCores = 40)',sep = "")
  eval(parse(text=a))
}




for(i in grep("5cap",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- tagClusters(',i,')[[1]]',sep = "")
  eval(parse(text=a))
  print(a)
  #a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- ',paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),'[',paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),'[,2]=="chr3",]',sep="")
  #eval(parse(text=a))
  #print(a)
}





for(i in grep("5cap",grep("tc_",objects(),value = TRUE),value = TRUE)) {
  a <- paste(paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tss',sep = ""),' <- data.frame(',i,'[,c(2,3,4,9,6,5)],stringsAsFactors = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)
}




for(i in grep("dominant",objects(),value = TRUE)) {
  a <- paste('write.table(',i,',"',paste('outdir/five_prime/peakfile/',i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
}

colnames(temp) <- c("V1","V2","V3","V4","V5","V6")

save.image(paste("outdir/five_prime/peakfile/temp_",strsplit(grep("dominant",objects(),value = TRUE),split="_5cap")[[1]][1],".RData",sep=""))

