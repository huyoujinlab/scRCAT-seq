############the first step is to call peaks:
########### on R:

options(stringsAsFactors = FALSE)
options(scipen = 100)

args=commandArgs(T)

args[1]

filename <- strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]

threshold=as.numeric(args[2])

library(CAGEr)


#### generate CAGEset 3'



tes <- read.table(args[1], header = FALSE, sep = "\t")
tes_plus <- tes[tes[,6]=="+",]
tes_minus <- tes[tes[,6]=="-",]
tes_plus <- tes_plus[,c(1,3,6)]
tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
tes <- rbind(tes_plus,tes_minus)
tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
write.table(tes,file = paste("tes_",strsplit(filename,split = "_TKD")[[1]][1],".bed",sep=""),quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")

a <- paste('myCAGEset',strsplit(filename,split = "_TKD")[[1]][1],'_3tail <- new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38",inputFiles = "',paste("tes_",strsplit(filename,split = "_TKD")[[1]][1],".bed",sep=""),'",inputFilesType = "bed",sampleLabels = "sample")',sep = "")
eval(parse(text=a))
print(a)

a <- paste('getCTSS(','myCAGEset',strsplit(filename,split = "_TKD")[[1]][1],'_3tail,removeFirstG = FALSE,correctSystematicG = FALSE,nrCores = 40)',sep = "")
eval(parse(text=a))
print(a)








####3' single cell call peak
for(i in objects()[grep("3tail",objects())]) {
  a <- paste('librarysize <- as.numeric(librarySizes(',i,'))/1000000*threshold',sep = "")
  eval(parse(text=a))
  a <- paste('normalizeTagCount(',i,', method = "none",  T = 10^6)',sep = "")
  eval(parse(text=a))
  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = librarysize, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 40)',sep = "")
  eval(parse(text=a))
}



for(i in grep("3tail",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- tagClusters(',i,')[[1]]',sep = "")
  eval(parse(text=a))
  print(a)
  #a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- ',paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),'[',paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),'[,2]=="chr3",]',sep="")
  #eval(parse(text=a))
  #print(a)
}





for(i in grep("3tail",grep("tc_",objects(),value = TRUE),value = TRUE)) {
  a <- paste(paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tes',sep = ""),' <- data.frame(',i,'[,c(2,3,4,9,6,5)],stringsAsFactors = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)
}




for(i in grep("dominant",objects(),value = TRUE)) {
  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
}

colnames(temp) <- c("V1","V2","V3","V4","V5","V6")

save.image(paste("temp_",strsplit(grep("dominant",objects(),value = TRUE),split="_3tail")[[1]][1],".RData",sep=""))
