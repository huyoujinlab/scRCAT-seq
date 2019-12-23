############the first step is to call peaks:
########### on R:

options(stringsAsFactors = FALSE)
options(scipen = 100)
library(CAGEr)


#### generate CAGEset 5'

for(i in grep("GGG",list.files(),value = T)) {
  a <- paste('myCAGEset',strsplit(i,split = "_TKD")[[1]][1],'_5cap <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "',i,'",inputFilesType = "bed",sampleLabels = "sample")',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('getCTSS(','myCAGEset',strsplit(i,split = "_TKD")[[1]][1],'_5cap,removeFirstG = FALSE,correctSystematicG = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)
}


#### generate CAGEset 3'

for(i in grep("A10",list.files(),value = T)) {
  tes <- read.table(i, header = FALSE, sep = "\t")
  tes_plus <- tes[tes[,6]=="+",]
  tes_minus <- tes[tes[,6]=="-",]
  tes_plus <- tes_plus[,c(1,3,6)]
  tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
  tes <- rbind(tes_plus,tes_minus)
  tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
  write.table(tes,file = "tes.bed",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
  a <- paste('myCAGEset',strsplit(i,split = "_TKD")[[1]][1],'_3tail <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "tes.bed",inputFilesType = "bed",sampleLabels = "sample")',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('getCTSS(','myCAGEset',strsplit(i,split = "_TKD")[[1]][1],'_3tail,removeFirstG = FALSE,correctSystematicG = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)

}

save.image("20191008lanecdef_CAT.RData")

####5' single cell call peak
for(i in objects()[grep("5cap",objects())]) {
  a <- paste('normalizeTagCount(',i,', method = "simpleTpm",  T = 10^6)',sep = "\t")
  eval(parse(text=a))
  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
  eval(parse(text=a))
}


for(i in objects()[grep("3tail",objects())]) {
  a <- paste('normalizeTagCount(',i,', method = "simpleTpm",  T = 10^6)',sep = "\t")
  eval(parse(text=a))
  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
  eval(parse(text=a))
}


for(i in grep("5cap",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- tagClusters(',i,')[[1]]',sep = "")
  eval(parse(text=a))
  print(a)
}


for(i in grep("3tail",grep("myCAGEset",objects(),value = T),value = T)) {
  a <- paste(paste('tc_',strsplit(i,split = "set")[[1]][2],sep = ""),' <- tagClusters(',i,')[[1]]',sep = "")
  eval(parse(text=a))
  print(a)
}



for(i in grep("5cap",grep("tc_",objects(),value = TRUE),value = TRUE)) {
  a <- paste(paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tss',sep = ""),' <- data.frame(',i,'[,c(2,3,4,9,6,5)],stringsAsFactors = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)
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
}




############the next step is to use bedtools to find peaks that locate in genes
######### on linux:
for i in `ls|grep "_dominant_tes.bed"`; do bedtools intersect -s -a ${i} -b gencode_mm10_all_gene_genebody_and_downstream2k.bed -wa -wb > ${i%%.*}_genebody_downstream2k.bed; done
for i in `ls|grep "_dominant_tss.bed"`; do bedtools intersect -s -a ${i} -b gencode_mm10_all_gene_ustream2k_and_genebody.bed -wa -wb > ${i%%.*}_upstream2k_and_genebody.bed; done


