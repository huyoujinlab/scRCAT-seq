#tes <- read.table("3-4_L1_805D75.R1.fastq-common.out_withA10_trim_Aligned.out.sam_extract_uniquely_map.sam_add_header_sorted.bed", header = FALSE, sep = "\t")
#tes_plus <- tes[tes[,6]=="+",]
#tes_minus <- tes[tes[,6]=="-",]
#tes_plus <- tes_plus[,c(1,3,6)]
#tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
#tes <- rbind(tes_plus,tes_minus)
#tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
#write.table(tes,file = "tes.bed",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
#myCAGEset3_4_3tail <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "tes.bed",inputFilesType = "bed",sampleLabels = "sample")
#getCTSS(myCAGEset3_4_3tail,removeFirstG = FALSE,correctSystematicG = FALSE)

#myCAGEset3_4_5cap <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "3-4_L1_805D75.R2.fastq_with_tag.fq.trimed.fq.trimed.nextera_Aligned.out.sam_extract_uniquely_map.sam_add_header_sorted.bed",inputFilesType = "bed",sampleLabels = "sample")
#getCTSS(myCAGEset3_4_5cap,removeFirstG = FALSE,correctSystematicG = FALSE)



#for(i in grep("myCAGEset",objects(),value = T)) {
#  a <- paste('normalizeTagCount(',i,', method = "simpleTpm")',sep = "")
#  print(a)
#  eval(parse(text=a))
#  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 0, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
#  print(a)
#  eval(parse(text=a))
#}



#for(i in grep("myCAGEset",objects(),value = T)) {
#  a <- paste('tc_',strsplit(i,split = "set")[[1]][2],' <- tagClusters(',i,')[[1]]',sep = "")
#  print(a)
#  eval(parse(text=a))
#}



#for(i in grep("tc",objects(),value = TRUE)) {
#  a <- paste('a',paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tss',sep = ""),' <- data.frame(',i,'[,2],',i,'[,8]-1,',i,'[,c(8,6,6,5)],stringsAsFactors = FALSE)',sep = "")
#  eval(parse(text=a))
#  print(a)
#}

#for(i in grep("dominant",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
#  eval(parse(text=a))
#  print(a)
#}



#yhrun -n 1 -p localdisk bedtools intersect -s -a a3_4_5cap_dominant_tss.bed -b gencode_mm10_all_gene_ustream2k_and_genebody.bed -wa -wb > a3_4_5cap_dominant_tss_ustream2k_and_genebody.bed
#yhrun -n 1 -p localdisk bedtools intersect -s -a a3_6_5cap_dominant_tss.bed -b gencode_mm10_all_gene_ustream2k_and_genebody.bed -wa -wb > a3_6_5cap_dominant_tss_ustream2k_and_genebody.bed
#yhrun -n 1 -p localdisk bedtools intersect -s -a a3_15_5cap_dominant_tss.bed -b gencode_mm10_all_gene_ustream2k_and_genebody.bed -wa -wb > a3_15_5cap_dominant_tss_ustream2k_and_genebody.bed


#yhrun -n 1 -p localdisk bedtools intersect -s -a a3_4_3tail_dominant_tss.bed -b gencode_mm10_all_gene_genebody_and_downstream2k.bed -wa -wb > a3_4_3tail_dominant_tss_genebody_and_downstream2k.bed
#yhrun -n 1 -p localdisk bedtools intersect -s -a a3_6_3tail_dominant_tss.bed -b gencode_mm10_all_gene_genebody_and_downstream2k.bed -wa -wb > a3_6_3tail_dominant_tss_genebody_and_downstream2k.bed
#yhrun -n 1 -p localdisk bedtools intersect -s -a a3_15_3tail_dominant_tss.bed -b gencode_mm10_all_gene_genebody_and_downstream2k.bed -wa -wb > a3_15_3tail_dominant_tss_genebody_and_downstream2k.bed




ENSM2symbol <- read.table("G:/annotation/gencode.vM18.annotation.sorted.all.gene.gtf",sep = "\t")

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



D3_4_5cap_dominant_tss_in_gene <- read.table("a3_4_5cap_dominant_tss_ustream2k_and_genebody.bed",header = FALSE)
D3_4_3tail_dominant_tss_in_gene <- read.table("a3_4_3tail_dominant_tss_genebody_and_downstream2k.bed",header = FALSE)


D3_4_5cap_dominant_tss_in_gene <- D3_4_5cap_dominant_tss_in_gene[,c(1,2,3,10,5,6)]
D3_4_3tail_dominant_tss_in_gene <- D3_4_3tail_dominant_tss_in_gene[,c(1,2,3,10,5,6)]



for(i in grep("in_gene",objects(),value = T)) {
  a <- paste(i,'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    b[1,5] <- sum(b[,5])
    a <- paste(i,'_major <- rbind(',i,'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}



OC1_1 <- read.table("G:/CAGEr/CAGEr20190821/OC_1_1_OC_1_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
OC1_2 <- read.table("G:/CAGEr/CAGEr20190821/OC_1_2_OC_1_2.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
OC1_3 <- read.table("G:/CAGEr/CAGEr20190821/OC_1_3_OC_1_3.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
DRG_1 <- read.table("G:/CAGEr/CAGEr20190821/DRG_1_DRG_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
DRG_2 <- read.table("G:/CAGEr/CAGEr20190821/DRG_2_DRG_2.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")


counts <- data.frame(gene=OC1_1[,1],oc1_1=OC1_1[,2],oc1_2=OC1_2[,2],oc1_3=OC1_3[,2],drg_1=DRG_1[,2],drg_2=DRG_2[,2])
counts <- counts[-c(54147:54151),]
counts <- counts[!(counts[,2]==0 & counts[,3]==0 & counts[,4]==0 & counts[,5]==0 & counts[,6]==0),]




for(i in 1:nrow(counts)) {
  if(counts[i,1] %in% ENSM2symbol[,1]) counts[i,1] <- ENSM2symbol[ENSM2symbol[,1]==counts[i,1],2]
}

colnames(counts)[1] <- "gene"


smart <- read.table("G:/CAGEr/CAGEr20190616/smart_seq2_single_cell_martix.txt",header = T)

smart <- data.frame(gene=rownames(smart),read_count=smart[,24])
smart[,2] <- smart[,2]/sum(smart[,2])*10^6
smart <- smart[1:54028,]
smart <- smart[!smart[,2]==0,]


tss_genelist <- D3_4_5cap_dominant_tss_in_gene_major[,4]
tes_genelist <- D3_4_3tail_dominant_tss_in_gene_major[,4]
iso_genelist <- counts[!counts[,2]==0,1]

df3 <- rbind(data.frame(V1=log10(smart[smart[,1] %in% tss_genelist,2]),V2="tss"),
             data.frame(V1=log10(smart[smart[,1] %in% tes_genelist,2]),V2="tes"),
             data.frame(V1=log10(smart[smart[,1] %in% iso_genelist,2]),V2="iso"))



my_comparisons <- list(c("iso", "tes"),c("iso", "tss"))
ggviolin(df3, x="V2", y="V1", fill = "V2", palette = "jco",ylab="smart_RPM_logscale",xlab="",legend = "right") +   ###boxplot可改成violin
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=45,hjust = 1,vjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  stat_compare_means(comparisons=my_comparisons,label = "p-value", method = "wilcox.test") 








####fig1 h


DRG_1 <- read.table("G:/CAGEr/CAGEr20190821/DRG_1_DRG_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")
DRG_2 <- read.table("G:/CAGEr/CAGEr20190821/DRG_2_DRG_2.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")


counts <- data.frame(gene=DRG_1[,1],drg_1=DRG_1[,2],drg_2=DRG_2[,2])
counts <- counts[-c(54147:54151),]
counts <- counts[!(counts[,2]==0 & counts[,3]==0),]


ENSM2symbol <- read.table("G:/annotation/gencode.vM18.annotation.sorted.all.gene.gtf",sep = "\t")

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


for(i in 1:nrow(counts)) {
  if(counts[i,1] %in% ENSM2symbol[,1]) counts[i,1] <- ENSM2symbol[ENSM2symbol[,1]==counts[i,1],2]
}

colnames(counts)[1] <- "gene"



#tes <- read.table("D47_L2_804D74.R1.fastq-common.out_withA10_remainA5_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.bed", header = FALSE, sep = "\t")
#tes_plus <- tes[tes[,6]=="+",]
#tes_minus <- tes[tes[,6]=="-",]
#tes_plus <- tes_plus[,c(1,3,6)]
#tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
#tes <- rbind(tes_plus,tes_minus)
#tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
#write.table(tes,file = "tes.bed",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
#myCAGEsetD47_L2_tes <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "tes.bed",inputFilesType = "bed",sampleLabels = "sample")
#getCTSS(myCAGEsetD47_L2_tes,removeFirstG = FALSE,correctSystematicG = FALSE)


#myCAGEsetD47_L2_tss <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = "D47_L2_804D74.R2.fastq_with_tag.fq.trimed.remainGGG.fq.trimed.nextera_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.bed",inputFilesType = "bed",sampleLabels = "sample")
#getCTSS(myCAGEsetD47_L2_tss,removeFirstG = FALSE,correctSystematicG = FALSE)





#for(i in grep("myCAGEset",objects(),value = T)) {
#  a <- paste('normalizeTagCount(',i,', method = "simpleTpm")',sep = "")
#  print(a)
#  eval(parse(text=a))
#  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 0, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
#  print(a)
#  eval(parse(text=a))
#}



#for(i in grep("myCAGEset",objects(),value = T)) {
#  a <- paste('tc_',strsplit(i,split = "set")[[1]][2],' <- tagClusters(',i,')[[1]]',sep = "")
#  print(a)
#  eval(parse(text=a))
#}



#for(i in grep("tc",objects(),value = TRUE)) {
#  a <- paste('a',paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tss',sep = ""),' <- data.frame(',i,'[,2],',i,'[,8]-1,',i,'[,c(8,6,6,5)],stringsAsFactors = FALSE)',sep = "")
#  eval(parse(text=a))
#  print(a)
#}

#for(i in grep("dominant",objects(),value = TRUE)) {
#  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
#  eval(parse(text=a))
#  print(a)
#}



#yhrun -n 1 -p localdisk bedtools intersect -s -a aD47_L2_tss_dominant_tss.bed -b gencode_mm10_all_gene_ustream2k_and_genebody.bed -wa -wb > aD47_L2_tss_dominant_tss_ustream2k_and_genebody.bed


#yhrun -n 1 -p localdisk bedtools intersect -s -a aD47_L2_tes_dominant_tss.bed -b gencode_mm10_all_gene_genebody_and_downstream2k.bed -wa -wb > aD47_L2_tes_dominant_tss_genebody_and_downstream2k.bed



D47_L2_tss_dominant_tss_in_gene <- read.table("aD47_L2_tss_dominant_tss_ustream2k_and_genebody.bed",header = FALSE)

D47_L2_tes_dominant_tss_in_gene <- read.table("aD47_L2_tes_dominant_tss_genebody_and_downstream2k.bed",header = FALSE)


D47_L2_tss_dominant_tss_in_gene <- D47_L2_tss_dominant_tss_in_gene[,c(1,2,3,10,5,6)]

D47_L2_tes_dominant_tss_in_gene <- D47_L2_tes_dominant_tss_in_gene[,c(1,2,3,10,5,6)]




for(i in grep("in_gene",objects(),value = T)) {
  a <- paste(i,'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    b[1,5] <- sum(b[,5])
    a <- paste(i,'_major <- rbind(',i,'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}


smart <- read.table("G:/CAGEr/CAGEr20190616/smart_seq2_single_cell_martix.txt",header = T)

smart <- data.frame(gene=rownames(smart),read_count=smart[,8]) #if drg

#smart <- data.frame(gene=rownames(smart),read_count=smart[,19]) #if oocyte
#smart[,2] <- smart[,2]/sum(smart[,2])*10^6
#smart[,3] <- smart[,3]/sum(smart[,3])*10^6
#smart[,4] <- smart[,4]/sum(smart[,4])*10^6
smart <- smart[1:54028,]
smart <- smart[!smart[,2]==0,]
smart <- smart[!smart[,1] %in% c("__no_feature","__ambiguous"),]

D49_tss <- read.table("D49_L2_804D76.R2.fastq_with_tag.fq.trimed.remainGGG.fq.trimed.nextera_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.count")
D49_tss <- D49_tss[!D49_tss[,2]==0,]
D49_tss <- D49_tss[!D49_tss[,1] %in% c("__no_feature","__ambiguous"),]

P47_tss <- read.table("45_L1_804D74.R2.fastq_with_tag.fq.trimed.remainGGG.fq.trimed.nextera_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.count")
P47_tss <- P47_tss[!P47_tss[,2]==0,]
P47_tss <- P47_tss[!P47_tss[,1] %in% c("__no_feature","__ambiguous"),]


D49_tes <- read.table("D49_L2_804D76.R1.fastq-common.out_withA10_remainA5_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.count")
D49_tes <- D49_tes[!D49_tes[,2]==0,]
D49_tes <- D49_tes[!D49_tes[,1] %in% c("__no_feature","__ambiguous"),]

P47_tes <- read.table("45_L1_804D74.R1.fastq-common.out_withA10_remainA5_Aligned.out.sam_extract_uniquely_map.sam_extractmismatch_add_header_sorted.count")
P47_tes <- P47_tes[!P47_tes[,2]==0,]
P47_tes <- P47_tes[!P47_tes[,1] %in% c("__no_feature","__ambiguous"),]

counts <- read.table("OC_1_1_OC_1_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")

counts <- counts[!counts[,2]==0,]

counts <- read.table("DRG_1_DRG_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")

counts <- counts[!counts[,2]==0,]


for(i in 1:nrow(counts)) {
  if(counts[i,1] %in% ENSM2symbol[,1]) counts[i,1] <- ENSM2symbol[ENSM2symbol[,1]==counts[i,1],2]
}

for(i in 1:nrow(D49_tss)) {
  if(D49_tss[i,1] %in% ENSM2symbol[,1]) D49_tss[i,1] <- ENSM2symbol[ENSM2symbol[,1]==D49_tss[i,1],2]
}

for(i in 1:nrow(D49_tes)) {
  if(D49_tes[i,1] %in% ENSM2symbol[,1]) D49_tes[i,1] <- ENSM2symbol[ENSM2symbol[,1]==D49_tes[i,1],2]
}



df_temp <- rbind(data.frame(method="ISO-seq",gene=counts[,1],RPM=log10(counts[,2][!counts[,2]==0])),
                 data.frame(method="Smart-seq2",gene=smart[,1],RPM=log10(smart[,2])),
                 data.frame(method="scCAT-seq_TSS",gene=D49_tss[,1],RPM=log10(D49_tss[,2])),
                 data.frame(method="scCAT-seq_TES",gene=D49_tes[,1],RPM=log10(D49_tes[,2])))

df_temp <- df_temp[df_temp[,2] %in% as.character(as.data.frame(table(df_temp[,2]))[as.data.frame(table(df_temp[,2]))[,2]==4,1]),]



ggplot(df_temp,aes(x=RPM))+
  
  geom_density(aes(x=RPM,y=..scaled..,color=method),alpha = 0.5,size=1,linetype="solid") +
  
  labs(x="log10 (read count)")


















######读取机器学习后的数据 

D_dominant_tss_in_gene_1 <- cbind(read.csv("C:/Users/zhong/Desktop/201907211815/20190711DRGtss_peak_new.csv"),
                                  read.csv("C:/Users/zhong/Desktop/201907211815/majority_vote_tss_test_DRG.csv")[,c(12,13)])
D_dominant_tss_in_gene_1 <- D_dominant_tss_in_gene_1[D_dominant_tss_in_gene_1[,17]==1,]
D_dominant_tss_in_gene_1 <- data.frame(V1=D_dominant_tss_in_gene_1[,2],V2=D_dominant_tss_in_gene_1[,7]-1,V3=D_dominant_tss_in_gene_1[,7],
                                       V10=D_dominant_tss_in_gene_1[,1],V5=D_dominant_tss_in_gene_1[,6],V6=D_dominant_tss_in_gene_1[,5])




D_dominant_tes_in_gene_1 <- cbind(read.csv("C:/Users/zhong/Desktop/201907211815/20190711DRGtes_peak_new.csv"),
                                  read.csv("C:/Users/zhong/Desktop/201907211815/majority_vote_tes_test_DRG.csv")[,c(12,13)])
D_dominant_tes_in_gene_1 <- D_dominant_tes_in_gene_1[D_dominant_tes_in_gene_1[,17]==1,]
D_dominant_tes_in_gene_1 <- data.frame(V1=D_dominant_tes_in_gene_1[,2],V2=D_dominant_tes_in_gene_1[,7]-1,V3=D_dominant_tes_in_gene_1[,7],
                                       V10=D_dominant_tes_in_gene_1[,1],V5=D_dominant_tes_in_gene_1[,6],V6=D_dominant_tes_in_gene_1[,5])






counts <- read.table("DRG_1_DRG_1.ccs.trim_5primer.trim_3primer.trim_polyA.sorted.count")

counts <- counts[!counts[,2]==0,]

for(i in 1:nrow(counts)) {
  if(counts[i,1] %in% ENSM2symbol[,1]) counts[i,1] <- ENSM2symbol[ENSM2symbol[,1]==counts[i,1],2]
}





df_temp <- data.frame()
for(i in unique(counts[,2])) {
  if(TRUE %in% (D47_L2_tss_dominant_tss_in_gene_major[,4] %in% counts[counts[,2]==i,1])) {
    df_temp <- rbind(df_temp,
                     data.frame(ccs=i,RPM=log10(D47_L2_tss_dominant_tss_in_gene_major[D47_L2_tss_dominant_tss_in_gene_major[,4] %in% counts[counts[,2]==i,1],5])))
    
  }
}



df_temp[,1][df_temp[,1]>10] <- ">10"
df_temp <- df_temp[!df_temp[,1]=="0",]
df_temp[,1] <- factor(df_temp[,1],levels = c("1","2","3","4","5","6","7","8","9","10",">10"))

ggboxplot(df_temp, x="ccs", y="RPM", fill = "ccs", palette = as.character(jdb_palette("brewer_blue",type="continuous"))[c(1:11)*90],ylab="scCAT-seq TSS RPM_logscale",xlab="ISO-seq ccs number",legend = "right",add = "jitter",add.params = list(binwidth = 0.02)) + 
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))+
  scale_y_continuous(limits = c(0,5.5))
