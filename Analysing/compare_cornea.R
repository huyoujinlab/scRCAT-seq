options(stringsAsFactors = FALSE)
options(scipen = 100)

setwd("G:/CAGEr/CAGEr20190719LEC_APA/")
library(CAGEr)
library(rstatix)


##--------------形成CAGEset
for(i in grep("CTRL",list.files(),value = T)) {
  tes <- read.table(i, header = FALSE, sep = "\t")
  tes_plus <- tes[tes[,6]=="+",]
  tes_minus <- tes[tes[,6]=="-",]
  tes_plus <- tes_plus[,c(1,3,6)]
  tes_minus <- data.frame(V1=tes_minus[,1],V3=tes_minus[,2]+1,V6=tes_minus[,6])
  tes <- rbind(tes_plus,tes_minus)
  tes <- data.frame(tes[,1],tes[,2]-1,tes[,2],rep(1,nrow(tes)),rep(1,nrow(tes)),tes[,3])
  write.table(tes,file = "tes.bed",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
  a <- paste('myCAGEset',strsplit(i,split = "[.]")[[1]][1],'_3tail <- new("CAGEset", genomeName = "BSgenome.Macfas.ZJW.1",inputFiles = "tes.bed",inputFilesType = "bed",sampleLabels = "',strsplit(i,split = "[.]")[[1]][1],'")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('getCTSS(','myCAGEset',strsplit(i,split = "[.]")[[1]][1],'_3tail,removeFirstG = FALSE,correctSystematicG = FALSE)',sep = "")
  print(a)
  eval(parse(text=a))
}



##########call peak

for(i in grep("myCAGEset",objects(),value = T)) {
  a <- paste('normalizeTagCount(',i,', method = "simpleTpm")',sep = "")
  print(a)
  eval(parse(text=a))
  a <- paste('clusterCTSS(object = ',i,', method = "distclu", threshold = 5, nrPassThreshold = 1, thresholdIsTpm = TRUE,removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = FALSE, nrCores = 8)',sep = "")
  print(a)
  eval(parse(text=a))
}



for(i in grep("myCAGEset",objects(),value = T)) {
  a <- paste('tc_',strsplit(i,split = "set")[[1]][2],' <- tagClusters(',i,')[[1]]',sep = "")
  print(a)
  eval(parse(text=a))
}

for(i in grep("tc",objects(),value = TRUE)) {
  a <- paste(paste(strsplit(i,split = "c_")[[1]][2],'_dominant_tes',sep = ""),' <- data.frame(',i,'[,2],',i,'[,8]-1,',i,'[,c(8,6,6,5)],stringsAsFactors = FALSE)',sep = "")
  eval(parse(text=a))
  print(a)
}

for(i in grep("dominant",objects(),value = TRUE)) {
  a <- paste('write.table(',i,',"',paste(i,'.bed',sep = ""),'",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")',sep = "")
  eval(parse(text=a))
}




#ref <- read.table("G:/revref_Macfas_sorted.gtf",sep = "\t")
#ref <- ref[ref[,3]=="transcript",]


#gene_id <- strsplit(ref[,9],split = "gene_id ")
#gene_id <- lapply(gene_id, function(x) x[2])
#gene_id <- unlist(gene_id)
#gene_id <- strsplit(gene_id,split = "; ")
#gene_id <- lapply(gene_id, function(x) x[1])
#gene_id <- unlist(gene_id)


#ref_bed <- data.frame(chr=ref[,1],start=ref[,4]-1,end=ref[,5],gene_id=gene_id,V5=gene_id,V6=ref[,7])
#ref_bed <- ref_bed[!ref_bed[,1] %in% grep("random",ref_bed[,1],value = T),]


#ref_bed_new <- data.frame()
#for(i in unique(ref_bed[,4])) {
#  temp <- ref_bed[ref_bed[,4]==i,]
#  temp[,7] <- temp[,3]-temp[,2]
#  temp <- temp[order(temp[,7],decreasing = T),]
#  ref_bed_new <- rbind(ref_bed_new,temp[1,1:6])
#}
#write.table(ref_bed_new,"revref_Macfas_gene.bed",quote = F,row.names = FALSE,col.names = F,sep = "\t")
#ref_bed_new <- read.table("revref_Macfas_gene.bed",sep = "\t")
#colnames(ref_bed_new) <- c("chr","start","end","gene_id","V5","strand")
#ref_bed_downstream2k <- data.frame()
#for(i in 1:nrow(ref_bed_new)) {
#  if(ref_bed_new[i,6]=="+") {
#    ref_bed_downstream2k <- rbind(ref_bed_downstream2k,
#                                  data.frame(ref_bed_new[i,c(1,2)],end=ref_bed_new[i,3]+2000,ref_bed_new[i,c(4,5,6)]))
#  } else {
#    ref_bed_downstream2k <- rbind(ref_bed_downstream2k,
#                                  data.frame(chr=ref_bed_new[i,1],start=ref_bed_new[i,2]-2000,ref_bed_new[i,c(3,4,5,6)]))
#  }
#  print(i)
#}
#ref_bed_downstream2k[,2][ref_bed_downstream2k[,2]<0] <- 0
#write.table(ref_bed_downstream2k,"revref_Macfas_gene_downstream2k.bed",quote = F,row.names = FALSE,col.names = F,sep = "\t")





for i in `ls|grep "tes.bed"`; do bedtools intersect -s -a ${i} -b revref_Macfas_gene_downstream2k.bed -wa -wb > ${i%%.*}_genebody_downstream2k.bed; done


rm(list = ls())



##### CEC and TAC
CTRL_0_3tail_dominant_tes_in_gene <- read.table("CTRL_0_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)
CTRL_2_3tail_dominant_tes_in_gene <- read.table("CTRL_2_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)


CTRL_0_3tail_dominant_tes_in_gene <- CTRL_0_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]
CTRL_2_3tail_dominant_tes_in_gene <- CTRL_2_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]



#取最高峰
for(i in grep("_3tail_dominant_tes_in_gene",objects(),value = T)) {
  a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- rbind(',strsplit(i,split = "_3tail")[[1]][1],'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}




CTRL_0_major[,7] <- rep("CTRL_0",nrow(CTRL_0_major))
CTRL_2_major[,7] <- rep("CTRL_2",nrow(CTRL_2_major))








############联合 CRTL0和CTRL1         从这里到下面的，把它们复制到一个新的是script里，然后改CTRL_0,CTRL_1等
# 1 2 LEC TAC
# 0 2 CEC TAC
union <- rbind(CTRL_2_major,CTRL_0_major)

union_2_diff <- data.frame()
union_2_same <- data.frame()
union_1 <- data.frame()

for(i in unique(union[,4])) {
  a <- union[union[,4] %in% i,]
  if(nrow(a)==2) {
    if(abs(a[1,3]-a[2,3])<21) union_2_same <- rbind(union_2_same,a)
    else union_2_diff <- rbind(union_2_diff,a)
  }
  if(nrow(a)==1) union_1 <- rbind(union_1,a)
}


union_2_diff_true <- data.frame()
for(i in 1:(nrow(union_2_diff)/2)) {
  one <- union_2_diff[2*i-1,]
  three <- union_2_diff[2*i,]
  two_tes <- CTRL_0_3tail_dominant_tes_in_gene[CTRL_0_3tail_dominant_tes_in_gene[,4]==one[1,4],]
  two_tes <- two_tes[order(two_tes[,5],decreasing = TRUE),]
  d <- two_tes[abs(two_tes[,3]-one[1,3])<21,]
  two <- data.frame(V1=one[1,1],V2=d[1,2],V3=d[1,3],V10=d[1,4],
                    V5=d[1,5],V6=one[1,6],V7="CTRL_0")
  four_tes <- CTRL_2_3tail_dominant_tes_in_gene[CTRL_2_3tail_dominant_tes_in_gene[,4]==three[1,4],]
  four_tes <- four_tes[order(four_tes[,5],decreasing = TRUE),]
  d <- four_tes[abs(four_tes[,3]-three[1,3])<21,]
  four <- data.frame(V1=three[1,1],V2=d[1,2],V3=d[1,3],V10=d[1,4],
                     V5=d[1,5],V6=three[1,6],V7="CTRL_2")
  
  union_2_diff_true <- rbind(union_2_diff_true,one,two,three,four)
}





for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(is.na(union_2_diff_true[4*i-2,3])) {
    union_2_diff_true[4*i-2,2] <- union_2_diff_true[4*i-3,2]
    union_2_diff_true[4*i-2,3] <- union_2_diff_true[4*i-3,3]
    union_2_diff_true[4*i-2,4] <- union_2_diff_true[4*i-3,4]
    union_2_diff_true[4*i-2,5] <- 1
  }
  if(is.na(union_2_diff_true[4*i,3])) {
    union_2_diff_true[4*i,2] <- union_2_diff_true[4*i-1,2]
    union_2_diff_true[4*i,3] <- union_2_diff_true[4*i-1,3]
    union_2_diff_true[4*i,4] <- union_2_diff_true[4*i-1,4]
    union_2_diff_true[4*i,5] <- 1
  }
}




for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(union_2_diff_true[4*i-3,6]=="+") {
    if(union_2_diff_true[4*i-3,3]>union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "long"
      union_2_diff_true[4*i-2,8] <- "long"
      union_2_diff_true[4*i-1,8] <- "short"
      union_2_diff_true[4*i,8] <- "short"
    }
    if(union_2_diff_true[4*i-3,3]<union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "short"
      union_2_diff_true[4*i-2,8] <- "short"
      union_2_diff_true[4*i-1,8] <- "long"
      union_2_diff_true[4*i,8] <- "long"
    }
  }
  if(union_2_diff_true[4*i-3,6]=="-") {
    if(union_2_diff_true[4*i-3,3]>union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "short"
      union_2_diff_true[4*i-2,8] <- "short"
      union_2_diff_true[4*i-1,8] <- "long"
      union_2_diff_true[4*i,8] <- "long"
    }
    if(union_2_diff_true[4*i-3,3]<union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "long"
      union_2_diff_true[4*i-2,8] <- "long"
      union_2_diff_true[4*i-1,8] <- "short"
      union_2_diff_true[4*i,8] <- "short"
    }
  }
}









for(i in 1:(nrow(union_2_diff_true)/4)) {
  a <- union_2_diff_true[(4*i-3):(4*i),]
  union_2_diff_true[(4*i-3):(4*i),9] <- fisher.test(matrix(c(ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5]),
                                                             ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5]),
                                                             ceiling(a[a[,7]=="CTRL_0" & a[,8]=="long",5]),
                                                             ceiling(a[a[,7]=="CTRL_0" & a[,8]=="short",5])),byrow = T,nrow = 2))$p.value
  if(ceiling(a[a[,7]=="CTRL_0" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_0" & a[,8]=="short",5])>ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5])) {
    union_2_diff_true[(4*i-3):(4*i),10] <- "longer"
  }
  if(ceiling(a[a[,7]=="CTRL_0" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_0" & a[,8]=="short",5])<ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5])) {
    union_2_diff_true[(4*i-3):(4*i),10] <- "shorter"
  }
  union_2_diff_true[(4*i-3),11] <- log10(ceiling(a[a[,7]=="CTRL_0" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5]))
  union_2_diff_true[(4*i-3),12] <- log10(ceiling(a[a[,7]=="CTRL_0" & a[,8]=="short",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5]))
}



for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(union_2_diff_true[4*i-3,9]<0.05) {
    if(union_2_diff_true[4*i-3,10]=="longer") {
      union_2_diff_true[4*i-3,13] <- "Longer"
    }
    else union_2_diff_true[4*i-3,13] <- "Shorter"
  }
  else union_2_diff_true[4*i-3,13] <- "No change"
}



colnames(union_2_diff_true)[9:13] <- c("TES_pvalue","TES_choose_longer?","TES_long_iso_later_devide_former","TES_short_iso_later_devide_former","TES_event")

a <- na.omit(union_2_diff_true)


nrow(a[a[,13]=="Longer",])
nrow(a[a[,13]=="Shorter",])
nrow(a[a[,13]=="No change",])
a[,13] <- factor(a[,13],levels = c("No change","Longer","Shorter"))

ggplot(a,aes(x=TES_long_iso_later_devide_former,y=TES_short_iso_later_devide_former,color=TES_event))+
  geom_point()+
  labs(x="Fold change of abundance of isoforms with distal TES (log10)",
       y="Fold change of abundance of isoforms with proximal TESs (log10)") +
  scale_color_manual(values = c("#BEBEBE","#00CDCD","#FF4500"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1,linetype = 1), #改外边框
        axis.text.x = element_text(face = "plain",size = 10.5), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 10.5)) +#改y轴字体 
  geom_abline(slope = 1, intercept=0, na.rm = FALSE, show.legend = NA,linetype="dashed",size=1)+
  theme(legend.position=c(10,10))+   ##去除图例
  xlim(c(-3.5,3.5))+
  ylim(c(-3.5,3.5))




CTRL_0_3tail_dominant_tes_in_gene <- read.table("CTRL_0_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)
CTRL_2_3tail_dominant_tes_in_gene <- read.table("CTRL_2_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)


CTRL_0_3tail_dominant_tes_in_gene <- CTRL_0_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]
CTRL_2_3tail_dominant_tes_in_gene <- CTRL_2_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]



#取总和峰
for(i in grep("_3tail_dominant_tes_in_gene",objects(),value = T)) {
  a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    b[1,5] <- sum(b[,5])
    a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- rbind(',strsplit(i,split = "_3tail")[[1]][1],'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}


a <- na.omit(union_2_diff_true)

longer_gene_RPM_CTRL_2 <- CTRL_2_major[CTRL_2_major[,4] %in% a[a[,13]=="Longer",4],]

longer_gene_RPM_CTRL_0 <- CTRL_0_major[CTRL_0_major[,4] %in% a[a[,13]=="Longer",4],]

temp_1 <- rbind(data.frame(longer_gene_RPM_CTRL_2[,c(4,5)],group="TAC"),   #1
                data.frame(longer_gene_RPM_CTRL_0[,c(4,5)],group="CEC"))   #1

temp_1[,2] <- log10(temp_1[,2]+1)  #1



shorter_gene_RPM_CTRL_2 <- CTRL_2_major[CTRL_2_major[,4] %in% a[a[,13]=="Shorter",4],]   

shorter_gene_RPM_CTRL_0 <- CTRL_0_major[CTRL_0_major[,4] %in% a[a[,13]=="Shorter",4],]   

temp_2 <- rbind(data.frame(shorter_gene_RPM_CTRL_2[,c(4,5)],group="TAC"),   #1
                data.frame(shorter_gene_RPM_CTRL_0[,c(4,5)],group="CEC"))   #1

temp_2[,2] <- log10(temp_2[,2]+1)   #1






my_comparisons <- list(c("TAC", "CEC"))


ggboxplot(temp_1,x="group",y="V5",fill = "group",palette = "jco",  line.size = 0) +
  stat_compare_means(comparisons = my_comparisons,paired = T, label = "p-value",label.y = 3.5)+
  scale_y_continuous(breaks = seq(0,3.5,0.5))

ggboxplot(temp_2,x="group",y="V5",fill = "group",palette = "jco",  line.size = 0) +
  stat_compare_means(comparisons = my_comparisons,paired = T, label = "p-value",label.y = 3.5)+
  scale_y_continuous(breaks = seq(0,3.5,0.5))








##### LEC and TAC

CTRL_1_3tail_dominant_tes_in_gene <- read.table("CTRL_1_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)
CTRL_2_3tail_dominant_tes_in_gene <- read.table("CTRL_2_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)


CTRL_1_3tail_dominant_tes_in_gene <- CTRL_1_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]
CTRL_2_3tail_dominant_tes_in_gene <- CTRL_2_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]



#取最高峰
for(i in grep("_3tail_dominant_tes_in_gene",objects(),value = T)) {
  a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- rbind(',strsplit(i,split = "_3tail")[[1]][1],'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}




CTRL_1_major[,7] <- rep("CTRL_1",nrow(CTRL_1_major))
CTRL_2_major[,7] <- rep("CTRL_2",nrow(CTRL_2_major))






union <- rbind(CTRL_2_major,CTRL_1_major)

union_2_diff <- data.frame()
union_2_same <- data.frame()
union_1 <- data.frame()

for(i in unique(union[,4])) {
  a <- union[union[,4] %in% i,]
  if(nrow(a)==2) {
    if(abs(a[1,3]-a[2,3])<21) union_2_same <- rbind(union_2_same,a)
    else union_2_diff <- rbind(union_2_diff,a)
  }
  if(nrow(a)==1) union_1 <- rbind(union_1,a)
}


union_2_diff_true <- data.frame()
for(i in 1:(nrow(union_2_diff)/2)) {
  one <- union_2_diff[2*i-1,]
  three <- union_2_diff[2*i,]
  two_tes <- CTRL_1_3tail_dominant_tes_in_gene[CTRL_1_3tail_dominant_tes_in_gene[,4]==one[1,4],]
  two_tes <- two_tes[order(two_tes[,5],decreasing = TRUE),]
  d <- two_tes[abs(two_tes[,3]-one[1,3])<21,]
  two <- data.frame(V1=one[1,1],V2=d[1,2],V3=d[1,3],V10=d[1,4],
                    V5=d[1,5],V6=one[1,6],V7="CTRL_1")
  four_tes <- CTRL_2_3tail_dominant_tes_in_gene[CTRL_2_3tail_dominant_tes_in_gene[,4]==three[1,4],]
  four_tes <- four_tes[order(four_tes[,5],decreasing = TRUE),]
  d <- four_tes[abs(four_tes[,3]-three[1,3])<21,]
  four <- data.frame(V1=three[1,1],V2=d[1,2],V3=d[1,3],V10=d[1,4],
                     V5=d[1,5],V6=three[1,6],V7="CTRL_2")
  
  union_2_diff_true <- rbind(union_2_diff_true,one,two,three,four)
}





for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(is.na(union_2_diff_true[4*i-2,3])) {
    union_2_diff_true[4*i-2,2] <- union_2_diff_true[4*i-3,2]
    union_2_diff_true[4*i-2,3] <- union_2_diff_true[4*i-3,3]
    union_2_diff_true[4*i-2,4] <- union_2_diff_true[4*i-3,4]
    union_2_diff_true[4*i-2,5] <- 1
  }
  if(is.na(union_2_diff_true[4*i,3])) {
    union_2_diff_true[4*i,2] <- union_2_diff_true[4*i-1,2]
    union_2_diff_true[4*i,3] <- union_2_diff_true[4*i-1,3]
    union_2_diff_true[4*i,4] <- union_2_diff_true[4*i-1,4]
    union_2_diff_true[4*i,5] <- 1
  }
}




for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(union_2_diff_true[4*i-3,6]=="+") {
    if(union_2_diff_true[4*i-3,3]>union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "long"
      union_2_diff_true[4*i-2,8] <- "long"
      union_2_diff_true[4*i-1,8] <- "short"
      union_2_diff_true[4*i,8] <- "short"
    }
    if(union_2_diff_true[4*i-3,3]<union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "short"
      union_2_diff_true[4*i-2,8] <- "short"
      union_2_diff_true[4*i-1,8] <- "long"
      union_2_diff_true[4*i,8] <- "long"
    }
  }
  if(union_2_diff_true[4*i-3,6]=="-") {
    if(union_2_diff_true[4*i-3,3]>union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "short"
      union_2_diff_true[4*i-2,8] <- "short"
      union_2_diff_true[4*i-1,8] <- "long"
      union_2_diff_true[4*i,8] <- "long"
    }
    if(union_2_diff_true[4*i-3,3]<union_2_diff_true[4*i-1,3]) {
      union_2_diff_true[4*i-3,8] <- "long"
      union_2_diff_true[4*i-2,8] <- "long"
      union_2_diff_true[4*i-1,8] <- "short"
      union_2_diff_true[4*i,8] <- "short"
    }
  }
}









for(i in 1:(nrow(union_2_diff_true)/4)) {
  a <- union_2_diff_true[(4*i-3):(4*i),]
  union_2_diff_true[(4*i-3):(4*i),9] <- fisher.test(matrix(c(ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5]),
                                                             ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5]),
                                                             ceiling(a[a[,7]=="CTRL_1" & a[,8]=="long",5]),
                                                             ceiling(a[a[,7]=="CTRL_1" & a[,8]=="short",5])),byrow = T,nrow = 2))$p.value
  if(ceiling(a[a[,7]=="CTRL_1" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_1" & a[,8]=="short",5])>ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5])) {
    union_2_diff_true[(4*i-3):(4*i),10] <- "longer"
  }
  if(ceiling(a[a[,7]=="CTRL_1" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_1" & a[,8]=="short",5])<ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5])) {
    union_2_diff_true[(4*i-3):(4*i),10] <- "shorter"
  }
  union_2_diff_true[(4*i-3),11] <- log10(ceiling(a[a[,7]=="CTRL_1" & a[,8]=="long",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="long",5]))
  union_2_diff_true[(4*i-3),12] <- log10(ceiling(a[a[,7]=="CTRL_1" & a[,8]=="short",5])/ceiling(a[a[,7]=="CTRL_2" & a[,8]=="short",5]))
}



for(i in 1:(nrow(union_2_diff_true)/4)) {
  if(union_2_diff_true[4*i-3,9]<0.05) {
    if(union_2_diff_true[4*i-3,10]=="longer") {
      union_2_diff_true[4*i-3,13] <- "Longer"
    }
    else union_2_diff_true[4*i-3,13] <- "Shorter"
  }
  else union_2_diff_true[4*i-3,13] <- "No change"
}



colnames(union_2_diff_true)[9:13] <- c("TES_pvalue","TES_choose_longer?","TES_long_iso_later_devide_former","TES_short_iso_later_devide_former","TES_event")

a <- na.omit(union_2_diff_true)


nrow(a[a[,13]=="Longer",])
nrow(a[a[,13]=="Shorter",])
nrow(a[a[,13]=="No change",])
a[,13] <- factor(a[,13],levels = c("No change","Longer","Shorter"))

ggplot(a,aes(x=TES_long_iso_later_devide_former,y=TES_short_iso_later_devide_former,color=TES_event))+
  geom_point()+
  labs(x="Fold change of abundance of isoforms with distal TES (log10)",
       y="Fold change of abundance of isoforms with proximal TESs (log10)") +
  scale_color_manual(values = c("#BEBEBE","#00CDCD","#FF4500"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1,linetype = 1), #改外边框
        axis.text.x = element_text(face = "plain",size = 10.5), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 10.5)) +#改y轴字体 
  geom_abline(slope = 1, intercept=0, na.rm = FALSE, show.legend = NA,linetype="dashed",size=1)+
  theme(legend.position=c(10,10))+   ##去除图例
  xlim(c(-3.5,3.5))+
  ylim(c(-3.5,3.5))





























CTRL_1_3tail_dominant_tes_in_gene <- read.table("CTRL_1_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)
CTRL_2_3tail_dominant_tes_in_gene <- read.table("CTRL_2_3tail_dominant_tes_genebody_downstream2k.bed",header = FALSE)


CTRL_1_3tail_dominant_tes_in_gene <- CTRL_1_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]
CTRL_2_3tail_dominant_tes_in_gene <- CTRL_2_3tail_dominant_tes_in_gene[,c(1,2,3,10,5,6)]



#取总和峰
for(i in grep("_3tail_dominant_tes_in_gene",objects(),value = T)) {
  a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- data.frame()',sep = "")
  eval(parse(text=a))
  print(a)
  a <- paste('temp <- ',i,sep = "")
  eval(parse(text=a))
  print(a)
  for(j in unique(temp[,4])) {
    b <- temp[temp[,4]==j,]
    b <- b[order(b[,5],decreasing = T),]
    b[1,5] <- sum(b[,5])
    a <- paste(strsplit(i,split = "_3tail")[[1]][1],'_major <- rbind(',strsplit(i,split = "_3tail")[[1]][1],'_major,b[1,])',sep = "")
    eval(parse(text=a))
  }
}


a <- na.omit(union_2_diff_true)

longer_gene_RPM_CTRL_2 <- CTRL_2_major[CTRL_2_major[,4] %in% a[a[,13]=="Longer",4],]

longer_gene_RPM_CTRL_1 <- CTRL_1_major[CTRL_1_major[,4] %in% a[a[,13]=="Longer",4],]

temp_1 <- rbind(data.frame(longer_gene_RPM_CTRL_2[,c(4,5)],group="TAC"),   #1
                data.frame(longer_gene_RPM_CTRL_1[,c(4,5)],group="LEC"))   #1

temp_1[,2] <- log10(temp_1[,2]+1)  #1



shorter_gene_RPM_CTRL_2 <- CTRL_2_major[CTRL_2_major[,4] %in% a[a[,13]=="Shorter",4],]   

shorter_gene_RPM_CTRL_1 <- CTRL_1_major[CTRL_1_major[,4] %in% a[a[,13]=="Shorter",4],]   

temp_2 <- rbind(data.frame(shorter_gene_RPM_CTRL_2[,c(4,5)],group="TAC"),   #1
                data.frame(shorter_gene_RPM_CTRL_1[,c(4,5)],group="LEC"))   #1

temp_2[,2] <- log10(temp_2[,2]+1)   #1






my_comparisons <- list(c("TAC", "LEC"))

ggboxplot(temp_1,x="group",y="V5",fill = "group",palette = "jco",  line.size = 0) +
  stat_compare_means(comparisons = my_comparisons,paired = T, label = "p-value",label.y = 3.5)+
  scale_y_continuous(breaks = seq(0,3.5,0.5))

ggboxplot(temp_2,x="group",y="V5",fill = "group",palette = "jco",  line.size = 0) +
  stat_compare_means(comparisons = my_comparisons,paired = T, label = "p-value",label.y = 3.5)+
  scale_y_continuous(breaks = seq(0,3.5,0.5))




