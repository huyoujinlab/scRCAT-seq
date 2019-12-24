options(stringsAsFactors = FALSE)
options(scipen = 100)


library(CAGEr)
library(splines)
library(data.table)
library(ggplot2)
library(ggalt)
library(ggsignif)
library(gridExtra)


setwd("C:/Users/zhong/Desktop/20190712peak/")
a <- read.csv("C:/Users/zhong/Desktop/20190712peak/20190711DRGtss_peak_new.csv")
a <- data.frame(a[,c(1,5,6,8,11)],peak_length=(a[,4]-a[,3]),dominant_tss_RPM.smart2_seq_RPM=a[,8]/a[,10],a[,c(14,15,13,2,3,4)])
a[is.na(a)] <- 0
a[,7][is.infinite(a[,7])] <- 0

logistic_DRG_tss <- read.csv("logistic_tss_test_DRG.csv")
rf_DRG_tss <- read.csv("rf_tss_test_DRG.csv")
svm_DRG_tss <- read.csv("SVM_tss_test_DRG.csv")

a <- cbind(a,logistic_DRG_tss[,12],rf_DRG_tss[,12],svm_DRG_tss[,12])

write.table(a[a[,14]==1,c(11,12,13,1,1,2)],"logistic_DRG_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,15]==1,c(11,12,13,1,1,2)],"rf_DRG_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,16]==1,c(11,12,13,1,1,2)],"svm_DRG_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[,c(11,12,13,1,1,2)],"origin_DRG_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
logistic_DRG_tss <- a[a[,14]==1,c(11,12,13,1,1,2)]
rf_DRG_tss <- a[a[,15]==1,c(11,12,13,1,1,2)]
svm_DRG_tss <- a[a[,16]==1,c(11,12,13,1,1,2)]
origin_DRG_tss <- a[,c(11,12,13,1,1,2)]



a <- read.csv("C:/Users/zhong/Desktop/20190712peak/20190711D3tss_peak_new.csv")
a <- data.frame(a[,c(1,5,6,8,11)],peak_length=(a[,4]-a[,3]),dominant_tss_RPM.smart2_seq_RPM=a[,8]/a[,10],a[,c(14,15,13,2,3,4)])
a[is.na(a)] <- 0
a[,7][is.infinite(a[,7])] <- 0

logistic_D3_tss <- read.csv("logistic_tss_test_D3.csv")
rf_D3_tss <- read.csv("rf_tss_test_D3.csv")
svm_D3_tss <- read.csv("SVM_tss_test_D3.csv")

a <- cbind(a,logistic_D3_tss[,12],rf_D3_tss[,12],svm_D3_tss[,12])

write.table(a[a[,14]==1,c(11,12,13,1,1,2)],"logistic_D3_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,15]==1,c(11,12,13,1,1,2)],"rf_D3_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,16]==1,c(11,12,13,1,1,2)],"svm_D3_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[,c(11,12,13,1,1,2)],"origin_D3_tss.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
logistic_D3_tss <- a[a[,14]==1,c(11,12,13,1,1,2)]
rf_D3_tss <- a[a[,15]==1,c(11,12,13,1,1,2)]
svm_D3_tss <- a[a[,16]==1,c(11,12,13,1,1,2)]
origin_D3_tss <- a[,c(11,12,13,1,1,2)]



a <- read.csv("C:/Users/zhong/Desktop/20190712peak/20190711DRGtes_peak_new.csv")
a <- data.frame(a[,c(1,5,6,8,11)],peak_length=(a[,4]-a[,3]),dominant_tes_RPM.smart2_seq_RPM=a[,8]/a[,10],a[,c(14,15,13,2,3,4)])
a[is.na(a)] <- 0
a[,7][is.infinite(a[,7])] <- 0

logistic_DRG_tes <- read.csv("logistic_tes_test_DRG.csv")
rf_DRG_tes <- read.csv("rf_tes_test_DRG.csv")
svm_DRG_tes <- read.csv("SVM_tes_test_DRG.csv")

a <- cbind(a,logistic_DRG_tes[,12],rf_DRG_tes[,12],svm_DRG_tes[,12])

write.table(a[a[,14]==1,c(11,12,13,1,1,2)],"logistic_DRG_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,15]==1,c(11,12,13,1,1,2)],"rf_DRG_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,16]==1,c(11,12,13,1,1,2)],"svm_DRG_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[,c(11,12,13,1,1,2)],"origin_DRG_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
logistic_DRG_tes <- a[a[,14]==1,c(11,12,13,1,1,2)]
rf_DRG_tes <- a[a[,15]==1,c(11,12,13,1,1,2)]
svm_DRG_tes <- a[a[,16]==1,c(11,12,13,1,1,2)]
origin_DRG_tes <- a[,c(11,12,13,1,1,2)]



a <- read.csv("C:/Users/zhong/Desktop/20190712peak/20190711D3tes_peak_new.csv")
a <- data.frame(a[,c(1,5,6,8,11)],peak_length=(a[,4]-a[,3]),dominant_tes_RPM.smart2_seq_RPM=a[,8]/a[,10],a[,c(14,15,13,2,3,4)])
a[is.na(a)] <- 0
a[,7][is.infinite(a[,7])] <- 0

logistic_D3_tes <- read.csv("logistic_tes_test_D3.csv")
rf_D3_tes <- read.csv("rf_tes_test_D3.csv")
svm_D3_tes <- read.csv("SVM_tes_test_D3.csv")

a <- cbind(a,logistic_D3_tes[,12],rf_D3_tes[,12],svm_D3_tes[,12])

write.table(a[a[,14]==1,c(11,12,13,1,1,2)],"logistic_D3_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,15]==1,c(11,12,13,1,1,2)],"rf_D3_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[a[,16]==1,c(11,12,13,1,1,2)],"svm_D3_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(a[,c(11,12,13,1,1,2)],"origin_D3_tes.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
logistic_D3_tes <- a[a[,14]==1,c(11,12,13,1,1,2)]
rf_D3_tes <- a[a[,15]==1,c(11,12,13,1,1,2)]
svm_D3_tes <- a[a[,16]==1,c(11,12,13,1,1,2)]
origin_D3_tes <- a[,c(11,12,13,1,1,2)]





gencode <- read.table("G:/CAGEr/CAGEr20190318tata_new/gencode_mm10_all_gene_all_transcript.bed",header = FALSE,sep = "\t")
gencode_plus <- gencode[gencode[,6]=="+",]
gencode_minus <- gencode[gencode[,6]=="-",]


tss_refer <- rbind(data.frame(gencode_plus[,c(1,2)],V3=gencode_plus[,2]+1,gencode_plus[,c(4,5,6)]),
                   data.frame(V1=gencode_minus[,1],V2=gencode_minus[,3]-1,gencode_minus[c(3,4,5,6)]))
tes_refer <- rbind(data.frame(V1=gencode_plus[,1],V2=gencode_plus[,3]-1,gencode_plus[,c(3,4,5,6)]),
                   data.frame(gencode_minus[,c(1,2)],V3=gencode_minus[,2]+1,gencode_minus[c(4,5,6)]))




logistic_D3_tss_ref <- c()
for(i in 1:nrow(logistic_D3_tss)) {
  b <- tss_refer[tss_refer[,1]==logistic_D3_tss[i,1] & tss_refer[,6]==logistic_D3_tss[i,6],]
  if(logistic_D3_tss[i,6]=="+") {
    logistic_D3_tss_ref <- c(logistic_D3_tss_ref,(logistic_D3_tss[i,3]-b[,3])[order(abs(logistic_D3_tss[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    logistic_D3_tss_ref <- c(logistic_D3_tss_ref,(b[,3]-logistic_D3_tss[i,3])[order(abs(b[,3]-logistic_D3_tss[i,3]),decreasing = FALSE)][1])
  }
}
#logistic_D3_tss_ref <- logistic_D3_tss_ref[logistic_D3_tss_ref<5001 & logistic_D3_tss_ref>-5001]

tss_df <- data.frame(V1=as.data.frame(table(logistic_D3_tss_ref))[,1],V2=as.data.frame(table(logistic_D3_tss_ref))[,2],V3=rep("test",nrow(as.data.frame(table(logistic_D3_tss_ref)))))

tss_df[,1] <- as.numeric(as.character(tss_df[,1]))
tss_df <- tss_df[tss_df[,1]<201 & tss_df[,1]>-201,] 



logistic_D3_tes_ref <- c()
for(i in 1:nrow(logistic_D3_tes)) {
  b <- tes_refer[tes_refer[,1]==logistic_D3_tes[i,1] & tes_refer[,6]==logistic_D3_tes[i,6],]
  if(logistic_D3_tes[i,6]=="+") {
    logistic_D3_tes_ref <- c(logistic_D3_tes_ref,(logistic_D3_tes[i,3]-b[,3])[order(abs(logistic_D3_tes[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    logistic_D3_tes_ref <- c(logistic_D3_tes_ref,(b[,3]-logistic_D3_tes[i,3])[order(abs(b[,3]-logistic_D3_tes[i,3]),decreasing = FALSE)][1])
  }
}
#logistic_D3_tses_ref <- logistic_D3_tes_ref[logistic_D3_tes_ref<5001 & logistic_D3_tes_ref>-5001]

tes_df <- data.frame(V1=as.data.frame(table(logistic_D3_tes_ref))[,1],V2=as.data.frame(table(logistic_D3_tes_ref))[,2],V3=rep("test",nrow(as.data.frame(table(logistic_D3_tes_ref)))))

tes_df[,1] <- as.numeric(as.character(tes_df[,1]))
tes_df <- tes_df[tes_df[,1]<201 & tes_df[,1]>-201,]   


tss_df <- rbind(data.frame(tss_df[,c(1,2)],V3="tss"),
                data.frame(tes_df[,c(1,2)],V3="tes"))

#tss_df[tss_df[,3]=="tss",2] <- tss_df[tss_df[,3]=="tss",2]/13249
#tss_df[tss_df[,3]=="tes",2] <- tss_df[tss_df[,3]=="tes",2]/20376

ggplot(tss_df,aes(V1,V2,color=V3))+ 
  
  scale_x_continuous(breaks=seq(-200,200,100))+
  
  scale_colour_manual(values=c("red", "blue"),
                      breaks=c("tes","tss"),
                      labels=c("tes","tss"))+ #改折线颜色  
  
  geom_xspline(data = tss_df,spline_shape = 1,size=1)+
  
  theme_minimal()+theme(text=element_text(size=17))+
  
  labs(x="Distance to the annotated tes",y="Number of cluster",col="cell type")+
  
  #guides(fill=guide_legend(title=NULL))+
  
  theme_set(theme_bw()) +
  
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  
  theme(panel.border = element_blank())+ #去除边框
  
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), #加上x轴
        axis.line.y=element_line(linetype=1,color="black",size=1), #加上y轴
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 10.5), #改y轴字体
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))




plot(c(1,1))
legend("topright",                                    #图例位置为右上角
       legend=c("Tokyo","NewYork","London","Berlin"),        #图例内容
       col=c("blue","red","green","orange"),                 #图例颜色
       lty=1,lwd=2)                                          #图例大小
