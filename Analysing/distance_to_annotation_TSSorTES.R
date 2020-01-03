options(stringsAsFactors = FALSE)
options(scipen = 100)


library(CAGEr)
library(splines)
library(data.table)
library(ggplot2)
library(ggalt)
library(ggsignif)
library(gridExtra)

load("distance_to_annotation_TSSorTES.RData")


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
                      labels=c("tes","tss"))+ 
  
  geom_xspline(data = tss_df,spline_shape = 1,size=1)+
  
  theme_minimal()+theme(text=element_text(size=17))+
  
  labs(x="Distance to the annotated tes",y="Number of cluster",col="cell type")+
  
  #guides(fill=guide_legend(title=NULL))+
  
  theme_set(theme_bw()) +
  
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  
  theme(panel.border = element_blank())+
  
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))

