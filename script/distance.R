options(stringsAsFactors = FALSE)
options(scipen = 100)

library(basicTrendline)
library(broom)
library(BuenColors)
library(data.table)
library(dplyr)
library(genefilter)
library(ggalt)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gplots)
library(gridExtra)
library(heatmap.plus)
library(nlstools)
library(pheatmap)
library(purrr)
library(RColorBrewer)
library(reshape)
library(rlist)
library(Rmisc)
library(rsample)
library(rstatix)
library(scales)
library(splines)
library(statmod)
library(stringr)

args = commandArgs(T)

tss_refer <- read.table("reference/gencode_hg38_all_gene_all_transcript_tss.bed",sep="\t")
tes_refer <- read.table("reference/gencode_hg38_all_gene_all_transcript_tes.bed",sep="\t")


hESC_tss <- read.csv(args[1])
hESC_tss <- hESC_tss[hESC_tss$model.prediction==1,]
hESC_tss <- data.frame(hESC_tss[,2],hESC_tss[,7]-1,hESC_tss[,c(7,1,6,5)])

hESC_tes <- read.csv(args[2])
hESC_tes <- hESC_tes[hESC_tes$model.prediction==1,]
hESC_tes <- data.frame(hESC_tes[,2],hESC_tes[,7]-1,hESC_tes[,c(7,1,6,5)])
outdir <- args[3]
if (!dir.exists(outdir)) dir.create(outdir)

hESC_tss_ref <- c()
for(i in 1:nrow(hESC_tss)) {
  b <- tss_refer[tss_refer[,1]==hESC_tss[i,1] & tss_refer[,6]==hESC_tss[i,6],]
  if(hESC_tss[i,6]=="+") {
    hESC_tss_ref <- c(hESC_tss_ref,(hESC_tss[i,3]-b[,3])[order(abs(hESC_tss[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    hESC_tss_ref <- c(hESC_tss_ref,(b[,3]-hESC_tss[i,3])[order(abs(b[,3]-hESC_tss[i,3]),decreasing = FALSE)][1])
  }
}


tss_df <- data.frame(V1=as.data.frame(table(hESC_tss_ref))[,1],V2=as.data.frame(table(hESC_tss_ref))[,2],V3=rep("test",nrow(as.data.frame(table(hESC_tss_ref)))))

tss_df[,1] <- as.numeric(as.character(tss_df[,1]))
tss_df <- tss_df[tss_df[,1]<201 & tss_df[,1]>-201,] 



hESC_tes_ref <- c()
for(i in 1:nrow(hESC_tes)) {
  b <- tes_refer[tes_refer[,1]==hESC_tes[i,1] & tes_refer[,6]==hESC_tes[i,6],]
  if(hESC_tes[i,6]=="+") {
    hESC_tes_ref <- c(hESC_tes_ref,(hESC_tes[i,3]-b[,3])[order(abs(hESC_tes[i,3]-b[,3]),decreasing = FALSE)][1])
  }else {
    hESC_tes_ref <- c(hESC_tes_ref,(b[,3]-hESC_tes[i,3])[order(abs(b[,3]-hESC_tes[i,3]),decreasing = FALSE)][1])
  }
}

tes_df <- data.frame(V1=as.data.frame(table(hESC_tes_ref))[,1],V2=as.data.frame(table(hESC_tes_ref))[,2],V3=rep("test",nrow(as.data.frame(table(hESC_tes_ref)))))

tes_df[,1] <- as.numeric(as.character(tes_df[,1]))
tes_df <- tes_df[tes_df[,1]<201 & tes_df[,1]>-201,]   


tss_df <- rbind(data.frame(tss_df[,c(1,2)],V3="tss"),
                data.frame(tes_df[,c(1,2)],V3="tes"))

#tss_df[tss_df[,3]=="tss",2] <- tss_df[tss_df[,3]=="tss",2]/13249
#tss_df[tss_df[,3]=="tes",2] <- tss_df[tss_df[,3]=="tes",2]/20376

p <- ggplot(tss_df,aes(V1,V2,color=V3))+ 
  
  scale_x_continuous(breaks=seq(-200,200,100))+
  
  scale_colour_manual(values=c("red", "blue"),
                      breaks=c("tes","tss"),
                      labels=c("tes","tss"))+ 
  
  geom_xspline(data = tss_df,spline_shape = 1,size=1)+
  
  theme_minimal()+theme(text=element_text(size=17))+
  
  labs(x="Distance to the annotated tss/tes",y="Number of cluster",col="cell type")+
  
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

ggsave(p, filename = file.path(outdir, "distance.pdf"), height = 4, width = 4)

