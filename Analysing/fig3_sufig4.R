### fig3a fig3b fig3f sufig6a sufig6b sufig6c sufig6d

options(stringsAsFactors = FALSE)
options(scipen = 100)


library(basicTrendline)
library(broom)
library(BuenColors)
library(CAGEr)
library(data.table)
library(DESeq2)
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

species <- "mouse"




load("fig3_sufig4.RData")




#fig3a

tss_rpm <- CTSStagCount(myCAGEsetERCC_5cap_filteryes)
tss_rpm <- tss_rpm[tss_rpm[,3]=="+",]

tss_rpm_new <- data.frame()
for(i in unique(tss_rpm[,1])) {
  tss_rpm_new <- rbind(tss_rpm_new,
                       data.frame(ERCC=i,tss_rpm1=sum(tss_rpm[tss_rpm[,1]==i,4])))
}
tss_rpm_new[,2] <- tss_rpm_new[,2]/sum(tss_rpm_new[,2])*10^6













df <- merge(concentration[,c(1,3)],tss_rpm_new,by = "ERCC",all = T)   ########change tes or tes
df[is.na(df)] <- 0

df[,2] <- log10(df[,2]+1)
df[,3] <- log10(df[,3]+1)




cor(df[,2:3],method = "spearman")
ggplot(df, aes(x=Mix1, y=tss_rpm1))+
  geom_point(size=0.9,color="blue")+
  labs(x="log10(ERCC+1)",y="Obversed level")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  geom_smooth(mapping = NULL, data = NULL, stat = "smooth",
              position = "identity", method = "lm", formula = y ~ x,
              se = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,color="red")+
  scale_x_continuous(breaks = seq(0,4,1)) +
  scale_y_continuous(limits = c(0,6),breaks = seq(0,6,1)) 



#fig3b

set.seed(15)
a <- sample(1:8,size = 3)
a
ctss_O_1 <- ctss_O_tss[,a]


a <- sample(setdiff(1:8,a),size = 3)
a
ctss_O_2 <- ctss_O_tss[,a]


ctss_O_1 <- data.frame(ctss_O_tss_first3col[,c(1,2,3)],as.numeric(rowSums(ctss_O_1)))
ctss_O_2 <- data.frame(ctss_O_tss_first3col[,c(1,2,3)],as.numeric(rowSums(ctss_O_2)))

ctss_O_1 <- ctss_O_1[!ctss_O_1[,4]==0,]
ctss_O_2 <- ctss_O_2[!ctss_O_2[,4]==0,]


write.table(ctss_O_1,"ctss_O_1.ctss",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(ctss_O_2,"ctss_O_2.ctss",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")



myCAGEsetOtss <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = c("ctss_O_1.ctss","ctss_O_2.ctss"), inputFilesType = "ctss",sampleLabels = c("O_1","O_2"))



getCTSS(myCAGEsetOtss,removeFirstG = FALSE,correctSystematicG = FALSE)


ctss_O <- CTSStagCount(myCAGEsetOtss)

print(cor(ctss_O[,c(4,5)],method = "pearson"))





ctss_O_new <- data.frame(V1=log10(ctss_O[,4]+1),V2=log10(ctss_O[,5]+1))

x <- densCols(ctss_O_new[,1],ctss_O_new[,2], colramp=colorRampPalette(c("black", "white")))
ctss_O_new$dens <- col2rgb(x)[1,] + 1L
## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)  
ctss_O_new$col <- cols[ctss_O_new$dens]
ctss_O_new <- unique(ctss_O_new)

cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(15)  


ggplot(ctss_O_new[order(ctss_O_new$dens),], aes(x=V1, y=V2,col=col))+
  geom_point(size=0.9)+
  labs(x="log10(RPM+1,3 cells pooled)",y="log10(RPM+1,3 cells pooled)")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  scale_x_continuous(limits = c(0,4.5)) +
  scale_y_continuous(limits = c(0,4.5)) +
  scale_color_manual(values = cols)




## fig3f

## fig3f 原来sufig4.RData就行了
counts <- data.frame(gene=OC1_1[,1],oc1_1=OC1_1[,2],oc1_3=OC1_3[,2],oc10_1=OC10_1[,2],oc10_3=OC10_3[,2])
counts <- counts[-c(54147:54151),]
counts <- counts[!(counts[,2]==0 & counts[,3]==0 & counts[,4]==0 & counts[,5]==0),]



# 去掉只在一个细胞表达的基因
#for (i in 1:nrow(counts)) {
#  if (!nrow(as.data.frame(table(c(counts[i,2:5]==0)))[as.data.frame(table(c(counts[i,2:5]==0)))[,1]==T,])==0) {
#    if (as.data.frame(table(c(counts[i,2:5]==0)))[as.data.frame(table(c(counts[i,2:5]==0)))[,1]==T,2]==3) {
#      counts[i,1]="."
#    }
#  }
#
#}
#counts <- counts[!counts[,1]==".",]


counts_tss <- data.frame(gene=OC1_1_tss[,1],oc1_1=OC1_1_tss[,2],oc1_2=OC1_2_tss[,2],oc1_3=OC1_3_tss[,2],oc10_1=OC10_1_tss[,2],oc10_2=OC10_2_tss[,2],oc10_3=OC10_3_tss[,2])
counts_tss <- counts_tss[-c(54147:54151),]
counts_tss <- counts_tss[!(counts_tss[,2]==0 & counts_tss[,3]==0 & counts_tss[,4]==0 & counts_tss[,5]==0 & counts_tss[,6]==0 & counts_tss[,7]==0),]

counts_tss <- data.frame(gene=OC1_1_tss[,1],oc1_1=OC1_1_tss[,2],oc1_3=OC1_3_tss[,2],oc10_1=OC10_1_tss[,2],oc10_3=OC10_3_tss[,2])
counts_tss <- counts_tss[-c(54147:54151),]
counts_tss <- counts_tss[!(counts_tss[,2]==0 & counts_tss[,3]==0 & counts_tss[,4]==0 & counts_tss[,5]==0),]







counts <- counts[counts[,1] %in% intersect(counts[,1],counts_tss[,1]),]
counts_tss <- counts_tss[counts_tss[,1] %in% intersect(counts[,1],counts_tss[,1]),]


colnames(counts)[1] <- "gene"
rownames(counts) <- counts[,1]
counts <- counts[,-1]

colnames(counts_tss)[1] <- "gene"
rownames(counts_tss) <- counts_tss[,1]
counts_tss <- counts_tss[,-1]



sfHeLa <- estimateSizeFactorsForMatrix( counts )  #size factor




nCountsHeLa <- t( t(counts) / sfHeLa )


colHeLa <- "#00207040"

meansHeLa <- rowMeans( nCountsHeLa )
meansHeLa2 <- rowMeans( counts )
varsHeLa <- rowVars( nCountsHeLa )
cv2HeLa <- varsHeLa / meansHeLa^2

#cv2HeLa2 <- data.frame(V1=rownames(counts),cv2HeLa)


minMeanForFit <- unname( quantile( meansHeLa[ which( cv2HeLa > .3 ) ], .95 ) )

minMeanForFit <- rep(F,nrow(counts))


useForFit <- meansHeLa >= minMeanForFit
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansHeLa[useForFit] ),
                   cv2HeLa[useForFit] )
fit$coefficients



xi <- mean( 1 / sfHeLa )


a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"] - xi )
c( a0, a1 )



###########################################################################################################
#
#
plot( NULL, xaxt="n", yaxt="n",                                                                           #
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),                                               #
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )          #
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",                                                   #
                       expression(10^4), expression(10^5) ) )                                             #
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )                                                #
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )                                                  #
#
#
###########################################################################################################




plot( NULL, xaxt="n", yaxt="n",
      log="x", xlim = c( 1e0, 3e3 ), ylim = c( 0, 4.2 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:3), c("0.1", "1", "10", "100", "1000" ) )
axis( 2, seq(0,4,1), c( "0", "1", "2", "3" , "4"), las=2 )


plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e3 ), ylim = c( .0005, 100 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:3), c("0.1", "1", "10", "100", "1000" ) )
axis( 2, 10^(-3:2), c("0.001","0.01", "0.1", "1", "10" ,100), las=2 )   

#abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Add the data points
#points( meansHeLa[useForFit], cv2HeLa[useForFit], pch=20, cex=.2, col="#B0E0E6" )
# Plot the fitted curve
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="blue", lwd=3 )
# Plot quantile lines around the fit
df <- ncol(counts) - 1
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df,
       col="#B0E0E6", lwd=2, lty=1 )
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df,
       col="#B0E0E6", lwd=2, lty=1 ) 




##

sfHeLa_tss <- estimateSizeFactorsForMatrix( counts_tss )  # size factor




nCountsHeLa_tss <- t( t(counts_tss) / sfHeLa_tss )


colHeLa <- "#00207040"

meansHeLa_tss <- rowMeans( nCountsHeLa_tss )
varsHeLa_tss <- rowVars( nCountsHeLa_tss )
cv2HeLa_tss <- varsHeLa_tss / meansHeLa_tss^2

#cv2HeLa2 <- data.frame(V1=rownames(counts),cv2HeLa)


minMeanForFit_tss <- unname( quantile( meansHeLa_tss[ which( cv2HeLa_tss > .3 ) ], .95 ) )
minMeanForFit_tss <- rep(F,nrow(counts))


useForFit_tss <- meansHeLa >= minMeanForFit

fit_tss <- glmgam.fit( cbind( a0_tss = 1, a1tilde = 1/meansHeLa[useForFit] ),
                       cv2HeLa_tss[useForFit] )  ##iso-seq的过滤

fit_tss$coefficients



xi_tss <- mean( 1 / sfHeLa )


a0_tss <- unname( fit_tss$coefficients["a0_tss"] )
a1_tss <- unname( fit_tss$coefficients["a1tilde"] - xi_tss )
c( a0_tss, a1_tss )


#points( meansHeLa[useForFit], cv2HeLa_tss[useForFit], pch=20, cex=.2, col="#ffb3a7" )

lines( xg, (xi_tss+a1_tss)/xg + a0_tss, col="red", lwd=3 )

df <- ncol(counts_tss) - 1
lines( xg, ( (xi_tss+a1_tss)/xg + a0_tss ) * qchisq( .975, df ) / df,
       col="#ffb3a7", lwd=2, lty=1 )

lines( xg, ( (xi_tss+a1_tss)/xg + a0_tss ) * qchisq( .025, df ) / df,
       col="#ffb3a7", lwd=2, lty=1 ) 







###sufig4a

#HEK293T_RC <- read.table("~/zjw/nc/figuregithub/HEK293T_RC.tsv",header = T)
#HEK293T_UMI <- read.table("~/zjw/nc/figuregithub/HEK293T_UMI.tsv",header = T)

HEK293T <- merge(HEK293T_UMI,HEK293T_RC,by="gene",all=F)
HEK293T[is.na(HEK293T)] <- 0


#ctss_O_new <- data.frame(V1=log10(HEK293T[,2]+1),V2=log10(HEK293T[,3]+1))

#ctss_O_new <- data.frame(V1=log10((HEK293T[,2])/sum(HEK293T[,2])*1000000+1),V2=log10(HEK293T[,3]/sum(HEK293T[,3])*1000000+1))

ctss_O_new <- data.frame(V1=log10(HEK293T[,2]+1),V2=log10(HEK293T[,3]+1))

#ctss_O_new <- data.frame(V1=HEK293T[,2],V2=HEK293T[,3])


x <- densCols(ctss_O_new[,1],ctss_O_new[,2], colramp=colorRampPalette(c("black", "white")))
ctss_O_new$dens <- col2rgb(x)[1,] + 1L
## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)  
ctss_O_new$col <- cols[ctss_O_new$dens]
ctss_O_new <- unique(ctss_O_new)

cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(63)   #length(unique(ctss_O_new[,3]))


ggplot(ctss_O_new[order(ctss_O_new$dens),], aes(x=V1, y=V2,col=col))+
  geom_point(size=0.9)+
  labs(x="log10(UMI count +1)",y="log10(Read count +1)")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  scale_x_continuous(limits = c(0,4)) +
  scale_y_continuous(limits = c(0,4)) +
  scale_color_manual(values = cols)







#sufig4b


tes_rpm <- CTSStagCount(myCAGEsetERCC_3tail_filteryes)
tes_rpm <- tes_rpm[tes_rpm[,3]=="+",]

tes_rpm_new <- data.frame()
for(i in unique(tes_rpm[,1])) {
  tes_rpm_new <- rbind(tes_rpm_new,
                       data.frame(ERCC=i,tes_rpm1=sum(tes_rpm[tes_rpm[,1]==i,4])))
}
tes_rpm_new[,2] <- tes_rpm_new[,2]/sum(tes_rpm_new[,2])*10^6





df <- merge(concentration[,c(1,3)],tes_rpm_new,by = "ERCC",all = T)   ########tes
df[is.na(df)] <- 0

df[,2] <- log10(df[,2]+1)
df[,3] <- log10(df[,3]+1)




cor(df[,2:3],method = "spearman")
ggplot(df, aes(x=Mix1, y=tes_rpm1))+
  geom_point(size=0.9,color="green")+
  labs(x="log10(ERCC+1)",y="Obversed level")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  geom_smooth(mapping = NULL, data = NULL, stat = "smooth",
              position = "identity", method = "lm", formula = y ~ x,
              se = F, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,color="red")+
  scale_x_continuous(breaks = seq(0,4,1)) +
  scale_y_continuous(limits = c(0,6),breaks = seq(0,6,1)) 
















#sufig4c

set.seed(14)
a <- sample(1:8,size = 3)
a
ctss_O_1 <- ctss_O_tes[,a]


a <- sample(setdiff(1:8,a),size = 3)
a
ctss_O_2 <- ctss_O_tes[,a]


ctss_O_1 <- data.frame(ctss_O_tes_first3col[,c(1,2,3)],as.numeric(rowSums(ctss_O_1)))
ctss_O_2 <- data.frame(ctss_O_tes_first3col[,c(1,2,3)],as.numeric(rowSums(ctss_O_2)))

ctss_O_1 <- ctss_O_1[!ctss_O_1[,4]==0,]
ctss_O_2 <- ctss_O_2[!ctss_O_2[,4]==0,]


write.table(ctss_O_1,"ctss_O_1.ctss",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
write.table(ctss_O_2,"ctss_O_2.ctss",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")



myCAGEsetOtss <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",inputFiles = c("ctss_O_1.ctss","ctss_O_2.ctss"), inputFilesType = "ctss",sampleLabels = c("O_1","O_2"))



getCTSS(myCAGEsetOtss,removeFirstG = FALSE,correctSystematicG = FALSE)


ctss_O <- CTSStagCount(myCAGEsetOtss)

cor(ctss_O[,c(4,5)],method = "pearson")





ctss_O_new <- data.frame(V1=log10(ctss_O[,4]+1),V2=log10(ctss_O[,5]+1))

x <- densCols(ctss_O_new[,1],ctss_O_new[,2], colramp=colorRampPalette(c("black", "white")))
ctss_O_new$dens <- col2rgb(x)[1,] + 1L
## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)  
ctss_O_new$col <- cols[ctss_O_new$dens]
ctss_O_new <- unique(ctss_O_new)

cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(17)  


ggplot(ctss_O_new[order(ctss_O_new$dens),], aes(x=V1, y=V2,col=col))+
  geom_point(size=0.9)+
  labs(x="log10(RPM+1,3 cells pooled)",y="log10(RPM+1,3 cells pooled)")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  scale_x_continuous(limits = c(0,4.5)) +
  scale_y_continuous(limits = c(0,4.5)) +
  scale_color_manual(values = cols)



### sufig4 d

#rm(list = ls())

#load("sufig4.RData")



counts <- data.frame(gene=DRG_1[,1],drg_1=DRG_1[,2],drg_2=DRG_2[,2])
counts <- counts[-c(54147:54151),]
counts <- counts[!(counts[,2]==0 & counts[,3]==0),]




for(i in 1:nrow(counts)) {
  if(counts[i,1] %in% ENSM2symbol[,1]) counts[i,1] <- ENSM2symbol[ENSM2symbol[,1]==counts[i,1],2]
}

colnames(counts)[1] <- "gene"




smart <- data.frame(gene=rownames(smart),read_count=smart[,8]) #if drg



smart <- smart[1:54028,]
smart <- smart[!smart[,2]==0,]
smart <- smart[!smart[,1] %in% c("__no_feature","__ambiguous"),]












counts <- counts[!counts[,2]==0,]




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












#D_dominant_tss_in_gene_1 <- cbind(read.csv("C:/Users/zhong/Desktop/201907211815/20190711DRGtss_peak_new.csv"),
#                                  read.csv("C:/Users/zhong/Desktop/201907211815/majority_vote_tss_test_DRG.csv")[,c(12,13)])
#D_dominant_tss_in_gene_1 <- D_dominant_tss_in_gene_1[D_dominant_tss_in_gene_1[,17]==1,]
#D_dominant_tss_in_gene_1 <- data.frame(V1=D_dominant_tss_in_gene_1[,2],V2=D_dominant_tss_in_gene_1[,7]-1,V3=D_dominant_tss_in_gene_1[,7],
#                                       V10=D_dominant_tss_in_gene_1[,1],V5=D_dominant_tss_in_gene_1[,6],V6=D_dominant_tss_in_gene_1[,5])




#D_dominant_tes_in_gene_1 <- cbind(read.csv("C:/Users/zhong/Desktop/201907211815/20190711DRGtes_peak_new.csv"),
#                                  read.csv("C:/Users/zhong/Desktop/201907211815/majority_vote_tes_test_DRG.csv")[,c(12,13)])
#D_dominant_tes_in_gene_1 <- D_dominant_tes_in_gene_1[D_dominant_tes_in_gene_1[,17]==1,]
#D_dominant_tes_in_gene_1 <- data.frame(V1=D_dominant_tes_in_gene_1[,2],V2=D_dominant_tes_in_gene_1[,7]-1,V3=D_dominant_tes_in_gene_1[,7],
#                                       V10=D_dominant_tes_in_gene_1[,1],V5=D_dominant_tes_in_gene_1[,6],V6=D_dominant_tes_in_gene_1[,5])











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

