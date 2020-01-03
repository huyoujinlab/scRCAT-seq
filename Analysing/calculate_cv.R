options(stringsAsFactors = FALSE)
options(scipen = 100)
library(CAGEr)
library(genefilter)
library(statmod)
library(BuenColors)
library(DESeq2)

load("calculate_cv.RData")



counts <- data.frame(gene=OC1_1[,1],oc1_1=OC1_1[,2],oc1_2=OC1_2[,2],oc1_3=OC1_3[,2],oc10_1=OC10_1[,2],oc10_2=OC10_2[,2],oc10_3=OC10_3[,2])
counts <- counts[-c(54147:54151),]
counts <- counts[!(counts[,2]==0 & counts[,3]==0 & counts[,4]==0 & counts[,5]==0 & counts[,6]==0 & counts[,7]==0),]

counts_tss <- data.frame(gene=OC1_1_tss[,1],oc1_1=OC1_1_tss[,2],oc1_2=OC1_2_tss[,2],oc1_3=OC1_3_tss[,2],oc10_1=OC10_1_tss[,2],oc10_2=OC10_2_tss[,2],oc10_3=OC10_3_tss[,2])
counts_tss <- counts_tss[-c(54147:54151),]
counts_tss <- counts_tss[!(counts_tss[,2]==0 & counts_tss[,3]==0 & counts_tss[,4]==0 & counts_tss[,5]==0 & counts_tss[,6]==0 & counts_tss[,7]==0),]







counts <- counts[counts[,1] %in% intersect(counts[,1],counts_tss[,1]),]
counts_tss <- counts_tss[counts_tss[,1] %in% intersect(counts[,1],counts_tss[,1]),]


colnames(counts)[1] <- "gene"
rownames(counts) <- counts[,1]
counts <- counts[,-1]

colnames(counts_tss)[1] <- "gene"
rownames(counts_tss) <- counts_tss[,1]
counts_tss <- counts_tss[,-1]



sfHeLa <- estimateSizeFactorsForMatrix( counts )  #选择1 size factor

#sfHeLa <- colSums(counts)/10^6  #选择2 RPM


#sfHeLa <- data.frame(oc1_1=c(1,0),oc1_2=c(1,0),oc1_3=c(1,0),oc10_1=c(1,0),oc10_2=c(1,0),oc10_3=c(1,0)) #选择3  单纯read count
#sfHeLa <- colSums(sfHeLa)  #选择3


nCountsHeLa <- t( t(counts) / sfHeLa )


colHeLa <- "#00207040"

meansHeLa <- rowMeans( nCountsHeLa )
meansHeLa2 <- rowMeans( counts )
varsHeLa <- rowVars( nCountsHeLa )
cv2HeLa <- varsHeLa / meansHeLa^2

#cv2HeLa2 <- data.frame(V1=rownames(counts),cv2HeLa)


minMeanForFit <- unname( quantile( meansHeLa[ which( cv2HeLa > .3 ) ], .95 ) )

minMeanForFit <- rep(F,6800)


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
      log="xy", xlim = c( 1e-1, 3e2 ), ylim = c( .005, 8 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:2), c( "0.1", "1", "10", "100" ) )
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
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

sfHeLa_tss <- estimateSizeFactorsForMatrix( counts_tss )  #选择1 size factor

#sfHeLa <- colSums(counts)/10^6  #选择2 RPM


#sfHeLa <- data.frame(oc1_1=c(1,0),oc1_2=c(1,0),oc1_3=c(1,0),oc10_1=c(1,0),oc10_2=c(1,0),oc10_3=c(1,0)) #选择3  单纯read count
#sfHeLa <- colSums(sfHeLa)  #选择3


nCountsHeLa_tss <- t( t(counts_tss) / sfHeLa_tss )


colHeLa <- "#00207040"

meansHeLa_tss <- rowMeans( nCountsHeLa_tss )
varsHeLa_tss <- rowVars( nCountsHeLa_tss )
cv2HeLa_tss <- varsHeLa_tss / meansHeLa_tss^2

#cv2HeLa2 <- data.frame(V1=rownames(counts),cv2HeLa)


minMeanForFit_tss <- unname( quantile( meansHeLa_tss[ which( cv2HeLa_tss > .3 ) ], .95 ) )
minMeanForFit_tss <- rep(F,6800)


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
