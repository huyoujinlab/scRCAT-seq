# fig1 sufig2 
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



setwd("~/zjw/scCAT-seq-master/Analysing/")
load("fig1_sufig2.RData")

### fig1 e


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




### fig1 f, sufig2 a

tss_rpm <- CTSStagCount(myCAGEsetERCC_5cap_filteryes)
tss_rpm <- tss_rpm[tss_rpm[,3]=="+",]

tss_rpm_new <- data.frame()
for(i in unique(tss_rpm[,1])) {
  tss_rpm_new <- rbind(tss_rpm_new,
                       data.frame(ERCC=i,tss_rpm1=sum(tss_rpm[tss_rpm[,1]==i,4])))
}
tss_rpm_new[,2] <- tss_rpm_new[,2]/sum(tss_rpm_new[,2])*10^6







tes_rpm <- CTSStagCount(myCAGEsetERCC_3tail_filteryes)
tes_rpm <- tes_rpm[tes_rpm[,3]=="+",]

tes_rpm_new <- data.frame()
for(i in unique(tes_rpm[,1])) {
  tes_rpm_new <- rbind(tes_rpm_new,
                       data.frame(ERCC=i,tes_rpm1=sum(tes_rpm[tes_rpm[,1]==i,4])))
}
tes_rpm_new[,2] <- tes_rpm_new[,2]/sum(tes_rpm_new[,2])*10^6






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



### fig1 g, sufig2 b

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



### fig1 h


catseq <- data.frame(number=c(647.5,1225.5,1647.0,2020.5,2365.0,2650.5,2941.5,3191.5,3400.0,3577.0,3790.5,3982.0,4148.0,4374.5,4486.5,4690.5,4794.0,4973.0,5117.5,5194.5,5395.5,5488.5,5609.5,5726.0,5762.5,5912.0,6065.5,6185.5,6261.0,6352.5,6422.0,6509.5,6567.5,6644.0),
                     money=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4)/0.04384678*1.82)   #######新一批测序 0.04384678的转化率



y1 <- catseq[,1]
x1 <- catseq[,2]






###############pacbio
#know_isoform <- read.table("G:/三代测序数据/单细胞三代测序/result_一致性分析/result/know_venn/trans_know.stat.xls",sep="\t",header=T)


#know_isoform <- know_isoform[-1,]


#union(know_isoform[know_isoform[,2]==!"o",])

#number <- c()

#for(i in 2:9) {
#  number <- c(number,length(unique(know_isoform[!know_isoform[,i]=="-",1])))
#}


money <- c(0.32,0.35,0.89,0.24,1.03,3.52,0.22,1.31)*7200/8/6.87


pacbio <- data.frame(number=c(798,423,785,282,811,1200,356,958), 
                     money=c(0.32,0.35,0.89,0.24,1.03,3.52,0.22,1.31)*7200/8/6.87)

y2 <- pacbio[,1]
x2 <- pacbio[,2]



##############smart-seq2


smartseq <- data.frame(number=c(5.5,10,19.5,45,146,313,478.5,653.5,800.5,966.5,1087.5,1213),
                       money=c(0.64,1.28,2.56,5,10,20,30,40,50,60,70,80)/0.36*1.82)  


y3 <- smartseq[,1]
x3 <- smartseq[,2]

#### http://blog.sciencenet.cn/blog-651374-1126673.html
x <- c(1, 3, 6,  9,  13,   17)
y <- c(5, 8, 11, 13, 13.2, 13.5)
trendline(x, y, model="power2P", eDigit = 3, eSize = 1.4, text.col = "blue")

####x1x2x3 y1y2y3
trendline(x1, y1, model="power2P", ePos.x = "topleft", summary=TRUE, 
          eDigit=5,show.Rpvalue = FALSE,xlab = "cost",ylab = "number",xlim=c(0,500),ylim=c(0,7000))

par(new=TRUE)


trendline(x2, y2, model="power2P", ePos.x = "topleft", summary=TRUE, 
          eDigit=5,show.Rpvalue = FALSE,xlab = "cost",ylab = "number",xlim=c(0,500),ylim=c(0,7000))

par(new=TRUE)

trendline(x3, y3, model="power2P", ePos.x = "topleft", summary=TRUE, 
          eDigit=5,show.Rpvalue = FALSE,xlab = "cost",ylab = "number",xlim=c(0,500),ylim=c(0,7000))













catseq <- data.frame(wt=c(647.5,1225.5,1647.0,2020.5,2365.0,2650.5,2941.5,3191.5,3400.0,3577.0,3790.5,3982.0,4148.0,4374.5,4486.5,4690.5,4794.0,4973.0,5117.5,5194.5,5395.5,5488.5,5609.5,5726.0,5762.5,5912.0,6065.5,6185.5,6261.0,6352.5,6422.0,6509.5,6567.5,6644.0),
                     mpg=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4)/0.04384678*1.82)   #######新一批测序 0.04384678的转化率




y <- catseq[,1]
x <- catseq[,2]

mtcars <- catseq

ggplot(mtcars, aes( wt,mpg)) + 
  geom_point()


nlsfit <- nls(mpg ~ a*wt^b, mtcars, start = list(a = 0.000002652, b = 1.5195279))
summary(nlsfit)


ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  geom_line(aes(y = predict(nlsfit)))




set.seed(15)






boots <- bootstraps(mtcars, times = 200)
boots



fit_nls_on_bootstrap <- function(split) {
  nls(mpg ~ a*wt^b, analysis(split), start = list(a = 0.000002652, b = 1.5195279))
}

boot_models <- boots %>% 
  mutate(model = map(splits, fit_nls_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- boot_models %>% 
  unnest(coef_info)



alpha <- .05
boot_coefs %>% 
  group_by(term) %>%
  summarize(low = quantile(estimate, alpha / 2),
            high = quantile(estimate, 1 - alpha / 2))




ggplot(boot_coefs, aes(estimate)) + 
  geom_histogram(binwidth = 2) + 
  facet_wrap(~ term, scales = "free")




boot_aug <- boot_models %>% 
  mutate(augmented = map(model, augment)) %>% 
  unnest(augmented)



ggplot(boot_aug, aes(wt, mpg)) +
  geom_point() +
  geom_line(aes(y = .fitted, group = id), alpha=.2)





fit_spline_on_bootstrap <- function(split) {
  data <- analysis(split)
  smooth.spline(data$wt, data$mpg, df = 4)
}

boot_splines <- boots %>% 
  mutate(spline = map(splits, fit_spline_on_bootstrap),
         aug_train = map(spline, augment))

splines_aug <- boot_splines %>% 
  unnest(aug_train)

ggplot(splines_aug, aes(x, y)) +
  geom_point() +
  geom_line(aes(y = .fitted, group = id), alpha = 0.2)


y=ax^b
df <- as.data.frame(boot_coefs)
df <- data.frame(a=df[df[,4]=="a",5],
                 b=df[df[,4]=="b",5])
df[,3] <- df[,1]*(1000^df[,2])
df


quantile(df[,3],probs = seq(0,1,0.025))



















###remove 1 outlier
pacbio <- data.frame(wt=c(423,785,282,811,1200,356,958),  
                     mpg=c(0.35,0.89,0.24,1.03,3.52,0.22,1.31)*7200/8/6.87)

y2 <- pacbio[,1]
x2 <- pacbio[,2]

mtcars <- pacbio

ggplot(mtcars, aes( wt,mpg)) + 
  geom_point()


nlsfit <- nls(mpg ~ a*wt^2.69295, mtcars, start = list(a = 1))
summary(nlsfit)


ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  geom_line(aes(y = predict(nlsfit)))




set.seed(15)






boots <- bootstraps(mtcars, times = 200)
boots



fit_nls_on_bootstrap <- function(split) {
  nls(mpg ~ a*wt^2.69295, analysis(split), start = list(a = 0.0002652))
}

boot_models <- boots %>% 
  mutate(model = map(splits, fit_nls_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- boot_models %>% 
  unnest(coef_info)



alpha <- .05
boot_coefs %>% 
  group_by(term) %>%
  summarize(low = quantile(estimate, alpha / 2),
            high = quantile(estimate, 1 - alpha / 2))




ggplot(boot_coefs, aes(estimate)) + 
  geom_histogram(binwidth = 2) + 
  facet_wrap(~ term, scales = "free")




boot_aug <- boot_models %>% 
  mutate(augmented = map(model, augment)) %>% 
  unnest(augmented)



ggplot(boot_aug, aes(wt, mpg)) +
  geom_point() +
  geom_line(aes(y = .fitted, group = id), alpha=.2)





fit_spline_on_bootstrap <- function(split) {
  data <- analysis(split)
  smooth.spline(data$wt, data$mpg, df = 4)
}

#boot_splines <- boots %>% 
#  mutate(spline = map(splits, fit_spline_on_bootstrap),
#         aug_train = map(spline, augment))

#splines_aug <- boot_splines %>% 
#  unnest(aug_train)

#ggplot(splines_aug, aes(x, y)) +
#  geom_point() +
#  geom_line(aes(y = .fitted, group = id), alpha = 0.2)


y=ax^2
df <- as.data.frame(boot_coefs)
df <- data.frame(a=df[df[,4]=="a",5])
df[,2] <- df[,1]*(1000^2.69295)
df


quantile(df[,2],probs = seq(0,1,0.025))







(2.913441+3.630661)/2
(199.5905+278.9891)/2

(2.913441-3.630661)/2
(199.5905-278.9891)/2


a <- data.frame(V1=c(3.272051,239.2898),V2=c("scCAT-seq","ISO-seq"),se=c(0.35861,39.6993))   
a[,2] <- factor(a[,2],levels = c("scCAT-seq","ISO-seq"),ordered = T)
a <- a[order(a[,2]),]
ggplot(a,aes(x=V2,y=V1,fill=V2)) +
  geom_bar(position=position_dodge(0.7),width=0.5,stat="identity") +
  geom_errorbar(aes(ymin=V1-se, ymax=V1+se), width=0.1,size=0.8) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ 
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="none") +  
  #labs(x=element_blank(),y="Number of transcriptst")+
  labs(x=element_blank(),y="Cost") +
  scale_y_continuous(breaks = seq(0,300,50),limits = c(0,300),expand=c(0,0),labels = c("0","50","100","150","200","250","300"))
#scale_y_continuous(breaks = seq(0,10000,2000),limits = c(0,11000),expand=c(0,0),labels = c("0","2,000","4,000","6,000","8,000","10,000"))+
#scale_y_continuous(breaks = seq(0,400,100),limits = c(0,400),expand=c(0,0),labels = c("0","100","200","300","400"))






