options(stringsAsFactors = FALSE)
options(scipen = 100)

library(CAGEr)
library(rlist)
library(ggplot2)
library(Rmisc)
library(scales)
library(basicTrendline)
library(broom)
library(dplyr)
library(nlstools)

#############成本



############cat-seq，  1M read=1.82美元   cat-seq的转化率是0.0054，这是DRG_deep从raw read到最后bed文件read来的。5000人民币=400M
############现在LEC的数据，是达到了0.197的转化率,但猴子的数据没有不是cat-seq，意味着它基本上不用考虑5' 3'有多少tag的问题，按照以前的算，3’带tag的占10%，所以转化率是2%





catseq <- data.frame(number=c(56.0,125.2,280.1,588.2,749.5,1133.4,1377.7,1926.2,2031.2,2413.1,2829.3,3204.7,3344.3,3542.1,3859.4,4143.5,4408.2,5033.1),
                     money=c(0.01,0.02,0.04,0.08,0.10,0.16,0.20,0.30,0.32,0.40,0.50,0.60,0.64,0.70,0.80,0.90,1.00,1.28)/0.02*1.82) ##count 2

catseq <- data.frame(number=c(647.5,1225.5,1647.0,2020.5,2365.0,2650.5,2941.5,3191.5,3400.0,3577.0,3790.5,3982.0,4148.0,4374.5,4486.5,4690.5,4794.0,4973.0,5117.5,5194.5,5395.5,5488.5,5609.5,5726.0,5762.5,5912.0,6065.5,6185.5,6261.0,6352.5,6422.0,6509.5,6567.5,6644.0),
                     money=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4)/0.04384678*1.82)   #######新一批测序 0.04384678的转化率



y1 <- catseq[,1]
x1 <- catseq[,2]






###############pacbio
know_isoform <- read.table("G:/三代测序数据/单细胞三代测序/result_一致性分析/result/know_venn/trans_know.stat.xls",sep="\t",header=T)
#know_isoform <- know_isoform[-1,c(-2,-3)]

know_isoform <- know_isoform[-1,]


union(know_isoform[know_isoform[,2]==!"o",])

number <- c()

for(i in 2:9) {
  number <- c(number,length(unique(know_isoform[!know_isoform[,i]=="-",1])))
}


money <- c(0.32,0.35,0.89,0.24,1.03,3.52,0.22,1.31)*7200/8/6.87


pacbio <- data.frame(number=c(798,423,785,282,811,1200,356,958),  ######已知转录本，1ccs支持就算
                     money=c(0.32,0.35,0.89,0.24,1.03,3.52,0.22,1.31)*7200/8/6.87)

y2 <- pacbio[,1]
x2 <- pacbio[,2]



##############smart-seq2


smartseq <- data.frame(number=c(5.5,10,19.5,45,146,313,478.5,653.5,800.5,966.5,1087.5,1213),
                       money=c(0.64,1.28,2.56,5,10,20,30,40,50,60,70,80)/0.36*1.82)  #######2 count



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








我们的数据

wt是自变量，mpg是因变量




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


library(dplyr)
library(rsample)
library(broom)
library(purrr)

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














我们的数据

wt是自变量，mpg是因变量


library(ggplot2)


###去掉一个离群点
pacbio <- data.frame(wt=c(423,785,282,811,1200,356,958),  ######已知转录本，1ccs支持就算
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


library(dplyr)
library(rsample)
library(broom)
library(purrr)

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

boot_splines <- boots %>% 
  mutate(spline = map(splits, fit_spline_on_bootstrap),
         aug_train = map(spline, augment))

splines_aug <- boot_splines %>% 
  unnest(aug_train)

ggplot(splines_aug, aes(x, y)) +
  geom_point() +
  geom_line(aes(y = .fitted, group = id), alpha = 0.2)


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
  theme(panel.border = element_blank())+ #去除边框)+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ #去除边框
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), #加上x轴
        axis.line.y=element_line(linetype=1,color="black",size=1), #加上y轴
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5), #改x轴字体
        axis.text.y = element_text(face = "plain",size = 11),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="none") +  ######去除图例
  #labs(x=element_blank(),y="Number of transcriptst")+
  labs(x=element_blank(),y="Cost") +
  scale_y_continuous(breaks = seq(0,300,50),limits = c(0,300),expand=c(0,0),labels = c("0","50","100","150","200","250","300"))
#scale_y_continuous(breaks = seq(0,10000,2000),limits = c(0,11000),expand=c(0,0),labels = c("0","2,000","4,000","6,000","8,000","10,000"))+
#scale_y_continuous(breaks = seq(0,400,100),limits = c(0,400),expand=c(0,0),labels = c("0","100","200","300","400"))
