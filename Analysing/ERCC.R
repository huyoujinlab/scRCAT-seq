options(stringsAsFactors = FALSE)
options(scipen = 100)

load("data.RData")

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






df <- merge(concentration[,c(1,3)],tss_rpm_new,by = "ERCC",all = T)   ########chahge tes or tes
df[is.na(df)] <- 0

df[,2] <- log10(df[,2]+1)
df[,3] <- log10(df[,3]+1)




cor(df[,2:3],method = "spearman")
ggplot(df, aes(x=Mix1, y=tes_rpm1))+
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
