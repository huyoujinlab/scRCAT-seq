options(stringsAsFactors = FALSE)
options(scipen = 100)

library(CAGEr)




load("20190316Dpool,Opool,Ppool,single-cell.RData")
load("20190310.RData")
ctss <- CTSStagCount(myCAGEsetSC5cap)

ctss_56 <- ctss[,c(-1,-2,-3)]
ctss_O <- ctss_56[,grep("O",colnames(ctss_56))]

ctss <- CTSStagCount(myCAGEsetSC5cap)
ctss_56 <- ctss[,c(-1,-2,-3)]
ctss_O <- ctss_56[,grep("O",colnames(ctss_56))]

set.seed(15)
a <- sample(1:8,size = 3)
a
ctss_O_1 <- ctss_O[,a]


a <- sample(setdiff(1:8,a),size = 3)
a
ctss_O_2 <- ctss_O[,a]


ctss_O_1 <- data.frame(ctss[,c(1,2,3)],as.numeric(rowSums(ctss_O_1)))
ctss_O_2 <- data.frame(ctss[,c(1,2,3)],as.numeric(rowSums(ctss_O_2)))

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

###20190428
ggplot(ctss_O_new[order(ctss_O_new$dens),], aes(x=V1, y=V2,col=col))+
  geom_point(size=0.9)+
  labs(x="log10(RPM+1,12 cells pooled)",y="log10(RPM+1,12 cells pooled)")+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ #去除边框
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1), 
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  scale_x_continuous(limits = c(0,4.5)) +
  scale_y_continuous(limits = c(0,4.5)) +
  scale_color_manual(values = cols)








ctss <- CTSStagCount(myCAGEsetSC3tail)
ctss_56 <- ctss[,c(-1,-2,-3)]
ctss_O <- ctss_56[,grep("O",colnames(ctss_56))]

set.seed(14)
a <- sample(1:8,size = 3)
a
ctss_O_1 <- ctss_O[,a]


a <- sample(setdiff(1:8,a),size = 3)
a
ctss_O_2 <- ctss_O[,a]


ctss_O_1 <- data.frame(ctss[,c(1,2,3)],as.numeric(rowSums(ctss_O_1)))
ctss_O_2 <- data.frame(ctss[,c(1,2,3)],as.numeric(rowSums(ctss_O_2)))

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
  labs(x="log10(RPM+1,12 cells pooled)",y="log10(RPM+1,12 cells pooled)")+
  theme_set(theme_bw()) +
  
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+ #去除边框
  theme(axis.line.x=element_line(linetype=1,color="black",size=1), 
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5), 
        axis.text.y = element_text(face = "plain",size = 10.5), 
        legend.position="none") +
  scale_x_continuous(limits = c(0,4.5)) +
  scale_y_continuous(limits = c(0,4.5)) +
  scale_color_manual(values = cols)
