options(stringsAsFactors = FALSE)
options(scipen = 100)

library(ggplot2)


args = commandArgs(T)

outdir <- args[5]


if (!dir.exists(outdir)) dir.create(outdir)

if (!dir.exists("outdir")) dir.create("outdir")

if (!dir.exists("outdir/novel")) dir.create("outdir/novel")


number <- c()

tss <- read.csv(args[1])
tss_intergenic <- read.csv(args[2])

number <- c(number,nrow(tss[tss$isindsc==0 & tss$model.prediction==1,])+nrow(tss_intergenic[tss_intergenic$isindsc==0 & tss_intergenic$model.prediction==1,]))

tes <- read.csv(args[3])
tes_intergenic <- read.csv(args[4])

number <- c(number,nrow(tes[tes$isinPAS==0 & tes$model.prediction==1,])+nrow(tes_intergenic[tes_intergenic$isinPAS==0 & tes_intergenic$model.prediction==1,]))





tss_intergenic <- tss_intergenic[tss_intergenic$model.prediction==1,]
tes_intergenic <- tes_intergenic[tes_intergenic$model.prediction==1,]

write.table(data.frame(tss_intergenic[,2],tss_intergenic[,7]-1,tss_intergenic[,c(7,1,1,5)]),"outdir/novel/tss_intergenic.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(data.frame(tes_intergenic[,2],tes_intergenic[,7]-1,tes_intergenic[,c(7,1,1,5)]),"outdir/novel/tes_intergenic.bed",quote = F,row.names = F,col.names = F,sep = "\t")




system("bedtools intersect -s -a outdir/novel/tss_intergenic.bed -b reference/gencode_hg38_all_gene_intergenic.bed -wa -wb > outdir/novel/tss_intergenic.bed.bed")
system("bedtools intersect -s -a outdir/novel/tes_intergenic.bed -b reference/gencode_hg38_all_gene_intergenic.bed -wa -wb > outdir/novel/tes_intergenic.bed.bed")


tss_intergenic <- read.table("outdir/novel/tss_intergenic.bed.bed")
tss_intergenic <- tss_intergenic[,c(1,2,3,10,5,6)]

tes_intergenic <- read.table("outdir/novel/tes_intergenic.bed.bed")
tes_intergenic <- tes_intergenic[,c(1,2,3,10,5,6)]


novel_gene <- merge(tss_intergenic,tes_intergenic,by = "V10")
novel_gene <- novel_gene[abs(novel_gene[,4]-novel_gene[,9])<100000,]

novel_gene_plus <- novel_gene[novel_gene[,6]=="+",]
novel_gene_minus <- novel_gene[novel_gene[,6]=="-",]


novel_gene_plus <- novel_gene_plus[novel_gene_plus[,4]<novel_gene_plus[,9],]
novel_gene_minus <- novel_gene_minus[novel_gene_minus[,4]>novel_gene_minus[,9],]

novel_gene_plus <- novel_gene_plus[,c(2,3,9,1,1,6)]
novel_gene_minus <- novel_gene_minus[,c(2,8,4,1,1,6)]

colnames(novel_gene_plus) <- colnames(novel_gene_minus)

novel_gene <- rbind(novel_gene_plus,novel_gene_minus)

write.table(novel_gene,"outdir/novel/novelgene.bed",quote = F,row.names = F,col.names = F,sep = "\t")



system("sort -k1,1 -k2,2n outdir/novel/novelgene.bed > outdir/novel/novelgene.bed.bed")
system("bedtools merge -s -i outdir/novel/novelgene.bed.bed -c 4 -o collapse > outdir/novel/novelgene.bed.bed.bed")



novel_gene <- read.table("outdir/novel/novelgene.bed.bed.bed")

for (i in 1:nrow(novel_gene)) {
  novel_gene[i,4] <- strsplit(novel_gene[i,4],split = "_")[[1]][length(strsplit(novel_gene[i,4],split = "_")[[1]])]
}



colnames(novel_gene) <- c("chr","start","end","strand")

number <- c(number,nrow(novel_gene))

novel <- data.frame(type=c("Novel TSS","Novel TES","Novel gene"),
                    number)

novel[,1] <- factor(novel[,1],levels = c("Novel TSS","Novel TES","Novel gene"))

p <- ggplot(novel,aes(x=type,y=number,fill=type)) +
  geom_bar(position=position_dodge(0.7),width=0.5,stat="identity") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  #theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(face = "plain",size = 11,angle=0,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="none") +
  labs(x=element_blank(),y="Number")+
  scale_y_continuous(expand=c(0,0))
#scale_y_continuous(breaks = seq(0,2500,500),limits = c(0,2500),expand=c(0,0),labels = c("0","500","1,000","1,500","2,000","2,500"))


ggsave(p, filename = file.path(outdir, "novel.pdf"), height = 4, width = 4)

write.table(novel_gene,paste(outdir,'/novelgene.tsv',sep=""),quote = F,row.names = F,col.names = T,sep = "\t")
