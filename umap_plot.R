library(ggplot2)

if (!dir.exists("figure")) dir.create("figure")

# fig4b
data <- read.csv("data/umap.csv")
p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, )) + 
  geom_point(aes(col=var, shape = var), size=1) + 
  theme_classic()
p
ggsave(p, filename = "figure/D80_uamp.png", height = 6, width = 6.2)
ggsave(p, filename = "figure/D80_umap.pdf", height = 6, width = 6.2)

# supplementary fig7i
p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(aes(col=cell_type), size=0.1) + 
  theme_classic()
p
ggsave(p, filename = "figure/D80_cell_type.png", height = 6, width = 6.2)
ggsave(p, filename = "figure/D80_cell_type.pdf", height = 6, width = 6.2)

monocle <- read.csv("data/barcode.csv")
monocle[,4] <- monocle[,4]/max(monocle[,4])*35

monocle[monocle[,4]<5,8] <- 1
monocle[monocle[,4]>=5 & monocle[,4]<10,8] <- 2
monocle[monocle[,4]>=10 & monocle[,4]<15,8] <- 3
monocle[monocle[,4]>=15 & monocle[,4]<20,8] <- 4
monocle[monocle[,4]>=20 & monocle[,4]<25,8] <- 5
monocle[monocle[,4]>=25 & monocle[,4]<30,8] <- 6
monocle[monocle[,4]>=30 & monocle[,4]<=35,8] <- 7

monocle[,8] <- as.character(monocle[,8])

p <- ggplot(monocle, aes(x=x, y=y, colour=cell_type)) + #Pseudotime
  geom_point(size=1)+ theme_classic()
p
ggsave(p, filename = "figure/monocle_cell_type.png", height = 6, width = 6.2)
ggsave(p, filename = "figure/monocle_cell_type.pdf", height = 6, width = 6.2)

p <- ggplot(monocle, aes(x=x, y=y, colour=V8)) + #Pseudotime
  geom_point(size=1)+ theme_classic()
p
ggsave(p, filename = "figure/monocle_time_point.png", height = 6, width = 6.2)
ggsave(p, filename = "figure/monocle_time_point.pdf", height = 6, width = 6.2)
