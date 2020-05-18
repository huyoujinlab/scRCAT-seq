library(ggplot2)

if (!dir.exists("figure")) dir.create("figure")

# fig4b
data <- read.csv("umap.csv")
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

