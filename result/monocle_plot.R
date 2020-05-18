library(ggplot2)
library(stringr)

cell_link_size = 0.75
x = 1
y = 2

# Fig 4b
data <- read.csv("data_df.csv")
edge <- read.csv("edge_df.csv")
data$State <- factor(data$State)
data$v7 <- factor(data$v7)
color_by <- "var"
p <- ggplot(data = data, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                      y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                      yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                      linetype = "solid", na.rm = TRUE, data = edge) + 
  geom_point(aes_string(color = color_by), 
                    size = I(1.5), na.rm = TRUE) +
  monocle:::monocle_theme_opts() + xlab(paste("Component", x)) + 
  ylab(paste("Component", y)) + 
  theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) +
  theme(legend.key = element_blank()) + 
  theme(panel.background = element_rect(fill = "white")) +
  scale_colour_manual(values=c("#31A46B", "#5F1885"))
ggsave(p, filename = "fig4b.png", height = 7, width = 6)

# Fig 4c
data <- read.csv("data_df.csv")
edge <- read.csv("edge_df.csv")
data$State <- factor(data$State)
data$v7 <- factor(data$cell_type)
color_by <- "cell_type"
p <- ggplot(data = data, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = cell_link_size, 
               linetype = "solid", na.rm = TRUE, data = edge) + 
  geom_point(aes_string(color = color_by), 
             size = I(1.5), na.rm = TRUE) +
  monocle:::monocle_theme_opts() + xlab(paste("Component", x)) + 
  ylab(paste("Component", y)) + 
  theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) +
  theme(legend.key = element_blank()) + 
  theme(panel.background = element_rect(fill = "white")) +
  scale_colour_manual(values=rev(c("#AE1515", "#B09A18", "#2919BF")))
ggsave(p, filename = "fig4c.png", height = 7, width = 6)

# supp Fig 8d
data <- read.csv("data_df.csv")
edge <- read.csv("edge_df.csv")
data$State <- factor(data$State)
data$v7 <- factor(data$cell_type)
color_by <- "Pseudotime"
p <- ggplot(data = data, aes(x = data_dim_1, y = data_dim_2)) + 
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = cell_link_size, 
               linetype = "solid", na.rm = TRUE, data = edge) + 
  geom_point(aes_string(color = color_by), 
             size = I(1.5), na.rm = TRUE) +
  monocle:::monocle_theme_opts() + xlab(paste("Component", x)) + 
  ylab(paste("Component", y)) + 
  theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) +
  theme(legend.key = element_blank()) + 
  theme(panel.background = element_rect(fill = "white")) 
ggsave(p, filename = "suppfig8d.png", height = 7, width = 6)
