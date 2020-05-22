library(ggplot2)
library(stringr)
library(dplyr)

args <- commandArgs(T)
input <- args[1]
output <- args[2]

input <- "output_threshold*"
for (dir in Sys.glob(file.path(input, "*"))){
  data <- data.frame()
  for (file in Sys.glob(file.path(dir, "model/*rf_prediction.csv"))){
    data <- rbind(data, read.csv(file))
  }
  file_names <- str_split(file, pattern="/", simplify=T) %>% 
    .[length(.)] %>% str_split("final", simplify=T) %>% 
    .[1] %>% paste0(., "prediction.tsv")
  data <- data[,c(1:2,7,5,ncol(data))]
  colnames(data) <- c("gene", "chr", "peak_position", "strand", "model_prediction")
  data <- data[data$model_prediction == 1,]
  write.table(data, file=file.path(output, file_names), quote = F, sep = "\t", row.name=F)
}

data <- data.frame()
for (file in Sys.glob(file.path(input, "*/*.csv"))){
  tmp <- read.csv(file)
  tmp_line <- tmp[1,]
  tmp_line[,"test_accuracy"] = tmp_line[, "accuracy_BM"]
  tmp_line[,"method"] = "Without ML"
  tmp <- rbind(tmp, tmp_line)
  data <- rbind(data, tmp)
}
data$method <- factor(data$method, levels = c("Without ML", "rf"))
p <- ggplot(data, aes(x=type, y = test_accuracy * 100, fill=method)) + 
  geom_text(data, mapping = aes(label=format(signif(test_accuracy * 100, 3), nsmall=1)), position=position_dodge(0.9),
            size=2, hjust=0.5, vjust=0) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") + theme_classic()
ggsave(p, filename = file.path(output, "accuracy.pdf"), height = 4, width = 5)
