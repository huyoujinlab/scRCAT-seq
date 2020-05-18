library(ggplot2)
library(stringr)
if (!dir.exists("figure")) dir.create("figure")

# figure1d
data <- read.csv("hESC_train.csv")
data$type <- factor(toupper(data$type), levels = c("TSS", "TES"))
tmp <- data[data$feature %in% c("all", "Without ML"),]
tmp$method <- factor(tmp$method, levels = c("Without ML", "knn", "lr", "rf", "svm"))
p <- ggplot(tmp, aes(x=type, y = test_accuracy * 100, fill=method)) + 
  geom_text(tmp, mapping = aes(label=format(signif(test_accuracy * 100, 4), nsmall=1)), position=position_dodge(0.9),
            size=4, hjust=0.5, vjust=0) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") +
  theme_classic() + xlab(("hESC"))
p
ggsave(p, filename = "./figure/hESC_method.png", height = 4, width = 5)
ggsave(p, filename = "./figure/hESC_method.pdf", height = 4, width = 5)




# supplement fig1b
data <- read.csv("ERCC_train.csv")
data$type <- factor(toupper(data$type), levels = c("TSS", "TES"))
tmp <- data[data$feature %in% c("all", "Without ML"),]
tmp$method <- factor(tmp$method, levels = c("Without ML", "knn", "lr", "rf", "svm"))
p <- ggplot(tmp, aes(x=type, y = test_accuracy * 100, fill=method)) + 
  geom_text(tmp, mapping = aes(label=format(signif(test_accuracy * 100, 4), nsmall=1)), position=position_dodge(0.9),
            size=4, hjust=0.5, vjust=0) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") +
  theme_classic() + xlab(("ERCC"))
p
ggsave(p, filename = "./figure/ERCC_method.png", height = 4, width = 5)
ggsave(p, filename = "./figure/ERCC_method.pdf", height = 4, width = 5)




# supplement fig1d
data <- read.csv("./otherERCC.csv")
tmp <- data[!(data$type %in% c("ERCC_tes", "ERCC_tss")),]
tmp$method <- factor(tmp$method, levels = rev(c("With ML", "Without ML")))
tmp <- tidyr::separate(tmp, col = type, into=c("cell", "tss"), sep="_")
p <- ggplot(tmp, aes(x=cell,y = test_accuracy * 100, fill = method)) + 
  geom_text(tmp, mapping = aes(label=format(signif(test_accuracy * 100, 4), nsmall=1)), position=position_dodge(0.9),
            size=4, hjust=0.5, vjust=-0.1) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") +
  theme_classic() + facet_grid(.~tss, scales = "free", space = "free") +
  theme(strip.background = element_rect(color = "white", fill = "white"), 
        panel.grid = element_blank()) + 
  scale_fill_manual(values = c("#F8766D", "#00B0F6"))
p
ggsave(p, filename = "./figure/ERCC_other_study.png", height = 5, width = 3.5)
ggsave(p, filename = "./figure/ERCC_other_study.pdf", height = 5, width = 3.5)



# supplement fig1e
tmp <- read.csv("hESC_sccat_as_model.csv")
tmp$method <- factor(tmp$method, levels = rev(c("With ML", "Without ML")))
tmp$cell <- factor(tmp$cell, levels = c("hESC", "HEK293T", "DRG", "oocyte", "scISOr-seq", "isoseq"))
tmp$tss <- factor(tmp$tss, levels = c("tss","tes"))
p <- ggplot(tmp, aes(x=cell,y = test_accuracy * 100, fill=method)) + 
  geom_text(tmp, mapping = aes(label=format(signif(test_accuracy * 100, 4), nsmall=1)), position=position_dodge(0.9),
            size=4, hjust=0.5, vjust=0) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") +
  theme_classic() + facet_grid(.~tss) +
  theme(strip.background = element_rect(color = "white", fill = "white"), 
        panel.grid = element_blank()) + 
  scale_fill_manual(values = c("#F8766D", "#00B0F6"))
p
ggsave(p, filename = "./figure/hESC_as_model_sccat.png", height = 5, width = 8.5)
ggsave(p, filename = "./figure/hESC_as_model_sccat.pdf", height = 5, width = 8.5)



# supplement fig1h
data <- read.csv("hESC_train.csv")
tmp <- data[data$method %in% c("rf", "Without ML"),]
tmp$method <- str_replace(tmp$method, "svc", "svm")
tmp$feature <- factor(tmp$feature, levels = c("Without ML", "basic", "internal", "motif", "all"))
tmp[,14] <- factor(tmp[,14],levels = c("tss","tes"))
p <- ggplot(tmp, aes(x=type, y = test_accuracy * 100, fill=feature)) + 
  geom_text(tmp, mapping = aes(label=format(signif(test_accuracy * 100, 4), nsmall=1)), position=position_dodge(0.9),
            size=4, hjust=0.5, vjust=0) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") +
  theme_classic() + xlab(("hESC"))
p
ggsave(p, filename = "./figure/hESC_feature_importance.png", height = 4, width = 5)
ggsave(p, filename = "./figure/hESC_feature_importance.pdf", height = 4, width = 5)





# supplement fig7a
tmp <- read.csv("hESC_10X_model_with_3.csv")
tmp$method <- factor(tmp$method, levels = rev(c("With ML", "Without ML")))
tmp$cell <- factor(tmp$cell, levels = c("hESC", "HEK293T", "ARPE", "organoid", "mESC"))
tmp$tss <- factor(tmp$tss, levels = c("tss","tes"))
p <- ggplot(tmp, aes(x=cell,y = test_accuracy * 100)) + 
  geom_text(tmp, mapping = aes(label=format(signif(test_accuracy * 100, 4), nsmall=1)), position=position_dodge(0.9),
            size=4, hjust=0.5, vjust=-0.1) + 
  geom_bar(stat = "identity", position = 'dodge', color = "#00B0F6", fill = "#00B0F6") + 
  ylab(label = "accuracy(%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10)) +
  theme(legend.position = "top") +
  theme_classic() + facet_grid(.~tss, scales = "free", space = "free") +
  theme(strip.background = element_rect(color = "white", fill = "white"), 
        panel.grid = element_blank())
p
ggsave(p, filename = "./figure/hESC_as_model_10x.png", height = 5, width = 3.5)
ggsave(p, filename = "./figure/hESC_as_model_10x.pdf", height = 5, width = 3.5)

