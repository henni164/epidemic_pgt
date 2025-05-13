options(scipen = 999)
library(ggplot2)
library(reshape2)
library(egg)


chr_dat <- read.delim("ETH2013_ITA2018_chrlengths.txt", sep = "\t", header = TRUE)
colnames(chr_dat) <- c("Haplotype", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")
chr_dat$Haplotype <- as.factor(chr_dat$Haplotype)
chr_dat$Haplotype <- factor(chr_dat$Haplotype, levels = c("hap01","hap02","hap03","hap04","hap05","hap06","hap07"))

chr_melted <- melt(chr_dat, value.name = "Length")
colnames(chr_melted)[2] <- "Chr"
chr_melted$Length <- chr_melted$Length/1000000

chr_lengths_plot <- ggplot(chr_melted) +
  geom_boxplot(aes(x = Chr, y = Length), color = "grey40", size = 0.3) +
  geom_jitter(aes(x = Chr, y = Length, fill = Haplotype), color = "black", shape = 21, size = 2, alpha = 0.9, width = 0.1, height = 0) +
  scale_fill_manual(values = c("#BDD7EE","#F5B184","#6fae47","#F08984","#AA9DCC","#EFB0D0","#D0E1A2"), breaks = c("hap01","hap02","hap03","hap04","hap05","hap06","hap07")) + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  labs(x = "Chromosome", y = "Length (Mbp)") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey90"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(color = "black"),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.key = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size =8),
        legend.justification = "top",
        legend.position = "right",
        legend.box.margin = margin(0,0,-1,0),
        legend.margin = margin(1,1.2,0,1.2),
        legend.box.spacing = unit(-1.5,"lines"))

chr_lengths_plot

ggsave("Pgt_chrlengths.tiff", 
       chr_lengths_plot, device = "tiff", height = 5, width = 5, units = "in", dpi = 600)
