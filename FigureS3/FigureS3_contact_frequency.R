library(ggplot2)

pgt_contacts <- read.csv("ETH_ITA_chr_contacts.csv", header = TRUE)

pgt_contacts$chr <- factor(pgt_contacts$chr, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"),
                           labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))

pgt_contacts$chrhap <- factor(pgt_contacts$chrhap, levels = c("hap04", "hap05", "hap06", "hap07"))
pgt_contacts$binhap <- factor(pgt_contacts$binhap, levels = c("hap04", "hap05", "hap06", "hap07"))
pgt_contacts$type <- factor(pgt_contacts$type, levels = c("trans", "cistrans"), labels = c("trans","cis + trans"))
pgt_contacts$isolate <- factor(pgt_contacts$isolate, levels = c("ETH2013", "ITA2018"))

contacts_plot <- ggplot(pgt_contacts) +
  geom_col(aes(x = chr, y = perc, group = binhap, fill = binhap), position = "stack", color = "black", linewidth = 0.2) + 
  scale_fill_manual(values = c("#F08984","#AA9DCC","#EFB0D0","#D0E1A2"), labels = c("hap04","hap05","hap06","hap07"), name = NULL) + 
  scale_y_continuous(expand = c(0,0.5)) +
  facet_grid(isolate + chrhap ~ type, scales = "free_x", drop = TRUE) + 
  labs(x = "Chromosome", y = "Hi-C contacts (%)") +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        axis.line = element_line(color = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.1, "in"),
        strip.text.y = element_text(size = 8, color = "black"),
        strip.text.x = element_text(size = 8, color = "black", face = "italic"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        panel.spacing.y = unit(1, "lines"))

ggsave("FigureS3.tiff",
       contacts_plot, device = "tiff", dpi = 600, width = 6.5, height = 5.5, units = "in")
