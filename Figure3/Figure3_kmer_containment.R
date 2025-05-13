library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

pgtmash <- read.csv("Figure3_kmer_containment.csv", header = TRUE, sep = ",")

pgtmelt1 <- melt(pgtmash, value.name = "var")

pgtmelt2 <- pgtmelt1 %>%
  separate(variable, c("hap", "type"))

pgtcast <- dcast(pgtmelt2, ... ~ type, drop = TRUE)


fig3_plot <- ggplot(data = subset(pgtcast, !(hap %in% c("hap06","hap07")))) +
  geom_point(aes(y = ident, x = shared), color = "grey30", size = 1.5, alpha = 0.6) +
  geom_point(data = subset(pgtcast, !(hap %in% c("hap06","hap07")) & ident > 99.85 & shared > 97.65), aes(y = ident, x = shared), color = "red", size = 1.5) +
  geom_vline(xintercept = 97.65, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 99.85, color = "black", linewidth = 0.5, linetype = "dashed") +
  scale_y_continuous(expand = c(0,0.07)) +
  scale_x_continuous(expand = c(0,1)) +
  facet_wrap(~ hap) +
  labs(y = expression(paste("% ", italic("k"), "-mer identity", sep = "")), x = expression(paste("% shared ", italic("k"), "-mers", sep = ""))) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "grey70"),
        panel.grid.minor = element_line(color = "grey90"),
        axis.line = element_line(color = "black"),
        strip.background = element_blank())

fig3_plot

ggsave("Figure3_kmer_containment.tiff", 
       fig3_plot, device = "tiff", dpi = 600, width = 7.5, height = 5.5, units = "in")
