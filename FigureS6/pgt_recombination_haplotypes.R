options(scipen = 999)
library(ggplot2)
library(dplyr)
library(reshape2)
library(purrr)
library(ggchicklet)


chr_labels <- c(
  chr1 = "1",
  chr2 = "2",
  chr3 = "3",
  chr4 = "4",
  chr5 = "5",
  chr6 = "6",
  chr7 = "7",
  chr8 = "8",
  chr9 = "9",
  chr10 = "10",
  chr11 = "11",
  chr12 = "12",
  chr13 = "13",
  chr14 = "14",
  chr15 = "15",
  chr16 = "16",
  chr17 = "17",
  chr18 = "18"
)

haps_chrs <- read.delim("pangenome_chromosome_lengths.txt", sep = "\t", header = TRUE)
colnames(haps_chrs) <- c("HAP","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18")
chrlen_melt <-melt(haps_chrs, id.vars = "HAP", variable.name = "CHROM", value.name = "len")
chrlen_melt$CHROM <- factor(chrlen_melt$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

merged_intervals <- read.delim("all_intervals.txt", sep = "\t", header = FALSE)

colnames(merged_intervals) <- c("CHROM","start","end","HAP","dens","REFHAP")
merged_intervals <- merged_intervals[,c(4,1,2,3,6)]

final_joined <- full_join(merged_intervals, chrlen_melt, by = c("HAP","CHROM"))
final_joined$HAP <- factor(final_joined$HAP)
keepme <- final_joined

keepme$HAP <- factor(keepme$HAP, levels = c("hap01","hap02","hap03","hap04","hap05","hap06","hap07"))
keepme$CHROM <- factor(keepme$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
chrlen_melt$HAP <- factor(chrlen_melt$HAP, levels = c("hap01","hap02","hap03","hap04","hap05","hap06","hap07"))
keepme$REFHAP <- factor(keepme$REFHAP)

plot2 <- ggplot(data = keepme, aes(y = HAP)) +
  geom_point(aes(x = len+15000), y = -0.5, color = "white") +
  geom_point(aes(x = len+15000), y = 1.5, color = "white") +
  geom_point(x = -15000, y = -0.5, color = "white") +
  geom_point(x = -15000, y = 1.5, color = "white") +
  ggchicklet:::geom_rrect(data = chrlen_melt, aes(xmax = len, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.1, 'npc')) +
  geom_segment(data = na.omit(keepme), aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 6, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= chrlen_melt, aes(xmax = len+15000), xmin = -14999, ymin = -0.01, ymax = 1.01, r = unit(0.1, 'npc'), color = "black", fill = NA, size = 0.3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "fixed", space = "free_x", switch = "y", margins = FALSE) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap01","hap02","hap03","hap04","hap05","hap06","hap07"),
                     values = c("#BDD7EE","#F5B184","#6FAE47","#F08984","#AA9DCC","#EFB0D0","#D0E1A2"), name = "Color key", drop = FALSE) +
  scale_fill_manual(breaks = c("hap01","hap02","hap03","hap04","hap05","hap06","hap07"),
                    values = c("#BDD7EE","#F5B184","#6FAE47","#F08984","#AA9DCC","#EFB0D0","#D0E1A2"), name = "Color key", drop = TRUE, guide = "none") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
        strip.text.x.top = element_text(size = 10, vjust = 1, hjust = 0.5, margin = margin(t = 1, b = 0, l = 0, r = 0)),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.1), units = "in"),
        plot.margin = unit(c(0,0,1.1,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = c(0.15,-0.04),
        legend.background = element_blank(),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1, fill = NA, linetype = 0)))

ggsave("FigureS6_recombination.tiff", plot2, device = "tiff", width = 7.5, height = 4, units = "in", dpi = 600)


