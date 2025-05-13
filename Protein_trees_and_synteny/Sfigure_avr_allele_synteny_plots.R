library(gggenomes)
library(svglite)
library(tibble)
library(dplyr)

palette <- c("#F8766D","#d3ff00","#7CAE00","#31e08f","#00BFC4","#0031ff","#C77CFF","#ff00b1","#00ffce","#d5d5d5","#000000")

## AvrSr13

seqs <- read_seqs("AvrSr13_intervals.fasta")
links <- read_paf("AvrSr13_intervals.paf")
gff <- read_gff3("AvrSr13.gff")
interval_coords <- read.table("AvrSr13_intervals_for_gggenomes.txt", header = FALSE, sep="\t")
colnames(interval_coords) <- c("seq_id","fullstart","fullend")


info <- left_join(gff, interval_coords, by = "seq_id")
info$newstart <- info$start - info$fullstart
info$newend <- info$end - info$fullstart

gff$start <- info$newstart
gff$end <- info$newend

bin_id <- c("hap01","hap03","hap04","hap05","hap06","hap07","hap02")
seq_id <- c("chr1_hap01","chr1_hap03","chr1_hap04","chr1_hap05","chr1_hap06","chr1_hap07","chr1_hap02")
length <- c(seqs$length[1],seqs$length[3], seqs$length[4], seqs$length[5], seqs$length[6], seqs$length[7], seqs$length[2])
bins <- data.frame(bin_id, seq_id, length)

links_1000bp <- filter(links, map_length >= 1000)

AvrSr13_plot <- gggenomes(
  genes = gff,
  seqs = bins,
  links = links_1000bp,
  .id = "file_id",
  spacing = 0.05,
  wrap = NULL,
  adjacent_only = TRUE, #will also plot links b/w non-adjacent seqs
  infer_bin_id = seq_id,
  infer_start = min(start, end),
  infer_end = max(start, end),
  infer_length = max(start, end),
  theme = c("clean", NULL),
  .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_link(color = "#d6f0eb", fill = "#d6f0eb") +
  geom_gene(aes(fill=name), color = "black", size = 3, show.legend = FALSE) + geom_gene_label(aes(label = name), size = 3, angle = 0, hjust = -0.05, vjust = -0.05, fontface = "italic") + 
  geom_bin_label(size = 3, hjust = 0.4) + scale_fill_manual(values = palette[1:nrow(gff)])

#AvrSr13_plot

ggsave(file = "AvrSr13_synteny_plot.tiff",
       AvrSr13_plot, device = "tiff", dpi = 600, height = 10, width = 18, units = "cm")



## AvrSr22

seqs <- read_seqs("AvrSr22_intervals.fasta")
links <- read_paf("AvrSr22_intervals.paf")
gff <- read_gff3("AvrSr22.gff")
interval_coords <- read.table("AvrSr22_intervals_for_gggenomes.txt", header = FALSE, sep="\t")
colnames(interval_coords) <- c("seq_id","fullstart","fullend")


info <- left_join(gff, interval_coords, by = "seq_id")
info$newstart <- info$start - info$fullstart
info$newend <- info$end - info$fullstart

gff$start <- info$newstart
gff$end <- info$newend

bin_id <- c("hap01","hap02","hap03","hap04","hap05","hap06","hap07")
seq_id <- c("chr16_hap01","chr16_hap02","chr16_hap03","chr16_hap04","chr16_hap05","chr16_hap06","chr16_hap07")
length <- c(seqs$length[1],seqs$length[2], seqs$length[3], seqs$length[4], seqs$length[5], seqs$length[6], seqs$length[7])
bins <- data.frame(bin_id, seq_id, length)

links_1000bp <- filter(links, map_length >= 1000)

AvrSr22_plot <- gggenomes(
  genes = gff,
  seqs = bins,
  links = links_1000bp,
  .id = "file_id",
  spacing = 0.05,
  wrap = NULL,
  adjacent_only = TRUE, #will also plot links b/w non-adjacent seqs
  infer_bin_id = seq_id,
  infer_start = min(start, end),
  infer_end = max(start, end),
  infer_length = max(start, end),
  theme = c("clean", NULL),
  .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_link(color = "#d6f0eb", fill = "#d6f0eb") +
  geom_gene(aes(fill=name), color = "black", size = 1.5, show.legend = FALSE) + geom_gene_label(aes(label = name), size = 3, angle = 0, hjust = -0.05, vjust = -0.05, fontface = "italic") + 
  geom_bin_label(size = 3, hjust = 0.4) + scale_fill_manual(values = palette[1:nrow(gff)])

#AvrSr22_plot

ggsave(file = "AvrSr22_synteny_plot.tiff",
       AvrSr22_plot, device = "tiff", dpi = 600, height = 10, width = 18, units = "cm")



## AvrSr27

seqs <- read_seqs("AvrSr27_intervals.fasta")
links <- read_paf("AvrSr27_intervals.paf")
gff <- read_gff3("AvrSr27.gff")
interval_coords <- read.table("AvrSr27_intervals_for_gggenomes.txt", header = FALSE, sep="\t")
colnames(interval_coords) <- c("seq_id","fullstart","fullend")


info <- left_join(gff, interval_coords, by = "seq_id")
info$newstart <- info$start - info$fullstart
info$newend <- info$end - info$fullstart

gff$start <- info$newstart
gff$end <- info$newend

bin_id <- c("hap01","hap05","hap02","hap03","hap04","hap06","hap07")
seq_id <- c("chr2_hap01","chr2_hap05","chr2_hap02","chr2_hap03","chr2_hap04","chr2_hap06","chr2_hap07")
length <- c(seqs$length[1],seqs$length[5], seqs$length[2], seqs$length[3], seqs$length[4], seqs$length[6], seqs$length[7])
bins <- data.frame(bin_id, seq_id, length)

links_700bp <- filter(links, map_length >= 700)

AvrSr27_plot <- gggenomes(
  genes = gff,
  seqs = bins,
  links = links_700bp,
  .id = "file_id",
  spacing = 0.05,
  wrap = NULL,
  adjacent_only = TRUE, #will also plot links b/w non-adjacent seqs
  infer_bin_id = seq_id,
  infer_start = min(start, end),
  infer_end = max(start, end),
  infer_length = max(start, end),
  theme = c("clean", NULL),
  .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_link(color = "#d6f0eb", fill = "#d6f0eb") +
  geom_gene(aes(fill=name), color = "black", size = 3, show.legend = FALSE) + geom_gene_label(aes(label = name), size = 2, angle = 0, hjust = -0.05, vjust = -0.05, fontface = "italic") + 
  geom_bin_label(size = 2.5, hjust = 0.2) + scale_fill_manual(values = palette[1:nrow(gff)])

#AvrSr27_plot

ggsave(file = "AvrSr27_synteny_plot.tiff",
       AvrSr27_plot, device = "tiff", dpi = 600, height = 10, width = 18, units = "cm")



## AvrSr50/35

seqs <- read_seqs("AvrSr50_35_intervals.fasta")
links <- read_paf("AvrSr50_35_intervals.paf")
gff <- read_gff3("AvrSr50_35.gff")
interval_coords <- read.table("AvrSr50_35_intervals_for_gggenomes.txt", header = FALSE, sep="\t")
colnames(interval_coords) <- c("seq_id","fullstart","fullend")


info <- left_join(gff, interval_coords, by = "seq_id")
info$newstart <- info$start - info$fullstart
info$newend <- info$end - info$fullstart

gff$start <- info$newstart
gff$end <- info$newend

bin_id <- c("hap01","hap04","hap05","hap06","hap07","hap02","hap03")
seq_id <- c("chr14_hap01","chr14_hap04","chr14_hap05","chr14_hap06","chr14_hap07","chr14_hap02","chr14_hap03")
length <- c(seqs$length[1],seqs$length[4], seqs$length[5], seqs$length[6], seqs$length[7], seqs$length[2], seqs$length[3])
bins <- data.frame(bin_id, seq_id, length)

links_500bp <- filter(links, map_length >= 500)

AvrSr50_35_plot <- gggenomes(
  genes = gff,
  seqs = bins,
  links = links_500bp,
  .id = "file_id",
  spacing = 0.05,
  wrap = NULL,
  adjacent_only = TRUE, #will also plot links b/w non-adjacent seqs
  infer_bin_id = seq_id,
  infer_start = min(start, end),
  infer_end = max(start, end),
  infer_length = max(start, end),
  theme = c("clean", NULL),
  .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_link(color = "#d6f0eb", fill = "#d6f0eb") +
  geom_gene(aes(fill=name), color = "black", size = 2, show.legend = FALSE) + geom_gene_label(aes(label = name), size = 2.5, angle = 0, hjust = -0.05, vjust = -0.15, fontface = "italic") + 
  geom_bin_label(size = 3, hjust = 0.4) + scale_fill_manual(values = c("#F8766D","#d3ff00","#F8766D","#d3ff00","#7CAE00","#31e08f","#00BFC4","grey40"), 
                                                            breaks = c("AvrSr35-01","AvrSr35-02","AvrSr50-01","AvrSr50-02","AvrSr50-03","AvrSr50-04","AvrSr50-05","AvrSr35-MITE"))

#AvrSr50_35_plot

ggsave(file = "AvrSr50_35_synteny_plot.tiff",
       AvrSr50_35_plot, device = "tiff", dpi = 600, height = 10, width = 18, units = "cm")



## AvrSr62

seqs <- read_seqs("AvrSr62_intervals.fasta")
links <- read_paf("AvrSr62_intervals.paf")
gff <- read_gff3("AvrSr62.gff")
interval_coords <- read.table("AvrSr62_intervals_for_gggenomes.txt", header = FALSE, sep="\t")
colnames(interval_coords) <- c("seq_id","fullstart","fullend")


info <- left_join(gff, interval_coords, by = "seq_id")
info$newstart <- info$start - info$fullstart
info$newend <- info$end - info$fullstart

gff$start <- info$newstart
gff$end <- info$newend

bin_id <- c("hap01","hap05","hap02","hap03","hap04","hap06","hap07")
seq_id <- c("chr5_hap01","chr5_hap05","chr5_hap02","chr5_hap03","chr5_hap04","chr5_hap06","chr5_hap07")
length <- c(seqs$length[1],seqs$length[5], seqs$length[2], seqs$length[3], seqs$length[4], seqs$length[6], seqs$length[7])
bins <- data.frame(bin_id, seq_id, length)

links_2000bp <- filter(links, map_length >= 2000)

AvrSr62_plot <- gggenomes(
  genes = gff,
  seqs = bins,
  links = links_2000bp,
  .id = "file_id",
  spacing = 0.05,
  wrap = NULL,
  adjacent_only = TRUE, #will also plot links b/w non-adjacent seqs
  infer_bin_id = seq_id,
  infer_start = min(start, end),
  infer_end = max(start, end),
  infer_length = max(start, end),
  theme = c("clean", NULL),
  .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_link(color = "#d6f0eb", fill = "#d6f0eb") +
  geom_gene(aes(fill=name), color = "black", size = 2, show.legend = FALSE) + geom_gene_label(aes(label = name), size = 2.5, angle = 30, hjust = -0.05, vjust = -0.05, fontface = "italic") + 
  geom_bin_label(size = 3, hjust = 0.4) + scale_fill_manual(values = c("#F8766D","#d3ff00","#7CAE00","#31e08f","#00BFC4","#0031ff","#C77CFF","#ff00b1","#00ffce","grey40","#000000","white"), 
                                                            breaks = c("AvrSr62-01","AvrSr62-02","AvrSr62-03","AvrSr62-04","AvrSr62-05","AvrSr62-06","AvrSr62-07","AvrSr62-08","AvrSr62-09","AvrSr62-10","AvrSr62-11","AvrSr62-pseudo"))

#AvrSr62_plot

ggsave(file = "AvrSr62_synteny_plot.tiff",
       AvrSr62_plot, device = "tiff", dpi = 600, height = 12, width = 19, units = "cm")



## AvrSr33

seqs <- read_seqs("AvrSr33_intervals.fasta")
links <- read_paf("AvrSr33_intervals.paf")
gff <- read_gff3("AvrSr33.gff")
interval_coords <- read.table("AvrSr33_intervals_for_gggenomes.txt", header = FALSE, sep="\t")
colnames(interval_coords) <- c("seq_id","fullstart","fullend")


info <- left_join(gff, interval_coords, by = "seq_id")
info$newstart <- info$start - info$fullstart
info$newend <- info$end - info$fullstart

gff$start <- info$newstart
gff$end <- info$newend

bin_id <- c("hap04","hap07","hap01","hap02","hap05","hap06","hap03")
seq_id <- c("chr17_hap04","chr17_hap07","chr17_hap01","chr17_hap02","chr17_hap05","chr17_hap06","chr17_hap03")
length <- c(seqs$length[4], seqs$length[7], seqs$length[1],seqs$length[2], seqs$length[5], seqs$length[6], seqs$length[3])
bins <- data.frame(bin_id, seq_id, length)

links_200bp <- filter(links, map_length >= 200)

AvrSr33_plot <- gggenomes(
  genes = gff,
  seqs = bins,
  links = links_200bp,
  .id = "file_id",
  spacing = 0.05,
  wrap = NULL,
  adjacent_only = TRUE, #will also plot links b/w non-adjacent seqs
  infer_bin_id = seq_id,
  infer_start = min(start, end),
  infer_end = max(start, end),
  infer_length = max(start, end),
  theme = c("clean", NULL),
  .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_link(color = "#d6f0eb", fill = "#d6f0eb") +
  geom_gene(aes(fill=name), color = "black", size = 2, show.legend = FALSE) + scale_fill_manual(values = c("#F8766D","#d3ff00","#7CAE00","#31e08f","#00BFC4","#0031ff","#C77CFF","#ff00b1","grey50"), 
                                                                                                breaks = c("AvrSr33-01","AvrSr33-02","AvrSr33-03","AvrSr33-04","AvrSr33-05","AvrSr33-06","AvrSr33-07","AvrSr33-08","AvrSr33-pseudo")) + 
  geom_gene_label(aes(label = name), size = 2.5, angle = 15, hjust = -0.05, vjust = -0.15, fontface = "italic") + 
  geom_bin_label(size = 3, hjust = 0.4)

#AvrSr33_plot

ggsave(file = "AvrSr33_synteny_plot.tiff",
       AvrSr33_plot, device = "tiff", dpi = 600, height = 10, width = 18, units = "cm")
