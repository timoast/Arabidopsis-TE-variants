library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# TEPID TE calls
tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = T)

# 1. Plot minor allele frequency
pdf("../Plots/MAF.pdf", height = 4, width = 4)
hist(tepav$MAF*100, breaks = 100, col = "grey", xlab = "MAF", main = "Minor allele frequency distribution")
dev.off()

# 2. Plot TE variant length distribution
all_te <- read_tsv("../RawData/TAIR9_TE.bed.gz", col_names = F) %>%
  mutate(l = X3-X2)

tepav_tes <- unlist(strsplit(tepav$TE, ","))
lengths <- all_te %>%
  mutate(pav = X5 %in% tepav_tes) %>%
  filter(pav == T)

pdf(file = "../Plots/length_distribution_histogram.pdf", width = 4, height = 4)
hist(all_te$l, breaks = 100, col = "grey", main = "Length distribution of all Col-0 TEs", xlab = "Length (bp)")
rug(all_te$l)

hist(lengths$l, breaks = 100, col = "grey", main = "Length distribution of population-variable TEs", xlab = "Length (bp")
rug(lengths$l)
dev.off()

# plot logs
pdf("../Plots/TE_length_log.pdf", width = 5, height = 5)
plot(density(log(lengths$l)), col=rgb(1, 0, 0,0.5), main = "TE length distribution", xlab = "log10 TE length")
polygon(density(log(lengths$l)), col=rgb(1, 0, 0,0.5))
polygon(density(log(all_te$l)), col = rgb(0, 0, 1,0.5))
legend(8,0.2,legend = c("All", "TEPAV"), lty=c(1,1), col=c("blue", "red"), bty = "n")
dev.off()

# 3. Plot genomic distribution of TE variants, for rare, common, TE deletions vs no insertions, and all together

rare_vs_common_colors <- c("#FC9272", "#67000D")
ins_vs_del_colors <- c("#9ECAE1", "#08306B")

ggplot(tepav, aes(start, color = Frequency_classification)) +
  geom_density(adjust = 1/8) + scale_color_manual(values = rare_vs_common_colors) +
  facet_wrap(~chromosome, nrow = 1) + theme_bw() +
  ggsave("chomosomal_dist_rare_vs_common_line_chart.pdf", path = "../Plots", height = 3, width = 20)

ggplot(tepav, aes(start, color = Absence_classification)) +
  geom_density(adjust = 1/8) + scale_color_manual(values = ins_vs_del_colors) +
  facet_wrap(~chromosome, nrow = 1) + theme_bw() +
  ggsave("chomosomal_dist_ins_vs_del_line_chart.pdf", path = "../Plots", height = 3, width = 20)

# show as heatmap

get_distribution <- function(d, n) {
  h <- hist(d, breaks = n)
  o <- h$counts / max(h$counts)
  return(o)
}

chr1 <- filter(tepav, chromosome == "chr1")
all_tepav <- get_distribution(chr1$start, 1000)[1:608]
rare_te <- get_distribution(filter(chr1, Frequency_classification == "Rare")$start, 1000)[1:608]
common_te <- get_distribution(filter(chr1,  Frequency_classification == "Common")$start, 1000)[1:608]
te_ins <- get_distribution(filter(chr1, Absence_classification == "No insertion")$start, 1000)[1:608]
te_del <- get_distribution(filter(chr1, Absence_classification == "True deletion")$start, 1000)[1:608]

# # group together and scale each to it's maximum value
# all_hist_scaled <- cbind(all_tepav, rare_te, common_te, te_ins, te_del)
# 
# color <- (colorRampPalette(brewer.pal(9,"GnBu"))(100))
# pdf("../Plots/heatmap_te_insertions_chr1.pdf", height = 3, width = 5)
# image(all_hist_scaled, col = color, xlab = "Chromosome 1", ylab = c("TE variants"))
# dev.off()

# Plot side-by-side with all Col-0 TEs, C-DMRs, CG-DMRs, genes
all_te_chr1 <- get_distribution(filter(all_te, X1 == "chr1")$X2, 500)[1:608]
all_genes <- read_tsv("../RawData/TAIR10_GFF3_genes.gff.gz", col_names = F)
all_genes_chr1 <- get_distribution(filter(all_genes, X1 == "chr1")$X4, 500)[1:608]
all_c_dmrs <- read_tsv("../RawData/c_dmrs.tsv.gz", col_names = F)
c_dmrs_chr1 <- get_distribution(filter(all_c_dmrs, X1 == "chr1")$X2, 500)[1:608]
all_cg_dmrs <- read_tsv("../RawData/cg_dmrs.tsv.gz", col_names = F)
cg_dmrs_chr1 <- get_distribution(filter(all_cg_dmrs, X1 == "chr1")$X2, 500)[1:608]

all_features_chr1 <- cbind(te_ins, te_del, common_te, rare_te, c_dmrs_chr1, all_tepav, all_te_chr1,
                           cg_dmrs_chr1, all_genes_chr1)

pdf("../Plots/heatmap_te_insertions_chr1_all_features.pdf", height = 4, width = 5)
image(all_features_chr1, col = color, xlab = "Chromosome 1")
dev.off()

# order top to bottom:
# Col-0 genes
# CG-DMRS
# Col-0 TEs
# TEPAVs
# C-DMRs
# Rare TEPAVs
# Common TEPAVs
# True TE deletions
# No insertion variants

# scale bar

pdf("../Plots/chr1_heatmap_scale.pdf", height=2, width=5)
image(data.matrix(seq(100)), col = color)
dev.off()
