library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(grid)

pvals <- read_tsv("../ProcessedData/cdmr_sig_pvals.tsv", col_names = c("dmr_chr", "dmr_start", "te_chr", "te_start", "indel"))

pvals$te_chr <- factor(pvals$te_chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5"))	
pvals$dmr_chr <- factor(pvals$dmr_chr, levels = c("chr5", "chr4", "chr3", "chr2", "chr1"))	

insertion_col <- brewer.pal(3, "Set2")[1]
deletion_col <- brewer.pal(3, "Set2")[3]

pvals %>%
 filter(indel=="Deletion") %>%
ggplot(., aes(te_start, dmr_start)) +
  geom_point(size=0.01, alpha=0.1, color=deletion_col) +
  facet_wrap(dmr_chr~te_chr, scales = "free") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.margin = unit(0, "lines"))

ggsave("../Plots/pval_del.png", height=5, width = 5, dpi = 600)

pvals %>%
  filter(indel=="Insertion") %>%
  ggplot(., aes(te_start, dmr_start)) +
  geom_point(size=0.01, alpha=0.1, color=insertion_col) +
  facet_wrap(dmr_chr~te_chr, scales = "free") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.margin = unit(0, "lines"))

ggsave("../Plots/pval_ins.png", height=5, width = 5, dpi = 600)

dot_plot <- function(pval, chr, y, z) {
  filtered <- subset(pval, te_chr == chr)
  h <- hist(filtered$te_start, breaks = 10000, plot=F)
  d <- data.frame(x = h$breaks[1:length(h$breaks) - 1], y = h$counts, z = NA)
  d <- d[d$y > 0,]
  # add column with color
  d[d$y > z, 3] <- 1
  d[d$y <= z, 3] <- 2
  plot(d$x, d$y, pch=19, cex=0.5, ylim = c(0, y), col=c("red", "grey")[d$z])
  abline(h=z, col="grey")
}

ins <- filter(pvals, indel == "Insertion")
del <- filter(pvals, indel == "Deletion")

pdf("../Plots/trans_dotplot_insertions.pdf", height=1.5, width=10, useDingbats = F)
par(mfrow=c(1,5), mar = c(1,1,1,1))
dot_plot(ins, "chr1", 400, 200)
dot_plot(ins, "chr2", 400, 200)
dot_plot(ins, "chr3", 400, 200)
dot_plot(ins, "chr4", 400, 200)
dot_plot(ins, "chr5", 400, 200)
dev.off()

pdf("../Plots/trans_dotplot_deletions.pdf", height=1.5, width=10, useDingbats = F)
par(mfrow=c(1,5), mar = c(1,1,1,1))
dot_plot(del, "chr1", 700, 200)
dot_plot(del, "chr2", 700, 200)
dot_plot(del, "chr3", 700, 200)
dot_plot(del, "chr4", 700, 200)
dot_plot(del, "chr5", 700, 200)
dev.off()

# draw heatmaps of TE density to go along the edges of the plots
te <- read_tsv("../RawData/TAIR9_TE.bed.gz", col_names = c("chr", "start", "stop", "strand", "TE", "Family", "Superfamily"))

color <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

plot_heatmap <- function(d, chr, color) {
  c <- d[d$chr == chr,]
  h <- hist(c$start, breaks = 200, plot = F)
  image(as.matrix(h$counts), col = color)
}

pdf("../Plots/te_density_heatmaps.pdf", height = 10, width = 5)
par(mfrow = c(5,1))
plot_heatmap(te, "chr1", color)
plot_heatmap(te, "chr2", color)
plot_heatmap(te, "chr3", color)
plot_heatmap(te, "chr4", color)
plot_heatmap(te, "chr5", color)
dev.off()
