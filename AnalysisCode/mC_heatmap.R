library(dplyr)
library(RColorBrewer)
library(readr)
library(ggplot2)

mc_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC_no_insertion.tsv.gz", col_names = T)
mc_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC_true_deletions.tsv.gz", col_names = T)

# sum rows for sort order
mc_no_insertion_2kb.sums <- mutate(mc_no_insertion_2kb, sums=rowSums(mc_no_insertion_2kb[,2:41], na.rm = T))
mc_true_deletion_2kb.sums <- mutate(mc_true_deletion_2kb, sums=rowSums(mc_true_deletion_2kb[,2:41], na.rm = T))

# Sort by total mC level
mc_no_insertion_2kb.sort <- arrange(mc_no_insertion_2kb.sums, sums) %>% select(-sums)
mc_true_deletion_2kb.sort <- arrange(mc_true_deletion_2kb.sums, sums) %>% select(-sums)

scale_max <- function(m, l){
  m.floor <- m > l
  m[m.floor] <- l
  return(t(data.matrix(m)))
}
both <- rbind(mc_true_deletion_2kb.sort, mc_no_insertion_2kb.sort)
ins <- scale_max(mc_no_insertion_2kb.sort[,2:40], 0.5)
del <- scale_max(mc_true_deletion_2kb.sort[,2:40], 0.5)
both_scaled <- scale_max(both[,2:40], 0.5)
color <- colorRampPalette(brewer.pal(9,"Reds"))(100)

png("../Plots/heatmap_insertion.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(ins, col = color, main = "TE insertions")
dev.off()

png("../Plots/heatmap_deletion.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(del, col = color, main = "TE deletions")
dev.off()

png("../Plots/heatmap_ins_del.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(both_scaled, col=color)
dev.off()

# scale bar
pdf("../Plots/mc_heatmap_scale.pdf", height=2, width=5)
image(data.matrix(seq(100)), col = color)
dev.off()


# make line charts of each mC context

mCG_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CG_no_insertion.tsv.gz", col_names = T)
mCG_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CG_true_deletions.tsv.gz", col_names = T)

mCHG_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHG_no_insertion.tsv.gz", col_names = T)
mCHG_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHG_true_deletions.tsv.gz", col_names = T)

mCHH_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHH_no_insertion.tsv.gz", col_names = T)
mCHH_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHH_true_deletions.tsv.gz", col_names = T)


# get colMeans
ins_means_cg <- colMeans(mCG_no_insertion_2kb[,2:21], na.rm = T)
ins_means_cg_absent <- colMeans(mCG_no_insertion_2kb[,21:41], na.rm = T)

ins_means_chg <- colMeans(mCHG_no_insertion_2kb[,2:21], na.rm = T)
ins_means_chg_absent <- colMeans(mCHG_no_insertion_2kb[,21:41], na.rm = T)

ins_means_chh <- colMeans(mCHH_no_insertion_2kb[,2:21], na.rm = T)
ins_means_chh_absent <- colMeans(mCHH_no_insertion_2kb[,21:41], na.rm = T)

del_means_cg <- colMeans(mCG_true_deletion_2kb[,2:21], na.rm = T)
del_means_cg_absent <- colMeans(mCG_true_deletion_2kb[,21:41], na.rm = T)

del_means_chg <- colMeans(mCHG_true_deletion_2kb[,2:21], na.rm = T)
del_means_chg_absent <- colMeans(mCHG_true_deletion_2kb[,21:41], na.rm = T)

del_means_chh <- colMeans(mCHH_true_deletion_2kb[,2:21], na.rm = T)
del_means_chh_absent <- colMeans(mCHH_true_deletion_2kb[,21:41], na.rm = T)

# define colours
cg <- '#B4B464'
chg <- '#6665AD'
chh <- '#B29492'

line_plot <- function(d1, d2, d3, scale) {
  plot(d1,
       type = 'l',
       main = "TE insertion",
       ylab = "Methylation level",
       xlab = "Position",
       las = 1,
       lwd = 2,
       cex.main=0.6,
       cex.axis=0.6,
       cex.lab=0.6,
       col = cg,
       ylim = c(0,scale)
  )
  lines(d2,
        col=chg,
        lwd = 2
  )
  lines(d3,
        col=chh,
        lwd = 2
  )
}

pdf("../Plots/flanking_mc_line_charts.pdf", height=4, width=4, useDingbats = F)
par(mfrow=c(2,2), mar=c(2,2,2,2))
line_plot(ins_means_cg,ins_means_chg,ins_means_chh, 0.8)
line_plot(ins_means_cg_absent, ins_means_chg_absent, ins_means_chh_absent, 0.8)
line_plot(del_means_cg,del_means_chg,del_means_chh, 0.8)
line_plot(del_means_cg_absent, del_means_chg_absent, del_means_chh_absent, 0.8)
dev.off()

# correlation between flanking DNA methylation levels and TE presence/absence
flanking_mc_corr <- read_tsv("../ProcessedData/DNAmethylation/flanking_mc_correlations.txt") %>%
  mutate(id = paste(chr, start, stop)) %>%
  select(-(chr:stop))

tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz") %>%
  mutate(id = paste(chromosome, start, end))

data <- left_join(tepav, flanking_mc_corr, by="id") %>%
  select(-id)

data$LD <- factor(data$LD, levels = c("high", "mid", "low", NA))

ggplot(data, aes(Absence_classification, cor, fill=LD)) +
  geom_boxplot(outlier.shape = 1) + theme_bw() +
  scale_fill_brewer(type="seq", palette = 11, direction = -1) +
  geom_hline(yintercept = 0) +
  coord_flip() + ylab("Correlation") + xlab("TE variant") +
  ggsave("../Plots/mc_flanking_correlation.pdf", height = 8, width = 15, units = "cm")

ggplot(data, aes(MAF, cor, fill=LD)) +
  geom_boxplot(aes(group = cut_width(MAF, 0.04)), outlier.shape = 1, outlier.size = 0.8) +
  theme_bw() + geom_hline(yintercept = 0) + facet_wrap(~LD, ncol = 1) +
  scale_fill_brewer(type="seq", palette = 11, direction = -1) +
  theme(text = element_text(size=8)) +
  ggsave("../Plots/mc_flanking_correlation_by_ld.pdf", height = 20, width = 8, units = "cm", useDingbats=F)

# Bud vs lead methylation

bud_leaf_ins <- read_tsv("../ProcessedData/DNAmethylation/bud_vs_leaf_allC_insertions.tsv.gz", col_names = T, na = c("NaN"))
bud_leaf_del <- read_tsv("../ProcessedData/DNAmethylation/bud_vs_leaf_allC_deletions.tsv.gz", col_names = T, na = c("NaN"))

bud_leaf_ins.sums <- mutate(bud_leaf_ins, sums=rowSums(bud_leaf_ins[,2:ncol(bud_leaf_ins)], na.rm = T))
bud_leaf_del.sums <- mutate(bud_leaf_del, sums=rowSums(bud_leaf_del[,2:ncol(bud_leaf_del)], na.rm = T))

bud_leaf_ins.sort <- arrange(bud_leaf_ins.sums, sums) %>% select(-sums)
bud_leaf_del.sort <- arrange(bud_leaf_del.sums, sums) %>% select(-sums)

both_bud <- rbind(bud_leaf_ins.sort, bud_leaf_del.sort)

bl <- scale_max(both_bud[,2:ncol(both_bud)], 0.5)

png("../Plots/heatmap_bud_leaf_indel.png", height = 6, width = 3, units = "in", res = 600, bg = "transparent")
image(bl, col = color)
dev.off()
