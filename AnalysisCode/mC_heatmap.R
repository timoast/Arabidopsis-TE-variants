library(dplyr)
library(RColorBrewer)
library(readr)
library(ggplot2)

mc_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC_no_insertion.tsv.gz",
                                col_names = T, na = c("nan")) %>%
  mutate(class = "Insertion")

mc_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC_true_deletions.tsv.gz",
                                 col_names = T, na = c("nan")) %>%
  mutate(class = "Deletion")

# mc_all <- read_tsv("../ProcessedData/DNAmethylation/flanking_allC.tsv.gz", col_names = T, na = c("nan"))

all <- rbind(mc_true_deletion_2kb, mc_no_insertion_2kb)

# Add column with distance to centromere
get_centro_distance <- function(crds, centromeres) {
  crds <- unlist(strsplit(crds, ","))
  cent <- filter(centromeres, chrom == crds[[1]])
  midpoint <- (cent[[3]] + cent[[2]]) / 2
  d <- abs(as.numeric(crds[[2]]) - midpoint)
  return(d)
}

centromeres <- read_tsv("../RawData/centromere_positions.txt", col_names = c("chrom", "start", "stop"))

# separate insertions into pericentromeric / euchromatic
# sort within by mC
# make line charts separately for each

# sum rows for sort order
all.sort <- all %>%
  mutate(mn = rowMeans(all[,2:41], na.rm = T)) %>%
  rowwise() %>%
  mutate(distance = get_centro_distance(coords, centromeres),
         pericentromeric = distance < 3 * 10^6) %>%
  arrange(class, pericentromeric, mn)

dels <- filter(all.sort, class == "Deletion")
ins <- filter(all.sort, class == "Insertion")

distance_dels <- select(dels, distance)
distance_ins <- select(ins, distance)

scale_max <- function(m, l){
  m.floor <- m > l
  m[m.floor] <- l
  return(t(data.matrix(m)))
}


both_scaled <- scale_max(all.sort[,2:40], 0.5)
ins_scaled <- scale_max(ins[,2:40], 0.5)
del_scaled <- scale_max(dels[,2:40], 0.5)

color <- colorRampPalette(brewer.pal(9,"Reds"))(100)
d_color <- colorRampPalette(brewer.pal(9, "YlOrBr"))(100)

png("../Plots/heatmap_mc.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(both_scaled, col = color)
dev.off()

png("../Plots/heatmap_insertion.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(ins_scaled, col = color, main = "TE insertions")
dev.off()

png("../Plots/heatmap_deletion.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(del_scaled, col = color, main = "TE deletions")
dev.off()

png("../Plots/heatmap_ins_del.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(both_scaled, col=color)
dev.off()

png("../Plots/heatmap_centromere_distances_deletions.png", height = 6, width = 3, units = "in", res = 600)
image(scale_max(distance_dels, 15 * 10^6), col = d_color)
dev.off()

png("../Plots/heatmap_centromere_distances_insertion.png", height = 6, width = 3, units = "in", res = 600)
image(scale_max(distance_ins, 15 * 10^6), col = d_color)
dev.off()

# scale bar
pdf("../Plots/mc_heatmap_scale.pdf", height=2, width=5)
image(data.matrix(seq(100)), col = color)
dev.off()

pdf("../Plots/scale_ylorbr.pdf", height=2, width=5)
image(data.matrix(seq(100)), col = d_color)
dev.off()

# make line charts of each mC context
# Load context-specific methylation data
# join with pericentromeric T/F and methylation sum column from allC context data and sort to get same sort order
# Separate by pericentromeric T/F and plot line chart of column means

ins_data <- select(ins, coords, mn, pericentromeric)
del_data <- select(dels, coords, mn, pericentromeric)

mCG_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CG_no_insertion.tsv.gz", col_names = T, na = c("nan")) %>%
  left_join(., ins_data, by="coords") %>%
  arrange(pericentromeric, mn)

mCG_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CG_true_deletions.tsv.gz", col_names = T, na = c("nan")) %>%
  left_join(., del_data, by="coords") %>%
  arrange(pericentromeric, mn)

mCHG_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHG_no_insertion.tsv.gz", col_names = T, na = c("nan")) %>%
  left_join(., ins_data, by="coords") %>%
  arrange(pericentromeric, mn)

mCHG_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHG_true_deletions.tsv.gz", col_names = T, na = c("nan")) %>%
  left_join(., del_data, by="coords") %>%
  arrange(pericentromeric, mn)

mCHH_no_insertion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHH_no_insertion.tsv.gz", col_names = T, na = c("nan")) %>%
  left_join(., ins_data, by="coords") %>%
  arrange(pericentromeric, mn)

mCHH_true_deletion_2kb <- read_tsv("../ProcessedData/DNAmethylation/flanking_CHH_true_deletions.tsv.gz", col_names = T, na = c("nan")) %>%
  left_join(., del_data, by="coords") %>%
  arrange(pericentromeric, mn)

# get colMeans
ins_means_cg_eu <- colMeans(filter(mCG_no_insertion_2kb, pericentromeric == F)[,2:21], na.rm = T)
ins_means_cg_het <- colMeans(filter(mCG_no_insertion_2kb, pericentromeric == T)[,2:21], na.rm = T)

ins_means_cg_absent_eu <- colMeans(filter(mCG_no_insertion_2kb, pericentromeric == F)[,21:41], na.rm = T)
ins_means_cg_absent_het <- colMeans(filter(mCG_no_insertion_2kb, pericentromeric == T)[,21:41], na.rm = T)

ins_means_chg_eu <- colMeans(filter(mCHG_no_insertion_2kb, pericentromeric == F)[,2:21], na.rm = T)
ins_means_chg_het <- colMeans(filter(mCHG_no_insertion_2kb, pericentromeric == T)[,2:21], na.rm = T)

ins_means_chg_absent_eu <- colMeans(filter(mCHG_no_insertion_2kb, pericentromeric == F)[,21:41], na.rm = T)
ins_means_chg_absent_het <- colMeans(filter(mCHG_no_insertion_2kb, pericentromeric == T)[,21:41], na.rm = T)

ins_means_chh_eu <- colMeans(filter(mCHH_no_insertion_2kb, pericentromeric == F)[,2:21], na.rm = T)
ins_means_chh_het <- colMeans(filter(mCHH_no_insertion_2kb, pericentromeric == T)[,2:21], na.rm = T)

ins_means_chh_absent_eu <- colMeans(filter(mCHH_no_insertion_2kb, pericentromeric == F)[,21:41], na.rm = T)
ins_means_chh_absent_het <- colMeans(filter(mCHH_no_insertion_2kb, pericentromeric == T)[,21:41], na.rm = T)

del_means_cg_eu <- colMeans(filter(mCG_true_deletion_2kb, pericentromeric == F)[,2:21], na.rm = T)
del_means_cg_het <- colMeans(filter(mCG_true_deletion_2kb, pericentromeric == T)[,2:21], na.rm = T)

del_means_cg_absent_eu <- colMeans(filter(mCG_true_deletion_2kb, pericentromeric == F)[,21:41], na.rm = T)
del_means_cg_absent_het <- colMeans(filter(mCG_true_deletion_2kb, pericentromeric == T)[,21:41], na.rm = T)

del_means_chg_eu <- colMeans(filter(mCHG_true_deletion_2kb, pericentromeric == F)[,2:21], na.rm = T)
del_means_chg_het <- colMeans(filter(mCHG_true_deletion_2kb, pericentromeric == T)[,2:21], na.rm = T)

del_means_chg_absent_eu <- colMeans(filter(mCHG_true_deletion_2kb, pericentromeric == F)[,21:41], na.rm = T)
del_means_chg_absent_het <- colMeans(filter(mCHG_true_deletion_2kb, pericentromeric == T)[,21:41], na.rm = T)

del_means_chh_eu <- colMeans(filter(mCHH_true_deletion_2kb, pericentromeric == F)[,2:21], na.rm = T)
del_means_chh_het <- colMeans(filter(mCHH_true_deletion_2kb, pericentromeric == T)[,2:21], na.rm = T)

del_means_chh_absent_eu <- colMeans(filter(mCHH_true_deletion_2kb, pericentromeric == F)[,21:41], na.rm = T)
del_means_chh_absent_het <- colMeans(filter(mCHH_true_deletion_2kb, pericentromeric == T)[,21:41], na.rm = T)

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

pdf("../Plots/flanking_mc_line_charts.pdf", height=5, width=3, useDingbats = F)
par(mfrow=c(4, 2), mar=c(2,2,2,2), cex = 0.8)
line_plot(ins_means_cg_het,ins_means_chg_het,ins_means_chh_het, 1)
line_plot(ins_means_cg_absent_het,ins_means_chg_absent_het,ins_means_chh_absent_het, 1)

line_plot(ins_means_cg_eu,ins_means_chg_eu,ins_means_chh_eu, 0.5)
line_plot(ins_means_cg_absent_eu, ins_means_chg_absent_eu, ins_means_chh_absent_eu, 0.5)

line_plot(del_means_cg_het, del_means_chg_het, del_means_chh_het , 1)
line_plot(del_means_cg_absent_het, del_means_chg_absent_het, del_means_chh_absent_het, 1)

line_plot(del_means_cg_eu, del_means_chg_eu, del_means_chh_eu , 0.5)
line_plot(del_means_cg_absent_eu, del_means_chg_absent_eu, del_means_chh_absent_eu, 0.5)
dev.off()

# Bud vs lead methylation

bud_leaf_ins <- read_tsv("../ProcessedData/DNAmethylation/bud_vs_leaf_allC_insertions.tsv.gz", col_names = T, na = c("nan")) %>%
  mutate(class = "Insertion")
bud_leaf_del <- read_tsv("../ProcessedData/DNAmethylation/bud_vs_leaf_allC_deletions.tsv.gz", col_names = T, na = c("nan")) %>%
  mutate(class = "Deletion")

# bud_leaf_ins.sums <- mutate(bud_leaf_ins, sums=rowSums(bud_leaf_ins[,2:ncol(bud_leaf_ins)], na.rm = T))
# bud_leaf_del.sums <- mutate(bud_leaf_del, sums=rowSums(bud_leaf_del[,2:ncol(bud_leaf_del)], na.rm = T))

bud_all <- rbind(bud_leaf_del, bud_leaf_ins)

# sum rows for sort order
bud_all_sort <- bud_all %>%
  mutate(mn = rowMeans(bud_all[,2:41], na.rm = T)) %>%
  rowwise() %>%
  mutate(distance = get_centro_distance(coords, centromeres),
         pericentromeric = distance < 3 * 10^6) %>%
  arrange(class, pericentromeric, mn)

dels_bud <- filter(bud_all_sort, class == "Deletion")
ins_bud <- filter(bud_all_sort, class == "Insertion")

distance_dels <- select(dels_bud, distance)
distance_ins <- select(ins_bud, distance)

bud_all_sort <- select(bud_all_sort, -(class:pericentromeric))
bl <- scale_max(bud_all_sort[,2:ncol(bud_all_sort)], 0.5)

png("../Plots/heatmap_bud_leaf_indel.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(bl, col = color)
dev.off()

png("../Plots/heatmap_bud_leaf_indel_distance_dels.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(scale_max(distance_dels, 15 * 10^6), col = d_color)
dev.off()

png("../Plots/heatmap_bud_leaf_indel_distance_ins.png", height = 6, width = 3, units = "in", res = 600, bg = "grey")
image(scale_max(distance_ins, 15 * 10^6), col = d_color)
dev.off()
