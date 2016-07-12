library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

te_c_dmr <- read_tsv("../ProcessedData/TE_C_DMR_distances.bed.gz", col_names = F) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15, LD=X16,
         distance = X17)

te_cg_dmr <- read_tsv("../ProcessedData/TE_CG_DMR_distances.bed.gz", col_names = F) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15, LD=X16,
         distance = X17)

random_c_dmr <- read_tsv("../ProcessedData/random_coords_closest_c_dmr.bed.gz", col_names = F)
random_cg_dmr <- read_tsv("../ProcessedData/random_coords_closest_cg_dmr.bed.gz", col_names = F)

cg <- "#B4B464"
chg <- "#6665AD"
chh <- "#B29492"

rare_col <- brewer.pal(6, "Set3")[4]
common_col <- brewer.pal(6, "Set3")[6]

insertion_col <- brewer.pal(3, "Set2")[1]
deletion_col <- brewer.pal(3, "Set2")[3]
na_col <- brewer.pal(9, "Blues")[2]

pdf("../Plots/distances_TEPAV_DMR.pdf", height = 3, width = 3)
hist(te_c_dmr$distance / 1000,
     breaks = 200, xlim = c(0, 20),
     col = chh, xlab = "Distance (kb)", main = "Distance from TE variant to nearest C-DMR")
hist(te_cg_dmr$distance / 1000,
     breaks = 300, xlim = c(0, 20),
     col = cg, xlab = "Distance (kb)", main = "Distance from TE variant to nearest CG-DMR")
dev.off()

te_c_dmr$AbsenceClassification <- factor(te_c_dmr$AbsenceClassification, levels = c("True deletion", "No insertion"))	
te_cg_dmr$AbsenceClassification <- factor(te_cg_dmr$AbsenceClassification, levels = c("True deletion", "No insertion"))	

# C-DMRs
ggplot(te_c_dmr, aes(distance / 1000, fill = FrequencyClassification)) +
  geom_density(alpha = 0.8) + facet_wrap(~FrequencyClassification) + xlim(0, 5) +
  scale_fill_manual(values = c(rare_col, common_col)) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  ggsave("../Plots/CDMR/c_dmr_distances_common_vs_rare.pdf", width = 10, height = 6, units = "cm", useDingbats=F)	

ggplot(na.omit(te_c_dmr), aes(distance / 1000, fill = AbsenceClassification)) +
  geom_density(alpha = 0.8) + facet_wrap(~AbsenceClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  scale_fill_manual(values =c(deletion_col, insertion_col)) + 
  ggsave("../Plots/CDMR/c_dmr_distances_insertion_vs_deletion.pdf", width = 10, height = 6, units = "cm", useDingbats=F)

# Percentage of insertion or deletions near C-DMRs
te_c_dmr %>%
  group_by(AbsenceClassification, distance < 1000) %>%
  summarise(count = n()) %>%
  mutate(tot = sum(count), perc = count / tot * 100)

# CG-DMRs
ggplot(te_cg_dmr, aes(distance / 1000, fill = FrequencyClassification)) +
  geom_density(alpha = 0.8) + facet_wrap(~FrequencyClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  scale_fill_manual(values = c(rare_col, common_col)) +
  ggsave("../Plots/CGDMR/cg_dmr_distances_common_vs_rare.pdf", width = 10, height = 6, units = "cm", useDingbats=F)	

ggplot(na.omit(te_cg_dmr), aes(distance / 1000, fill = AbsenceClassification)) +
  geom_density(alpha = 0.8) + facet_wrap(~AbsenceClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  scale_fill_manual(values =c(deletion_col, insertion_col)) + 
  ggsave("../Plots/CGDMR/cg_dmr_distances_insertion_vs_deletion.pdf", width = 10, height = 6, units = "cm", useDingbats=F)

# Percentage of insertion or deletions near C-DMRs
te_cg_dmr %>%
  group_by(AbsenceClassification, distance < 1000) %>%
  summarise(count = n()) %>%
  mutate(tot = sum(count), perc = count / tot * 100)

# percentage of DMRs with TE variant within 1 kb, compared to random regions replicated 10k times
# see ProcessingCode/bootstrap_te_variants.sh for source of random overlap percentages
c_dmrs <- read_tsv("../RawData/c_dmrs.tsv.gz", col_names = F)
cg_dmrs <- read_tsv("../RawData/cg_dmrs.tsv.gz", col_names = F)

bootstrap_c_dmr_insertions_1kb <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k_insertions.txt", col_names = F)
bootstrap_c_dmr_deletions_1kb <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k_deletions.txt", col_names = F)

bootstrap_cg_dmr_insertions_1kb <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k_insertions.txt", col_names = F)
bootstrap_cg_dmr_deletions_1kb <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k_deletions.txt", col_names = F)

counts <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = T) %>%
  mutate(l = end-start) %>%
  group_by(Absence_classification) %>%
  summarise(ml = mean(l), c = n())

pdf("../Plots/DMR_intersections_bootstrap.pdf", height = 5, width = 5)
par(mfrow = (c(2,2)))
hist(bootstrap_c_dmr_deletions_1kb$X1/counts$c[3] * 100,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "TE deletions near C-DMRs", border = chh)
intersections <- nrow(filter(te_c_dmr, AbsenceClassification == "True deletion" & distance < 1000))
total_class <- nrow(filter(te_c_dmr, AbsenceClassification == "True deletion"))
abline(v=intersections / total_class * 100)

hist(bootstrap_c_dmr_insertions_1kb$X1/counts$c[2] * 100,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "TE insertions near C-DMRs", border = chh)
intersections <- nrow(filter(te_c_dmr, AbsenceClassification == "No insertion" & distance < 1000))
total_class <- nrow(filter(te_c_dmr, AbsenceClassification == "No insertion"))
abline(v=intersections / total_class * 100)

hist(bootstrap_cg_dmr_deletions_1kb$X1/counts$c[3] * 100,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "TE deletions near CG-DMRs", border = cg)
intersections <- nrow(filter(te_cg_dmr, AbsenceClassification == "True deletion" & distance < 1000))
total_class <- nrow(filter(te_cg_dmr, AbsenceClassification == "True deletion"))
abline(v=intersections / total_class * 100)

hist(bootstrap_cg_dmr_insertions_1kb$X1/counts$c[2] * 100,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "TE insertions near CG-DMRs", border = cg)
intersections <- nrow(filter(te_cg_dmr, AbsenceClassification == "No insertion" & distance < 1000))
total_class <- nrow(filter(te_cg_dmr, AbsenceClassification == "No insertion"))
abline(v=intersections / total_class * 100)
dev.off()

# TE-DMR methylation level
c_dmr_mc <- read_tsv("../ProcessedData/c_dmr_allC.tsv.gz", col_names = T, na = c("None"))
cg_dmr_mc <- read_tsv("../ProcessedData/cg_dmrs_allC.tsv.gz", col_names = T, na = c("None"))

# Add common ID (DMR)
c_dmr_mc <-  c_dmr_mc %>%
  mutate(ID = paste(chr, start, stop)) %>%
  select(-(chr:stop))

te_c_dmr <-  te_c_dmr %>%
  mutate(ID = paste(dmr_chr, dmr_start, dmr_stop)) %>%
  select(-(dmr_chr:dmr_stop))

cg_dmr_mc <-  cg_dmr_mc %>%
  mutate(ID = paste(chr, start, stop)) %>%
  select(-(chr:stop))

te_cg_dmr <-  te_cg_dmr %>%
  mutate(ID = paste(dmr_chr, dmr_start, dmr_stop)) %>%
  select(-(dmr_chr:dmr_stop))

# function to get mC values from mC dataframe given DMR ID and list of accessions
lookupMeth <- function(df, acc, TEID) {
  x <- vector()
  for(i in acc) {
    x <- c(x, as.numeric(df[df$ID == TEID, i])) # match TE ID, get mC value for accession, append to list
  }
  if(length(x) == 0) { return(NA) } else { return(x) }
}

# function that takes list of pos/neg accession, returns r2 if length of each list > x
calc_r2 <- function(pos_acc, neg_acc, data, TEID, minimum) {
  # first filter each list of accessions to include only those that have DNA methylation data
  pos <- intersect(unlist(strsplit(pos_acc, ",")), colnames(data))
  neg <- intersect(unlist(strsplit(neg_acc, ",")), colnames(data))
  l_pos <- length(pos)
  l_neg <- length(neg)
  # if we have enough in each group, find the correlation
  if(l_pos > minimum && l_neg > minimum) {
    data_pos <- na.omit(lookupMeth(data, pos, TEID))
    data_neg <- na.omit(lookupMeth(data, neg, TEID))
    lp <- length(data_pos)
    ln <- length(data_neg)
    if(lp > minimum && ln > minimum) {
      # run correlation
      p <- rep(1, lp)
      n <- rep(0, ln)
      r2 <- cor(c(p, n), c(data_pos, data_neg))
      mC_pos <- mean(data_pos)
      mC_neg <- mean(data_neg)
      return(paste(r2, mC_pos, mC_neg, sep = "_"))
    } else {
      return(paste(NA, NA, NA, sep = "_"))
    }
  } else {
    return(paste(NA, NA, NA, sep = "_"))
  }
}

# make column names match
colnames(c_dmr_mc) <- gsub("_", "-", colnames(c_dmr_mc))
colnames(cg_dmr_mc) <- gsub("_", "-", colnames(cg_dmr_mc))

te_c_dmr <- te_c_dmr %>%
  rowwise() %>%
  mutate(d = calc_r2(pos_accessions, neg_accessions, c_dmr_mc, ID, 3)) %>%
  separate(d, c("r2", "mC_pos", "mC_neg"), sep = "_", remove = TRUE)

te_cg_dmr <- te_cg_dmr %>%
  rowwise() %>%
  mutate(d = calc_r2(pos_accessions, neg_accessions, cg_dmr_mc, ID, 3)) %>%
  separate(d, c("r2", "mC_pos", "mC_neg"), sep = "_", remove = TRUE)

# save data
write_tsv(te_c_dmr, "../ProcessedData/c_dmr_correlations.tsv")
write_tsv(te_cg_dmr, "../ProcessedData/cg_dmr_correlations.tsv")

### C-DMRs ###
# Make dataframe
c_dmr_info <- te_c_dmr %>%
  mutate(class = ifelse(distance < 1000, "TE-DMR", "Non-TE-DMR")) %>% 
  select(AbsenceClassification, class, mC_pos, mC_neg, r2, distance) %>% 
  gather(TE_present, mC, mC_pos:mC_neg)

c_dmr_info$TE_present[c_dmr_info$TE_present == "mC_pos"] <- TRUE
c_dmr_info$TE_present[c_dmr_info$TE_present == "mC_neg"] <- FALSE

c_dmr_info$AbsenceClassification[c_dmr_info$AbsenceClassification == "No insertion"] <- "Insertion"
c_dmr_info$AbsenceClassification[c_dmr_info$AbsenceClassification == "True deletion"] <- "Deletion"

# sapply(c_dmr_info, mode)

# convert to numeric
c_dmr_info <- transform(c_dmr_info, r2 = as.numeric(r2), mC = as.numeric(mC))

filtered_cdmr <- filter(c_dmr_info, AbsenceClassification == "Deletion" | AbsenceClassification == "Insertion")

# Plots
# r2 vs distance to DMR
ggplot(filtered_cdmr, aes(r2, distance/1000)) +
  geom_point(alpha=0.1) + theme_bw() + facet_wrap(~AbsenceClassification) +
  ylab("Distance to C-DMR (kb)") + xlab("r2") +
  ggsave("../Plots/CDMR/c_dmr_distance_vs_r2.pdf", height=4, width = 6, useDingbats=F)

# on log scale with linear regression
ggplot(filtered_cdmr, aes(r2, distance)) +
  geom_point(alpha=0.1, size=0.4) + theme_bw() + facet_wrap(~AbsenceClassification) +
  scale_y_log10() + geom_smooth(method = "lm") +
  xlim(-1,1) +
  ylab("Distance to C-DMR") + xlab("r2") +
  ggsave("../Plots/CDMR/c_dmr_distance_vs_r2_log.png", height=6, width = 10, units = "cm", dpi = 600)

# Distribution of r2 values for TE-DMRs vs non-TE-DMRs, for insertions and deletions
ggplot(filtered_cdmr, aes(r2, col=class)) +
  ylab(expression(italic("F"['n']*"(x)"))) + ggtitle("Pearson correlation values") +
  stat_ecdf() + theme_bw() + facet_wrap(~AbsenceClassification) +
  scale_color_brewer(palette = "Set1") +
  ggsave("../Plots/CDMR/r2_distribution_insertions_deletions_te_dmrs.pdf", height=3, width = 6, useDingbats=F)

# mC density plots for deletions vs insertions, TE-DMRs vs non-TE-DMRs
ggplot(filtered_cdmr, aes(fill=TE_present, mC)) + theme_bw() +
  geom_density(alpha = 0.5) + facet_wrap(AbsenceClassification~class) +
  scale_fill_brewer(palette = "Set1") +
  ggsave("../Plots/CDMR/mC_distribution_te_cdmr_insertion_deletion.pdf", height=4, width = 6, useDingbats=F)

# ecdf of mC values for TE-DMRs vs non-TE-DMRs, insertions vs deletions
ggplot(filtered_cdmr, aes(color=TE_present, mC)) + theme_bw() +
  stat_ecdf() + facet_wrap(AbsenceClassification~class) +
  scale_color_brewer(palette = "Set1") +
  ylab(expression(italic("F"['n']*"(x)"))) + ggtitle("DNA methylation") +
  ggsave("../Plots/CDMR/ecdf_mC_te_c_dmrs_insertions_deletions.pdf", width = 6, height = 6, useDingbats=F)

# boxplots of mC values for TE-DMRs vs non-TE-DMRs, insertions vs deletions
ggplot(filtered_cdmr, aes(TE_present, mC, fill=AbsenceClassification)) + theme_bw() +
  geom_boxplot() + facet_wrap(AbsenceClassification~class) +
  xlab("TE present") + ylab("mC / C") +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(legend.position = "none") +
  ggsave("../Plots/CDMR/c_dmr_mc_boxplots.pdf", height=6, width = 3, useDingbats=F)

### CG-DMRs ###
# Make dataframe
cg_dmr_info <- te_cg_dmr %>%
  mutate(class = ifelse(distance < 1000, "TE-DMR", "Non-TE-DMR")) %>%
  select(AbsenceClassification, class, mC_pos, mC_neg, r2, distance) %>%
  gather(TE_present, mC, mC_pos:mC_neg)

cg_dmr_info$TE_present[cg_dmr_info$TE_present == "mC_pos"] <- TRUE
cg_dmr_info$TE_present[cg_dmr_info$TE_present == "mC_neg"] <- FALSE

cg_dmr_info$AbsenceClassification[cg_dmr_info$AbsenceClassification == "No insertion"] <- "Insertion"
cg_dmr_info$AbsenceClassification[cg_dmr_info$AbsenceClassification == "True deletion"] <- "Deletion"

# sapply(cg_dmr_info, mode)

# convert to numeric
cg_dmr_info <- transform(cg_dmr_info, r2 = as.numeric(r2), mC = as.numeric(mC))

filtered_cgdmr <- filter(cg_dmr_info, AbsenceClassification == "Deletion" | AbsenceClassification == "Insertion")

# Plots
# r2 vs distance to DMR
ggplot(filtered_cgdmr, aes(r2, distance/1000)) +
  geom_point(alpha=0.1) + theme_bw() + facet_wrap(~AbsenceClassification) +
  ylab("Distance to CG-DMR (kb)") + xlab("r2") +
  ggsave("../Plots/CGDMR/cg_dmr_distance_vs_r2.pdf", height=3, width = 5, useDingbats=F)

# on log scale with linear regression
ggplot(filtered_cgdmr, aes(r2, distance)) +
  geom_point(alpha=0.1, size=0.4) + theme_bw() + facet_wrap(~AbsenceClassification) +
  scale_y_log10() + geom_smooth(method = "lm") +
  xlim(-1,1) +
  ylab("Distance to CG-DMR") + xlab("r2") +
  ggsave("../Plots/CGDMR/cg_dmr_distance_vs_r2_log.png", height=6, width = 10, units = "cm", dpi = 600)

# Distribution of r2 values for TE-DMRs vs non-TE-DMRs, for insertions and deletions
ggplot(filtered_cgdmr, aes(r2, col=class)) +
  ylab(expression(italic("F"['n']*"(x)"))) + ggtitle("Pearson correlation values") +
  stat_ecdf() + theme_bw() + facet_wrap(~AbsenceClassification) +
  scale_color_brewer(palette = "Set1") +
  ggsave("../Plots/CGDMR/r2_distribution_insertions_deletions_te_cgdmrs.pdf", height=3, width = 6, useDingbats=F)

# mC density plots for deletions vs insertions, TE-DMRs vs non-TE-DMRs
ggplot(filtered_cgdmr, aes(fill=TE_present, mC)) + theme_bw() +
  geom_density(alpha = 0.5) + facet_wrap(AbsenceClassification~class) +
  scale_fill_brewer(palette = "Set1") +
  ggsave("../Plots/CGDMR/mC_distribution_te_cgdmr_insertion_deletion.pdf", height=4, width = 6, useDingbats=F)

# ecdf of mC values for TE-DMRs vs non-TE-DMRs, insertions vs deletions
ggplot(filtered_cgdmr, aes(color=TE_present, mC)) + theme_bw() +
  stat_ecdf() + facet_wrap(AbsenceClassification~class) +
  ylab(expression(italic("F"['n']*"(x)"))) + ggtitle("DNA methylation") +
  scale_color_brewer(palette = "Set1") +
  ggsave("../Plots/CGDMR/ecdf_mC_te_cg_dmrs_insertions_deletions.pdf", width = 6, height = 6, useDingbats=F)

# boxplots of mC values for TE-DMRs vs non-TE-DMRs, insertions vs deletions
ggplot(filtered_cgdmr, aes(TE_present, mC, fill=AbsenceClassification)) + theme_bw() +
  geom_boxplot() + facet_wrap(AbsenceClassification~class) +
  xlab("TE present") + ylab("mC / C") +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(legend.position = "none") +
  ggsave("../Plots/CGDMR/cg_dmr_mc_boxplots.pdf", height=6, width = 3, useDingbats=F)
