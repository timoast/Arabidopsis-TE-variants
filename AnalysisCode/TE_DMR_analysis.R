library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Distance from each C-DMR to the closest TEPAV
# Contains duplicates where a DMR was equally close to 2 or more TE variants
te_c_dmr <- read_tsv("../ProcessedData/TE_C_DMR_distances.bed.gz", col_names = F) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3, dmr_feature = X4,
         pos_accessions = X9, neg_accessions = X10,
         AbsenceClassification = X14, MAF=X15, FrequencyClassification = X16, LD=X17,
         distance = X18) %>%
  mutate(TE_close = distance < 1000)

# Distance from each CG-DMR to the closest TEPAV
# Contains duplicates where a DMR was equally close to 2 or more TE variants
te_cg_dmr <- read_tsv("../ProcessedData/TE_CG_DMR_distances.bed.gz", col_names = F) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3, dmr_feature = X4,
         pos_accessions = X9, neg_accessions = X10,
         AbsenceClassification = X14, MAF=X15, FrequencyClassification = X16, LD=X17,
         distance = X18) %>%
  mutate(TE_close = distance < 1000)

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
     col = chh, xlab = "Distance (kb)", main = "Distance from C-DMR to nearest TE variant")
hist(te_cg_dmr$distance / 1000,
     breaks = 300, xlim = c(0, 20),
     col = cg, xlab = "Distance (kb)", main = "Distance from CG-DMR to nearest TE variant")
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

# Of the CG-DMRs with a TE variant within 1 kb, what percentage of those are TE deletion, what percentage are insertions
# te_cg_dmr %>%
#   group_by(AbsenceClassification, distance < 1000) %>%
#   summarise(count = n()) %>%
#   mutate(tot = sum(count), perc = count / tot * 100)

# percentage of DMRs with TE variant within 1 kb, compared to random regions replicated 10k times
# see ProcessingCode/bootstrap_te_variants.sh for source of random overlap percentages
c_dmrs <- read_tsv("../RawData/c_dmrs.tsv.gz", col_names = F)
cg_dmrs <- read_tsv("../RawData/cg_dmrs.tsv.gz", col_names = F)

n_cdmr <- nrow(c_dmrs)
n_cgdmr <- nrow(cg_dmrs)

bootstrap_c_dmr <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k.txt", col_names = F) %>%
  mutate(c_dmr_all = X1 / n_cdmr * 100) %>% select(-X1)

bootstrap_c_dmr_insertions_1kb <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k_insertions.txt", col_names = F) %>%
  mutate(c_dmr_ins = X1 / n_cdmr * 100) %>% select(-X1)

bootstrap_c_dmr_deletions_1kb <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k_deletions.txt", col_names = F) %>%
  mutate(c_dmr_dels = X1 / n_cdmr * 100) %>% select(-X1)

bootstrap_c_dmr_na_1kb <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k_intermediate.txt", col_names = F) %>%
  mutate(c_dmr_na = X1 / n_cdmr * 100) %>% select(-X1)

bootstrap_cg_dmr <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k.txt", col_names = F) %>%
  mutate(cg_dmr_all = X1 / n_cgdmr * 100) %>% select(-X1)

bootstrap_cg_dmr_insertions_1kb <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k_insertions.txt", col_names = F) %>%
  mutate(cg_dmr_ins = X1 / n_cgdmr * 100) %>% select(-X1)

bootstrap_cg_dmr_deletions_1kb <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k_deletions.txt", col_names = F) %>%
  mutate(cg_dmr_dels = X1 / n_cgdmr * 100) %>% select(-X1)

bootstrap_cg_dmr_na_1kb <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k_intermediate.txt", col_names = F) %>%
  mutate(cg_dmr_na = X1 / n_cgdmr * 100) %>% select(-X1)

bootstrapping <- cbind(bootstrap_c_dmr, bootstrap_c_dmr_insertions_1kb, bootstrap_c_dmr_deletions_1kb, bootstrap_c_dmr_na_1kb,
                       bootstrap_cg_dmr, bootstrap_cg_dmr_insertions_1kb, bootstrap_cg_dmr_deletions_1kb, bootstrap_cg_dmr_na_1kb)

# What is the percentage of all C-DMRs with a TE variant within 1 kb
close_c_dmrs <- te_c_dmr %>% 
  group_by(AbsenceClassification) %>%
  filter(distance < 1000) %>%
  select(dmr_chr, dmr_start, dmr_stop, AbsenceClassification) %>% 
  unique() %>%
  summarise(count = n())

close_cg_dmrs <- te_cg_dmr %>% 
  group_by(AbsenceClassification) %>%
  filter(distance < 1000) %>%
  select(dmr_chr, dmr_start, dmr_stop, AbsenceClassification) %>%
  unique() %>%
  summarise(count = n())

# 95% CI = qnorm(0.975) * sd(sample) / sqrt(length(sample))
ci <- function(d) {
  m <- mean(d)
  s <- sd(d)
  n <- length(d)
  z <- qnorm(0.975)
  ci <- (z * s) / sqrt(n)
  return(ci)
}

conf <- lapply(bootstrapping, ci)
means <- lapply(bootstrapping, mean)

overlaps <- data_frame(TEPAV = c("TE deletions", "TE insertions", "NA calls", "Total"),
                       observeC = c(close_c_dmrs[1,2][[1]]/n_cdmr * 100,
                                close_c_dmrs[2,2][[1]]/n_cdmr * 100,
                                close_c_dmrs[3,2][[1]]/n_cdmr * 100,
                                sum(close_c_dmrs[,2]/n_cdmr * 100)),
                       expectC = c(means$c_dmr_dels, means$c_dmr_ins, means$c_dmr_na, means$c_dmr_all),
                       CIC = c(conf$c_dmr_dels, conf$c_dmr_ins, conf$c_dmr_na, conf$c_dmr_all),
                       observeCG = c(close_cg_dmrs[1,2][[1]]/n_cgdmr * 100,
                                 close_cg_dmrs[2,2][[1]]/n_cgdmr * 100,
                                 close_cg_dmrs[3,2][[1]]/n_cgdmr * 100,
                                 sum(close_cg_dmrs[,2]/n_cgdmr * 100)),
                       expectCG = c(means$cg_dmr_dels, means$cg_dmr_ins, means$cg_dmr_na, means$cg_dmr_all),
                       CICG = c(conf$cg_dmr_dels, conf$cg_dmr_ins, conf$cg_dmr_na, conf$cg_dmr_all)
                       )

is.num <- sapply(overlaps, is.numeric)
overlaps[is.num] <- lapply(overlaps[is.num], signif, 2)

write_tsv(overlaps, "../ProcessedData/dmr_overlaps.tsv")

pdf("../Plots/TE-DMRs.pdf", height = 5, width = 5)
par(mfrow = (c(2,2)))
# plot percentage of C-DMRs near randomly selected genomcic coordinates
hist(bootstrapping$c_dmr_dels,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "C-DMRs near TE deletions", border = chh)
# add line showing pecentage near deletions
intersections <- close_c_dmrs[1,2][[1]] / n_cdmr * 100
abline(v=intersections)

hist(bootstrapping$c_dmr_ins,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "C-DMRs near TE insertions", border = chh)
intersections <- close_c_dmrs[2,2][[1]] / n_cdmr * 100
abline(v=intersections)

hist(bootstrapping$cg_dmr_dels,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "CG-DMRs near TE deletions", border = cg)
intersections <- close_cg_dmrs[1,2][[1]] / n_cgdmr * 100
abline(v=intersections)

hist(bootstrapping$cg_dmr_ins,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "CG-DMRs near TE insertions", border = cg)
intersections <- close_cg_dmrs[2,2][[1]] / n_cgdmr * 100
abline(v=intersections)
dev.off()

# TE-DMR methylation level
c_dmr_mc <- read_tsv("../ProcessedData/c_dmr_allC.tsv.gz", col_names = T, na = c("None"))
cg_dmr_mc <- read_tsv("../ProcessedData/cg_dmrs_allC.tsv.gz", col_names = T, na = c("None"))

# Add common ID (DMR)
c_dmr_mc <-  c_dmr_mc %>%
  mutate(ID = paste(chr, start, stop)) %>%
  select(-(chr:stop))

te_c_dmr <- mutate(te_c_dmr, ID = paste(dmr_chr, dmr_start, dmr_stop))

cg_dmr_mc <-  cg_dmr_mc %>%
  mutate(ID = paste(chr, start, stop)) %>%
  select(-(chr:stop))

te_cg_dmr <- mutate(te_cg_dmr, ID = paste(dmr_chr, dmr_start, dmr_stop))

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
  separate(d, c("r2", "mC_pos", "mC_neg"), sep = "_", remove = TRUE) %>%
  select(-ID)

te_cg_dmr <- te_cg_dmr %>%
  rowwise() %>%
  mutate(d = calc_r2(pos_accessions, neg_accessions, cg_dmr_mc, ID, 3)) %>%
  separate(d, c("r2", "mC_pos", "mC_neg"), sep = "_", remove = TRUE) %>%
  select(-ID)

# save data
write_tsv(te_c_dmr, "../ProcessedData/c_dmr_correlations.tsv")
write_tsv(te_cg_dmr, "../ProcessedData/cg_dmr_correlations.tsv")
system("gzip ../ProcessedData/c_dmr_correlations.tsv")
system("gzip ../ProcessedData/cg_dmr_correlations.tsv")

## WHAT GENOMIC FEATURES ARE TE-DMRS AND NON-TE-DMRS FOUND IN?
te_c_dmr %>%
  select(dmr_feature, TE_close, AbsenceClassification) %>%
  group_by(dmr_feature) %>%
  mutate(total = n()) %>%
  group_by(dmr_feature, TE_close) %>%
  mutate(count = n()) %>%
  rowwise() %>%
  filter(count > 200) %>%
  mutate(perc = count / total * 100) %>%
  select(dmr_feature, TE_close, perc, AbsenceClassification) %>%
  unique() %>%
  ggplot(., aes(dmr_feature, perc, fill=TE_close)) +
  geom_bar(stat="identity", position="dodge", color = "black") + theme_bw() +
  facet_wrap(~AbsenceClassification)

# How many C-DMRs are distant from a TE variant, but close to a fixed TE?
# Add distance of each C-DMR to closest TE. Find percentage of C-DMRs not close to a TE variant
# that are within 1 kb of a TE (will be a lot)
c_dmr_tepav_distance <- read_tsv("../ProcessedData/c_dmr_tepav_distance.tsv.gz",
                                 col_names = c("dmr_chr", "dmr_start", "dmr_stop", "dmr_tepav_distance"))
c_dmr_te_distance <- read_tsv("../ProcessedData/c_dmrs_closest_te.tsv.gz",
                              col_names = c("dmr_chr", "dmr_start", "dmr_stop", "dmr_te", "dmr_te_distance"))

dmr_distances <- left_join(c_dmr_tepav_distance, c_dmr_te_distance)

dmr_distances %>%
  filter(dmr_tepav_distance > 1000) %>%
  ggplot(., aes(dmr_te_distance)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 1000)

dmr_distances %>%
  filter(dmr_tepav_distance > 1000, dmr_te_distance < 1000) %>%
  count()

### C-DMRs ###
# Make dataframe
c_dmr_info <- te_c_dmr %>%
  mutate(class = ifelse(distance < 1000, "TE-DMR", "Non-TE-DMR")) %>% 
  select(AbsenceClassification, class, mC_pos, mC_neg, r2, distance, dmr_feature) %>% 
  gather(TE_present, mC, mC_pos:mC_neg)

c_dmr_info <- transform(c_dmr_info, r2 = as.numeric(r2), mC = as.numeric(mC), AbsenceClassification = as.character(AbsenceClassification))
c_dmr_info$TE_present[c_dmr_info$TE_present == "mC_pos"] <- TRUE
c_dmr_info$TE_present[c_dmr_info$TE_present == "mC_neg"] <- FALSE
c_dmr_info$AbsenceClassification[c_dmr_info$AbsenceClassification == "No insertion"] <- "Insertion"
c_dmr_info$AbsenceClassification[c_dmr_info$AbsenceClassification == "True deletion"] <- "Deletion"
c_dmr_info <- tbl_df(c_dmr_info)

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
  xlim(-1,1) +
  ggsave("../Plots/CDMR/r2_distribution_insertions_deletions_te_dmrs.pdf", height=3, width = 6, useDingbats=F)

# perform Kolmogorov-Smirnov tests
get_ks <- function(df, s) {
  filtered <- filter(df, AbsenceClassification == s)
  tedmr <- na.omit(filter(filtered, class == "TE-DMR")$r2)
  nontedmr <- na.omit(filter(filtered, class == "Non-TE-DMR")$r2)
  ks.test(tedmr, nontedmr)
}
get_ks(filtered_cdmr, "Insertion")
# D = 0.23019, p-value < 2.2e-16
get_ks(filtered_cdmr, "Deletion")
# D = 0.098078, p-value = 0.0005356

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

cg_dmr_info <- transform(cg_dmr_info, r2 = as.numeric(r2), mC = as.numeric(mC), AbsenceClassification = as.character(AbsenceClassification))

cg_dmr_info$TE_present[cg_dmr_info$TE_present == "mC_pos"] <- TRUE
cg_dmr_info$TE_present[cg_dmr_info$TE_present == "mC_neg"] <- FALSE
cg_dmr_info$AbsenceClassification[cg_dmr_info$AbsenceClassification == "No insertion"] <- "Insertion"
cg_dmr_info$AbsenceClassification[cg_dmr_info$AbsenceClassification == "True deletion"] <- "Deletion"
cg_dmr_info <- tbl_df(cg_dmr_info)

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
  ylab("Distance to CG-DMR") + xlab("r") +
  ggsave("../Plots/CGDMR/cg_dmr_distance_vs_r2_log.png", height=6, width = 10, units = "cm", dpi = 600)

# Distribution of r2 values for TE-DMRs vs non-TE-DMRs, for insertions and deletions
ggplot(filtered_cgdmr, aes(r2, col=class)) +
  ylab(expression(italic("F"['n']*"(x)"))) + ggtitle("Pearson correlation values") +
  stat_ecdf() + theme_bw() + facet_wrap(~AbsenceClassification) +
  scale_color_brewer(palette = "Set1") +
  xlim(-1,1) +
  ggsave("../Plots/CGDMR/r2_distribution_insertions_deletions_te_cgdmrs.pdf", height=3, width = 6, useDingbats=F)

# perform Kolmogorov-Smirnov tests
get_ks(filtered_cgdmr, "Insertion")
# D = 0.094844, p-value = 4.261e-08
get_ks(filtered_cgdmr, "Deletion")
# D = 0.073487, p-value = 0.001093

del_cg <- filter(filtered_cgdmr, AbsenceClassification == "Deletion")
del_cg_tedmr <- na.omit(filter(del_cg, class == "TE-DMR")$r2)
del_cg_nontedmr <- na.omit(filter(ins_cg, class == "Non-TE-DMR")$r2)
ks.test(ins_cg_tedmr, ins_cg_nontedmr)

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
