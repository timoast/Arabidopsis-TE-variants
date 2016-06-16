library(readr)
library(dplyr)
library(ggplot2)

te_c_dmr <- read_tsv("../ProcessedData/TE_C_DMR_distances.bed.gz", col_names = F) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15,
         distance = X16)

te_cg_dmr <- read_tsv("../ProcessedData/TE_CG_DMR_distances.bed.gz", col_names = F) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15,
         distance = X16)

random_c_dmr <- read_tsv("../ProcessedData/random_coords_closest_c_dmr.bed.gz", col_names = F)
random_cg_dmr <- read_tsv("../ProcessedData/random_coords_closest_cg_dmr.bed.gz", col_names = F)

cg <- "#B4B464"
chg <- "#6665AD"
chh <- "#B29492"

pdf("../Plots/distances_TEAPV_DMR.pdf", height = 3, width = 3)
hist(te_c_dmr$distance / 1000,
     breaks = 200, xlim = c(0, 20),
     col = chh, xlab = "Distance (kb)", main = "Distance from TE variant to nearest C-DMR")
hist(te_cg_dmr$distance / 1000,
     breaks = 300, xlim = c(0, 20),
     col = cg, xlab = "Distance (kb)", main = "Distance from TE variant to nearest CG-DMR")
dev.off()

# # need to break down by insertion / deletion and rare / common
# hist(filter(te_c_dmr, FrequencyClassification == "Rare")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = chh)
# hist(filter(te_c_dmr, FrequencyClassification == "Common")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = chh)
# 
# hist(filter(te_cg_dmr, FrequencyClassification == "Rare")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = cg)
# hist(filter(te_cg_dmr, FrequencyClassification == "Common")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = cg)

# C-DMRs
ggplot(te_c_dmr, aes(distance / 1000, fill = FrequencyClassification)) +
  geom_density(alpha = 0.5) + facet_wrap(~FrequencyClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  ggsave("../Plots/c_dmr_distances_common_vs_rare.pdf", width = 10, height = 6, units = "cm", useDingbats=F)	

ggplot(te_c_dmr, aes(distance / 1000, fill = AbsenceClassification)) +
  geom_density(alpha = 0.5) + facet_wrap(~AbsenceClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  ggsave("../Plots/c_dmr_distances_insertion_vs_deletion.pdf", width = 15, height = 6, units = "cm", useDingbats=F)

# CG-DMRs
ggplot(te_cg_dmr, aes(distance / 1000, fill = FrequencyClassification)) +
  geom_density(alpha = 0.5) + facet_wrap(~FrequencyClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  ggsave("../Plots/cg_dmr_distances_common_vs_rare.pdf", width = 10, height = 6, units = "cm", useDingbats=F)	

ggplot(te_cg_dmr, aes(distance / 1000, fill = AbsenceClassification)) +
  geom_density(alpha = 0.5) + facet_wrap(~AbsenceClassification) + xlim(0, 5) +
  theme_bw() + xlab("Distance (kb)") + ylab("Density") +
  theme(legend.position="none") +	
  ggsave("../Plots/cg_dmr_distances_insertion_vs_deletion.pdf", width = 15, height = 6, units = "cm", useDingbats=F)

# percentage of DMRs with TE variant within 1 kb, compared to random regions replicated 10k times
# see ProcessingCode/bootstrap_te_variants.sh for source of random overlap percentages
c_dmrs <- read_tsv("../RawData/c_dmrs.tsv.gz", col_names = F)
cg_dmrs <- read_tsv("../RawData/cg_dmrs.tsv.gz", col_names = F)

bootstrap_c_dmr_1kb <- read_table("../ProcessedData/random_selections_c_dmr_rep_10k.txt", col_names = F)
bootstrap_cg_dmr_1kb <- read_table("../ProcessedData/random_selections_cg_dmr_rep_10k.txt", col_names = F)

pdf("../Plots/DMR_intersections_bootstrap.pdf", height = 3, width = 3)
hist(bootstrap_c_dmr_1kb$X1/nrow(c_dmrs) * 100,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "TE variants near C-DMRs", border = chh)
abline(v=nrow(te_c_dmr %>% filter(distance < 1000)) / nrow(c_dmrs) * 100)

hist(bootstrap_cg_dmr_1kb$X1/nrow(cg_dmrs) * 100,
     breaks = 10, xlim = c(0, 100), xlab = "Percent within 1 kb",
     main = "TE variants near CG-DMRs", border = cg)
abline(v=nrow(te_cg_dmr %>% filter(distance < 1000)) / nrow(cg_dmrs) * 100)
dev.off()

# TE-DMR methylation level



