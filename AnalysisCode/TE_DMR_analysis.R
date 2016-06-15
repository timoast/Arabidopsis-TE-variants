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

cg <- "#B4B464"
chg <- "#6665AD"
chh <- "#B29492"

hist(te_c_dmr$distance / 1000, breaks = 1000, xlim = c(0, 20), col = chh)
hist(te_cg_dmr$distance / 1000, breaks = 1000, xlim = c(0, 20), col = cg)

# need to compare to random sample

# need to break down by insertion / deletion and rare / common
hist(filter(te_c_dmr, FrequencyClassification == "Rare")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = chh)
hist(filter(te_c_dmr, FrequencyClassification == "Common")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = chh)

hist(filter(te_cg_dmr, FrequencyClassification == "Rare")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = cg)
hist(filter(te_cg_dmr, FrequencyClassification == "Common")$distance / 1000, breaks = 1000, xlim = c(0, 20), col = cg)

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

# TODO
# make function that randomly selects a given number of regions of a given length from the genome
# automate the intersection of these regions with DMRs to give number within 1 kb
# plot distribution of overlaps and number of observed overlaps for the real data
# do for C-DMRs and CG-DMRs

chr <- data.frame(chromosome = c("chr1","chr2","chr3","chr4","chr5"),
                          length = c( 30427671, 19698289, 23459830, 18585056, 26975502))

getRandomCoords <- function(chromosomes, l) {
  random_chr <- sample(chromosomes$chromosome, 1)
  max_len <- chromosomes[chromosomes$chromosome == random_chr, 2]
  start <- sample(1:max_len-l, 1)
  stop <- start + l
  return(c(as.character(random_chr), as.character(start), as.character(stop)))
}

tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz")
av_size <- round(mean(tepav$end - tepav$start))

random_coordinates <- replicate(nrow(tepav), getRandomCoords(chr, av_size))
random_coordinates <- t(random_coordinates)

