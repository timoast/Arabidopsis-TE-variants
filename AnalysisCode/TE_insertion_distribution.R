library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# TEPID TE calls
tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = T)

tepav %>%
  mutate(Absence_call = grepl("Col-0", Accessions_TE_present)) %>%
  group_by(Absence_call, Absence_classification) %>%
  summarise(count = n()) %>%
  write_tsv("../ProcessedData/deletion_stats.tsv")

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

rare_col <- brewer.pal(6, "Set3")[4]
common_col <- brewer.pal(6, "Set3")[6]

insertion_col <- brewer.pal(3, "Set2")[1]
deletion_col <- brewer.pal(3, "Set2")[3]
na_col <- brewer.pal(9, "Blues")[2]

ggplot(tepav, aes(start, color = Frequency_classification)) +
  geom_density(adjust = 1/8) + scale_color_manual(values = c(rare_col, common_col)) +
  facet_wrap(~chromosome, nrow = 1) + theme_bw() +
  ggsave("chomosomal_dist_rare_vs_common_line_chart.pdf", path = "../Plots", height = 3, width = 20)

ggplot(tepav, aes(start, color = Absence_classification)) +
  geom_density(adjust = 1/8) + scale_color_manual(values = c(na_col, deletion_col, insertion_col)) +
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
high_ld <- get_distribution(filter(chr1, LD == "high")$start, 1000)[1:608]
mid_ld <- get_distribution(filter(chr1, LD == "mid")$start, 1000)[1:608]
low_ld <- get_distribution(filter(chr1, LD == "low")$start, 1000)[1:608]

# # group together and scale each to it's maximum value
# all_hist_scaled <- cbind(all_tepav, rare_te, common_te, te_ins, te_del)
# 
color <- (colorRampPalette(brewer.pal(9,"GnBu"))(100))
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
ld_features_chr1 <- cbind(high_ld, mid_ld, low_ld)

pdf("../Plots/heatmap_te_insertions_chr1_all_features.pdf", height = 4, width = 5)
image(all_features_chr1, col = color, xlab = "Chromosome 1")
dev.off()

pdf("../Plots/heatmap_te_insertions_ld.pdf", height = 4, width = 5)
image(ld_features_chr1, col = color, xlab = "Chromosome 1")
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

# Now compare insertion distribution of Type 1 and Type 2 elements (reviewer request)
# First add annotation of TE superfamily to each variant
te <- rename(all_te, TE = X5, family = X6, superfamily = X7) %>%
  select(TE, family, superfamily)

k <- select(te, superfamily) %>% rename(s = superfamily) %>% unique()
t <- c(1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, NA, 2, 1, 1, 1, 1)
names(t) <- k$s

annotations <- left_join(tepav, te) %>% mutate(class = t[superfamily])

annotations %>%
  filter(chromosome == "chr1", !is.na(Absence_classification), !is.na(class)) %>%
  ggplot(., aes(start, fill=Absence_classification)) +
  geom_density(adjust=1/7) +
  facet_wrap(Absence_classification~class) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../Plots/type_1_type_2_distribution.pdf", height = 5, width = 5, useDingbats=F)

# 3. Find frequency that insertions / deletions occur in different genomic features
# common header for most files
intersectHead <- c("chromsome", "start", "end", "TE", "accessions_te_present",
                   "accessions_te_absent", "present_count", "absent_count", "MinorAllele",
                   "AbsenceClassification", "MAF", "FrequencyClassification", "FeatureChrom",
                   "FeatureStart", "FeatureStop", "FeatureName")

# common header for upstream and downstream
geneHeader <- c("chromsome", "start", "end", "TE", "accessions_te_present",
                   "accessions_te_absent", "present_count", "absent_count", "MinorAllele",
                   "AbsenceClassification", "MAF", "FrequencyClassification", "FeatureChrom",
                   "FeatureStart", "FeatureStop", "GeneName", "Strand", "TAIR10", "Feature", "ID")


# upstream
upstreamIntersections <- read_tsv("../ProcessedData/GeneFeatures/gene_upstream_regions_intersections.bed.gz",
                                  col_names = geneHeader)
all_upstream <- read_tsv("../ProcessedData/GeneFeatures/gene_upstream_regions.bed.gz", col_names = F)

fractionUpstream <- nrow(upstreamIntersections) / nrow(all_upstream) * 100

# downstream
downstreamIntersections <- read_tsv("../ProcessedData/GeneFeatures/gene_downstream_regions_intersections.bed.gz",
                                    col_names = geneHeader)
all_downstream <- read_tsv("../ProcessedData/GeneFeatures/gene_downstream_regions.bed.gz",
                           col_names = F)

fractionDownstream <- nrow(downstreamIntersections) / nrow(all_downstream) * 100

# exons
exon_intersections <- read_tsv("../ProcessedData/GeneFeatures/exon_intersections.bed.gz",
                               col_names = intersectHead)
all_exons <- read_tsv("../ProcessedData/GeneFeatures/exons.bed.gz",
                      col_names = c("chromosome", "start", "end", "geneName"))

fractionExon <- nrow(exon_intersections) / nrow(all_exons) * 100

# introns
intron_intersections <- read_tsv("../ProcessedData/GeneFeatures/intron_intersections.bed.gz",
                                 col_names = intersectHead)
all_intron <- read_tsv("../ProcessedData/GeneFeatures/introns.bed.gz", 
                       col_names = F)
fractionIntrons <- nrow(intron_intersections) / nrow(all_intron) * 100

# utr5
utr5Intersections <- read_tsv("../ProcessedData/GeneFeatures/utr5_intersections.bed.gz",
                              col_names = intersectHead)
all_utr5 <- read_tsv("../ProcessedData/GeneFeatures/utr5.bed.gz",
                     col_names = F)
fractionUTR5 <- nrow(utr5Intersections) / nrow(all_utr5) * 100

# utr3
utr3Intersections <- read_tsv("../ProcessedData/GeneFeatures/utr3_intersections.bed.gz",
                              col_names = intersectHead)
all_utr3 <- read_tsv("../ProcessedData/GeneFeatures/utr3.bed.gz",
                     col_names = F)
fractionUTR3 <- nrow(utr3Intersections) / nrow(all_utr3) * 100

# pseudogene
pseudogeneIntersections <- read_tsv("../ProcessedData/GeneFeatures/pseudogene_intersections.bed.gz",
                                    col_names = intersectHead)
all_pseudo <- read_tsv("../ProcessedData/GeneFeatures/pseudogene.bed.gz",
                       col_names = F)
fractionPseudogenes <- nrow(pseudogeneIntersections) / nrow(all_pseudo) * 100

# TE
teIntersections <- read_tsv("../ProcessedData/GeneFeatures/te_intersections.bed.gz", 
                            col_names = intersectHead)
# all TEs already loaded
fractionTE <- nrow(teIntersections) / nrow(all_te) * 100

# DHS
dhsIntersections <- read_tsv("../ProcessedData/GeneFeatures/dhs_intersections.bed.gz",
                             col_names = F)
all_dhs <- read_tsv("../RawData/Sullivan_DHS_PE_peaks_control.bed.gz", col_names = F)
fractionDHS <- nrow(dhsIntersections) / nrow(all_dhs) * 100

# Intergenic
intergenicIntersections <- read_tsv("../ProcessedData/GeneFeatures/intergenic_intersections.bed.gz",
                                    col_names = FALSE)

# now look at rare/common insertion/deletion frequency in each feature
# make a single dataframe with all the features, absence classification and frequency classification

df <- mutate(exon_intersections, feature = "Exon") %>% select(feature, FrequencyClassification, AbsenceClassification)

df <- intron_intersections %>%
  mutate(feature = "Intron") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- upstreamIntersections %>%
  mutate(feature = "Upstream") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- downstreamIntersections %>%
  mutate(feature = "Downstream") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- utr3Intersections %>%
  mutate(feature = "3' UTR") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- utr5Intersections %>%
  mutate(feature = "5' UTR") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- pseudogeneIntersections %>%
  mutate(feature = "Pseudogene") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- teIntersections %>%
  mutate(feature = "TE") %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- dhsIntersections %>%
  mutate(feature = "DHS", FrequencyClassification = X12, AbsenceClassification = X10) %>%
  select(feature, FrequencyClassification, AbsenceClassification) %>%
  rbind(df)

df <- intergenicIntersections %>%
        mutate(feature = "Intergenic", FrequencyClassification = X12, AbsenceClassification = X10) %>%
        select(feature, FrequencyClassification, AbsenceClassification) %>%
        rbind(df)

df %>%
  group_by(feature) %>%
  mutate(rare_perc = sum(FrequencyClassification == "Rare") / n() * 100,
         common_perc = sum(FrequencyClassification == "Common") / n() * 100,
         deletion_perc = sum(AbsenceClassification == "True deletion", na.rm = TRUE) / n() * 100,
         insertion_perc = sum(AbsenceClassification == "No insertion", na.rm = TRUE) / n() * 100,
         del_na_perc = sum(is.na(AbsenceClassification)) / n() * 100) %>%
  select(feature, rare_perc:del_na_perc) %>% unique(.) -> insertionStats

# make plots
feature_order <-  c("DHS", "Upstream", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Pseudogene", "TE", "Intergenic")

insertionStats %>%
  select(feature, rare_perc, common_perc) %>%
  melt() %>%
  ggplot(., aes(feature, value, fill = variable)) +
  geom_bar(stat = "identity", color = "Black") +theme_bw() +
  ylab("Percentage") + ggtitle("Rare vs common") +
  scale_fill_manual(values=c(rare_col, common_col)) +
  scale_x_discrete(limits = feature_order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), legend.position="none") +
  ggsave(filename = "rare_vs_common_genomic_feaures.pdf",
         path = "../Plots", height = 8, width = 5, units = "cm")

insertionStats %>%
  select(feature, insertion_perc, deletion_perc, del_na_perc) %>%
  melt() %>%
  ggplot(., aes(feature, value, fill = variable)) +
  geom_bar(stat = "identity", color = "Black") + theme_bw() + 
  scale_fill_manual(values=c(insertion_col, deletion_col, na_col)) +
  ylab("Percentage") + ggtitle("Insertions vs deletions") +
  scale_x_discrete(limits = feature_order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), legend.position="none") +
  ggsave(filename = "insertion_vs_deletion_genomic_feaures.pdf",
       path = "../Plots", height = 8, width = 5, units = "cm")

tepav$LD <- factor(tepav$LD, levels = c("high", "mid", "low", NA))

low <- brewer.pal(n = 3, "PuRd")[1]
mid <- brewer.pal(n = 3, "PuRd")[2]
high <- brewer.pal(n = 3, "PuRd")[3]

ggplot(tepav, aes(LD, MAF, fill=LD)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  scale_fill_brewer(type="seq", palette = 11, direction = -1) +
  coord_flip() +
  theme(legend.position="none") +
  ggsave(filename = "MAF_vs_LD.pdf",
         path = "../Plots", height = 4, width = 10, units = "cm")

data <- select(tepav, Absence_classification, LD)

testable <- filter(tepav, !is.na(LD))
testable %>% group_by(LD) %>% summarise(count = n(), perc = count/nrow(testable)*100)

data <- na.omit(data) %>%
  group_by(Absence_classification, LD) %>% summarise(count = n()) %>%
  mutate(n = sum(count)) %>%
  rowwise() %>% mutate(perc = count / n * 100)

ggplot(data, aes(Absence_classification, perc, fill=LD)) +
  geom_bar(stat="identity", color="black") +
  theme_bw() + scale_fill_manual(values = c(high, low, mid)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  ylab("Percentage") +
  coord_flip() +
  ggsave(filename = "LD_vs_indel.pdf",
         path = "../Plots", height = 4, width = 15, units = "cm")

ld <- read_tsv("../RawData/snp_association.tsv.gz")

ggplot(ld, aes(ranks_over_median, fill=flag)) + geom_histogram(color="black") +
  theme_bw() + scale_fill_manual(values = c(mid, high, low)) +
  geom_vline(xintercept = c(200, 400)) +
  theme(legend.position = "none") +
  ggsave(filename = "LD_histogram.pdf",
         path = "../Plots", height = 4, width = 10, units = "cm")
