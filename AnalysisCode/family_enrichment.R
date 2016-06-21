library(readr)
library(dplyr)


all_te <- read_tsv("../RawData/TAIR9_TE.bed.gz", col_names = c("chrom", "start", "stop", "strand", "TE", "Family", "Superfamily"))
tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = T)

# get family and superfamily information
tepav_fam <- left_join(tepav, all_te, by="TE") %>%
  select(Absence_classification, MAF, Frequency_classification, Family, Superfamily)

all_te %>%
  group_by(Family) %>%
  mutate(family_count = n()) %>%
  ungroup() %>% group_by(Superfamily) %>%
  mutate(superfamily_count = n()) %>% 
  ungroup() %>%
  rowwise() %>%
  mutate(family_perc = family_count / nrow(all_te), superfamily_perc = superfamily_count / nrow(all_te))

tepav_fam_perc <- tepav_fam %>%
  group_by(Family) %>%
  summarise(count = n(), perc = count / nrow(tepav_fam) * 100)

all_fam_count <- all_te %>%
  group_by(Family) %>%
  summarise(count = n(), perc = count / nrow(all_te) * 100)

fam <- full_join(all_fam_count, tepav_fam_perc, by="Family") %>%
  mutate(Enrichment = perc.y - perc.x) %>%
  select(Family, Enrichment) %>%
  mutate(cat = ifelse(Enrichment > 0, "Up", "Down"))

fam$sample[fam$sample == "perc.x"] <- "All TEs"
fam$sample[fam$sample == "perc.y"] <- "Variable TEs"

ggplot(fam, aes(Family, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10),
        legend.position="none") +
  ggsave(height=5,width=28, units = "cm", filename="../Plots/family_enrichments.pdf", useDingbats=FALSE)

# superfam
tepav_sfam_perc <- tepav_fam %>%
  group_by(Superfamily) %>%
  summarise(count = n(), perc = count / nrow(tepav_fam) * 100)

all_sfam_count <- all_te %>%
  group_by(Superfamily) %>%
  summarise(count = n(), perc = count / nrow(all_te) * 100)

sfam <- full_join(all_sfam_count, tepav_sfam_perc, by="Superfamily") %>%
  mutate(Enrichment = perc.y - perc.x) %>%
  select(Superfamily, Enrichment) %>%
  mutate(cat = ifelse(Enrichment > 0, "Up", "Down"))

sfam$sample[sfam$sample == "perc.x"] <- "All TEs"
sfam$sample[sfam$sample == "perc.y"] <- "Variable TEs"

ggplot(sfam, aes(Superfamily, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10),
        legend.position="none") +
  ggsave(height=5,width=8, units = "cm", filename="../Plots/superfamily_enrichments.pdf", useDingbats=FALSE)

### TODO: insertions and deletions separately ###



  