library(readr)
library(dplyr)
library(xtable)
library(ggplot2)
library(RColorBrewer)


all_te <- read_tsv("../RawData/TAIR9_TE.bed.gz", col_names = c("chrom", "start", "stop", "strand", "TE", "Family", "Superfamily"))
tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = T)

# get family and superfamily information
tepav_fam <- left_join(tepav, all_te, by="TE") %>%
  select(Absence_classification, MAF, Frequency_classification, Family, Superfamily)

all_te <- all_te %>%
  group_by(Family) %>%
  mutate(family_count = n()) %>%
  ungroup() %>% group_by(Superfamily) %>%
  mutate(superfamily_count = n()) %>% 
  ungroup() %>%
  rowwise() %>%
  mutate(family_perc = family_count / nrow(all_te), superfamily_perc = superfamily_count / nrow(all_te))

# family
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

# fam$sample[fam$sample == "perc.x"] <- "All TEs"
# fam$sample[fam$sample == "perc.y"] <- "Variable TEs"

ggplot(fam, aes(Family, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10),
        legend.position="none") +
  ggsave(height=5,width=28, units = "cm", filename="../Plots/Enrichments/family_enrichments.pdf", useDingbats=FALSE)

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

# sfam$sample[sfam$sample == "perc.x"] <- "All TEs"
# sfam$sample[sfam$sample == "perc.y"] <- "Variable TEs"

ggplot(sfam, aes(Superfamily, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10),
        legend.position="none") +
  ggsave(height=5,width=8, units = "cm", filename="../Plots/Enrichments/superfamily_enrichments.pdf", useDingbats=FALSE)

### Deletions

tepav_del <- filter(tepav_fam, Absence_classification == "True deletion")

# superfam
tepav_del_sfam <- tepav_del %>%
  group_by(Superfamily) %>%
  summarise(count = n(), perc = count / nrow(tepav_del) * 100)

sfam_del <- full_join(all_sfam_count, tepav_del_sfam, by="Superfamily") %>%
  mutate(Enrichment = perc.y - perc.x) %>%
  select(Superfamily, Enrichment) %>%
  mutate(cat = ifelse(Enrichment > 0, "Up", "Down"))

ggplot(sfam_del, aes(Superfamily, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10),
        legend.position="none") +
  ggsave(height=5,width=8, units = "cm", filename="../Plots/Enrichments/superfamily_enrichments_deletions.pdf", useDingbats=FALSE)

# fam
tepav_del_fam <- tepav_del %>%
  group_by(Family) %>%
  summarise(count = n(), perc = count / nrow(tepav_del) * 100)

fam_del <- full_join(all_fam_count, tepav_del_fam, by="Family") %>%
  mutate(perc.y = ifelse(is.na(perc.y), 0, perc.y)) %>%
  mutate(Enrichment = perc.y - perc.x) %>%
  select(Family, Enrichment) %>%
  mutate(cat = ifelse(Enrichment > 0, "Up", "Down"))

del_pt1 <- fam_del[1:160,]
del_pt2 <- fam_del[160:320,]

ggplot(del_pt1, aes(Family, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6),
        legend.position="none") +
  ggsave(height=5,width=20, units = "cm", filename="../Plots/Enrichments/family_enrichments_deletions_pt1.pdf", useDingbats=FALSE)

ggplot(del_pt2, aes(Family, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6),
        legend.position="none") +
  ggsave(height=5,width=20, units = "cm", filename="../Plots/Enrichments/family_enrichments_deletions_pt2.pdf", useDingbats=FALSE)

### Insertions

tepav_ins <- filter(tepav_fam, Absence_classification == "No insertion")

# superfam
tepav_ins_sfam <- tepav_ins %>%
  group_by(Superfamily) %>%
  summarise(count = n(), perc = count / nrow(tepav_ins) * 100)

sfam_ins <- full_join(all_sfam_count, tepav_ins_sfam, by="Superfamily") %>%
  mutate(Enrichment = perc.y - perc.x) %>%
  select(Superfamily, Enrichment) %>%
  mutate(cat = ifelse(Enrichment > 0, "Up", "Down"))

ggplot(sfam_ins, aes(Superfamily, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10),
        legend.position="none") +
  ggsave(height=5,width=8, units = "cm", filename="../Plots/Enrichments/superfamily_enrichments_insertions.pdf", useDingbats=FALSE)

# fam
tepav_ins_fam <- tepav_ins %>%
  group_by(Family) %>%
  summarise(count = n(), perc = count / nrow(tepav_ins) * 100)

fam_ins <- full_join(all_fam_count, tepav_ins_fam, by="Family") %>%
  mutate(perc.y = ifelse(is.na(perc.y), 0, perc.y)) %>%
  mutate(Enrichment = perc.y - perc.x) %>%
  select(Family, Enrichment) %>%
  mutate(cat = ifelse(Enrichment > 0, "Up", "Down"))

ins_pt1 <- fam_ins[1:160,]
ins_pt2 <- fam_ins[160:320,]

ggplot(ins_pt1, aes(Family, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6),
        legend.position="none") +
  ggtitle("Insertions") +
  ggsave(height=5,width=20, units = "cm", filename="../Plots/Enrichments/family_enrichments_insertions_pt1.pdf", useDingbats=FALSE)

ggplot(ins_pt2, aes(Family, Enrichment, fill = cat)) + geom_bar(stat="identity", color="black") + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 6, direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6),
        legend.position="none") +
  ggtitle("Insertions") +
  ggsave(height=5,width=20, units = "cm", filename="../Plots/Enrichments/family_enrichments_insertions_pt2.pdf", useDingbats=FALSE)

### Both together ###

# superfamily
all <- left_join(sfam_ins, sfam_del, by="Superfamily") %>%
  gather(indel, Enrichment, Enrichment.x, Enrichment.y)

all$indel[all$indel == "Enrichment.x"] <- "Insertions"
all$indel[all$indel == "Enrichment.y"] <- "Deletions"

insertion_col <- brewer.pal(3, "Set2")[1]
deletion_col <- brewer.pal(3, "Set2")[3]

ggplot(all, aes(Superfamily, Enrichment, fill = indel)) + geom_bar(stat="identity", color="black", position="dodge") + theme_bw() +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=10)) +
  ggsave(height=6,width=12, units = "cm", filename="../Plots/Enrichments/superfamily_enrichments_indel.pdf", useDingbats=FALSE)

# family
all_fam <- left_join(fam_ins, fam_del, by="Family") %>%
  gather(indel, Enrichment, Enrichment.x, Enrichment.y)

all_fam$indel[all_fam$indel == "Enrichment.x"] <- "Insertions"
all_fam$indel[all_fam$indel == "Enrichment.y"] <- "Deletions"

all_fam <- arrange(all_fam, Family)

pt1 <- all_fam[1:160,]
pt2 <- all_fam[160:320,]
pt3 <- all_fam[320:480,]
pt4 <- all_fam[480:642,]

ggplot(pt1, aes(Family, Enrichment, fill = indel)) + geom_bar(stat="identity", color="black", position="dodge") + theme_bw() +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6)) +
  ggsave(height=6,width=25, units = "cm", filename="../Plots/Enrichments/family_enrichments_indel_pt1.pdf", useDingbats=FALSE)

ggplot(pt2, aes(Family, Enrichment, fill = indel)) + geom_bar(stat="identity", color="black", position="dodge") + theme_bw() +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6)) +
  ggsave(height=6,width=25, units = "cm", filename="../Plots/Enrichments/family_enrichments_indel_pt2.pdf", useDingbats=FALSE)

ggplot(pt3, aes(Family, Enrichment, fill = indel)) + geom_bar(stat="identity", color="black", position="dodge") + theme_bw() +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6)) +
  ggsave(height=6,width=25, units = "cm", filename="../Plots/Enrichments/family_enrichments_indel_pt3.pdf", useDingbats=FALSE)

ggplot(pt4, aes(Family, Enrichment, fill = indel)) + geom_bar(stat="identity", color="black", position="dodge") + theme_bw() +
  scale_fill_manual(values = c(deletion_col, insertion_col)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=6)) +
  ggsave(height=6,width=25, units = "cm", filename="../Plots/Enrichments/family_enrichments_indel_pt4.pdf", useDingbats=FALSE)

# save tables
all_fam %>%
  select(Family, indel, Enrichment) %>%
  spread(indel, Enrichment) %>%
  write_tsv("../ProcessedData/family_enrichments.tsv")

all %>%
  select(Superfamily, indel, Enrichment) %>%
  spread(indel, Enrichment) %>%
  write_tsv("../ProcessedData/superfamily_enrichments.tsv")
