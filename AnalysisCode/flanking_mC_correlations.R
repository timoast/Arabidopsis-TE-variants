library(dplyr)
library(readr)
library(ggplot2)

# Load TE variant data, add a unique ID that will allow tables to be joined later
tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz") %>%
  mutate(ID = paste(chromosome, start, end))

# Load DNA methylation data
# DNA methylation levels in TE variant site flanking +/- 200 bp regions, all C contexts
# Add a column with the unique TE variant ID, same as in the TEPAV table
mC_data <- read_tsv("../ProcessedData/DNAmethylation/flanking_allc_200_bp.tsv.gz",
                         na = c("None")) %>%
  mutate(ID = paste(chr, start, stop)) %>% select(-(c(chr, start, stop)))


# Convert column names to the format the accession names are in the TEPAV file
colnames(mC_data) <- gsub("mC_calls_", "", colnames(mC_data))
colnames(mC_data) <- gsub("_", "-", colnames(mC_data))

# remove any duplicate rows
mC_data <- mC_data[!duplicated(mC_data), ]

# Function to look up DNA methylation values and perform correlations
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
      return(r2)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}

data <- tepav %>%
  rowwise() %>%
  mutate(r2 = calc_r2(Accessions_TE_present, Accessions_TE_absent, mC_data, ID, 3))

data$LD <- factor(data$LD, levels = c("high", "mid", "low", NA))

# Plots
d <- filter(data, Absence_classification == "True deletion" | Absence_classification == "No insertion")
ggplot(d, aes(Absence_classification, r2, fill=LD)) +
  geom_boxplot(outlier.shape = 1) + theme_bw() +
  scale_fill_brewer(type="seq", palette = 11, direction = -1) +
  geom_hline(yintercept = 0) + ylim(c(-1, 1)) +
  coord_flip() + ylab("Correlation") + xlab("TE variant") + 
  theme(text = element_text(size=8)) +
  ggsave("../Plots/mc_flanking_correlation.pdf", height = 5, width = 10, units = "cm", useDingbats = F)

ggplot(data, aes(MAF, r2, fill=LD)) +
  geom_boxplot(aes(group = cut_width(MAF, 0.04)), outlier.shape = 1, outlier.size = 0.8) +
  theme_bw() + geom_hline(yintercept = 0) + facet_wrap(~LD, ncol = 1) +
  scale_fill_brewer(type="seq", palette = 11, direction = -1) +
  theme(text = element_text(size=8)) +
  ylim(-1,1) +
  ggsave("../Plots/mc_flanking_correlation_by_ld.pdf", height = 20, width = 8, units = "cm", useDingbats=F)
