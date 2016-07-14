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

get_centro_distance <- function(crds, centromeres) {
  cent <- filter(centromeres, chrom == crds[[1]])
  midpoint <- (cent[[3]] + cent[[2]]) / 2
  d <- abs(as.numeric(crds[[2]]) - midpoint)
  return(d)
}

centromeres <- read_tsv("../RawData/centromere_positions.txt", col_names = c("chrom", "start", "stop"))

data <- data %>% mutate(distance = get_centro_distance(c(chromosome, start, end), centromeres),
                        peri = distance < 3*10^6)

# Plot r vs distance to centromere
insertion_col <- brewer.pal(3, "Set2")[1]
deletion_col <- brewer.pal(3, "Set2")[3]

ggplot(data, aes(distance/1000000, r2, color=Absence_classification)) +
  geom_point(alpha=0.7, size=0.3) +
  facet_wrap(~Absence_classification, ncol=1) + theme_bw() +
  geom_hline(yintercept = 0) + ylim(c(-1, 1)) + 
  scale_color_manual(values = c("black", insertion_col, deletion_col)) +
  geom_vline(xintercept = 3) +
  theme(text = element_text(size=8), legend.position="none") +
  ggsave("../Plots/mc_flanking_correlation_indel_distance_point.pdf", height = 15, width = 6, units = "cm", useDingbats = F)