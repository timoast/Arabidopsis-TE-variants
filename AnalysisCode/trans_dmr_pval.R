library(readr)
library(dplyr)

# get the names of chromosomes we will compare
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop("Must supply c/cgdmr, two chromosomes and two output files", call.=FALSE)
} else {
  dmr_file <- args[1]
  chrom_te <- args[2]
  chrom_dmr <- args[3]
  outfile_sig <- args[4]
  outfile_cor <- args[5]
  outfile_insertions <- args[6]
  outfile_deletions <- args[7]
}

# Load data
# Will do insertions and deletions separately then combine and plot with different colors

dmr <- read_tsv(dmr_file, col_names = T, na = c("None")) %>%
  filter(chr == chrom_dmr) %>%
  select(-stop)

insertions <- read_tsv("../ProcessedData/TE_insertions_matrix.tsv.gz", col_names=T) %>%
  rowwise() %>%
  filter(min(accession_poly_count, accession_nopoly_count) > 7) %>% 
  filter(chrom == chrom_te) %>%
  select(-(stop:accession_nopoly_count))

deletions <- read_tsv("../ProcessedData/TE_deletions_matrix.tsv.gz", col_names=T) %>%
  rowwise() %>%
  filter(min(accession_poly_count, accession_nopoly_count) > 7) %>%
  filter(chrom == chrom_te) %>%
  select(-(stop:accession_nopoly_count))

# filter out DMRs with lots of NAs (>15%)
dmr <- dmr[rowSums(is.na(dmr[,2:ncol(dmr)])) < ((ncol(dmr) -1) * 0.15), ]

# save the position information. The order will stay the same, so can later match p values to chr positions
dmr_positions <- dmr[,1:2]
insertion_positions <- insertions[,1:2]
deletion_positions <- deletions[,1:2]

# Remove columns without mC data
colnames(insertions) <- gsub("-", "_", colnames(insertions))
colnames(deletions) <- gsub("-", "_", colnames(deletions))
havename <- intersect(names(dmr),names(insertions))

insertions <- insertions[,colnames(insertions) %in% havename]
deletions <- deletions[,colnames(deletions) %in% havename]
dmr <- dmr[,colnames(dmr) %in% havename]

# order the accessions the same way in both TE and DMR sets
insertions <- insertions[,order(names(insertions))]
deletions <- deletions[,order(names(deletions))]
dmr <- dmr[,order(names(dmr))]

# Perform correlations
cor_insertions_dmr <- cor(t(insertions), t(dmr), use="pairwise.complete.obs")
cor_deletions_dmr <- cor(t(deletions), t(dmr), use="pairwise.complete.obs")

# Now permute labels for each pairwise comparision 500 times to get distribution of expected values
# function returns a list of matricies
permute_labels <- function(dmr, te, n) {
  perm_matricies <- list()
  for(x in 1:n) {
    permuted <- dmr[, sample(ncol(dmr))]
    perm <- cor(t(te), t(permuted), use="pairwise.complete.obs")
    perm_matricies[[x]] <- perm
  }
  return(perm_matricies)
}

# Calculate p-value for each correlation based on above bootstrap
lookup_perm <- function(row, column, perm) {
  d <- c()
  for(x in 1:length(perm)) {
    d <- c(d, perm[[x]][row, column])
  }
  return(d)
}

estimate_pval <- function(real, perm) {
  pvals <- real[,]
  n <- length(perm)
  for(column in 1:ncol(real)) {
    for(row in 1:nrow(real)) {
      d <- lookup_perm(row, column, perm)
      p <- sum(abs(d) > abs(real[row, column])) / n
      pvals[row, column] <- p
    }
  }
  return(pvals)
}

# Get row and column number of each significant correlation, look up positions of TE and DMR
# Gives dataframe with coordinates of correlated DMR-TE pairs
get_points <- function(te_positions, dmr_positions, dat, cutoff, side) {
  df <- data.frame(dmr_chr = character(),
                   dmr_pos = integer(),
                   te_chr = character(),
                   te_pos = integer(),
                   stringsAsFactors = FALSE)
  x <- 1
  for(column in 1:ncol(dat)) {
    for(row in 1:nrow(dat)) {
      # te = row, dmr = column
      if(!is.na(dat[row, column])) {
        if(side == "lower") {
          if(dat[row, column] < cutoff) {
            df[x,] <- c(dmr_positions[column,1],
                        dmr_positions[column,2],
                        te_positions[row,1],
                        te_positions[row,2])
            x <- x + 1
          }
        } else if(side == "higher") {
          if(dat[row, column] > cutoff) {
            df[x,] <- c(dmr_positions[column,1],
                        dmr_positions[column,2],
                        te_positions[row,1],
                        te_positions[row,2])
            x <- x + 1
          }
        }
      }
    }
  }
  return(df)
}

# Uses a lot of memory so need to do one at a time and then delete the permuted dataframes
# insertions
perm_insertions_dmr <- permute_labels(dmr, insertions, 500)
pval_insertion_dmr <- estimate_pval(cor_insertions_dmr, perm_insertions_dmr)
rm(perm_insertions_dmr)
dmr_insertions <- get_points(insertion_positions, dmr_positions, pval_insertion_dmr, 0.001, "lower")
high_cor_insertions <- get_points(insertion_positions, dmr_positions, cor_insertions_dmr, 0.5, "higher")

# deletions
perm_deletions_dmr <- permute_labels(dmr, deletions, 500)
pval_deletions_dmr <- estimate_pval(cor_deletions_dmr, perm_deletions_dmr)
rm(perm_deletions_dmr)
dmr_deletions <- get_points(deletion_positions, dmr_positions, pval_deletions_dmr, 0.001, "lower")
high_cor_deletions <- get_points(deletion_positions, dmr_positions, pval_deletions_dmr, 0.5, "higher")

# We now have coordinates of all TE-DMRs pairs that were highly correlated
# and those that were significantly different than expected based on 500 bootstraps (2 separate dataframes)

# Join insertions and deletions, adding column with insertion/deletion
dmr_insertions <- mutate(dmr_insertions, class = "Insertion")
high_cor_insertions <- mutate(high_cor_insertions, class = "Insertion")

dmr_deletions <- mutate(dmr_deletions, class = "Deletion")
high_cor_deletions <- mutate(high_cor_deletions, class = "Deletion")

dmr_sig <- rbind(dmr_insertions, dmr_deletions)
dmr_cor <- rbind(high_cor_insertions, high_cor_deletions)

# Save the data
write_tsv(dmr_sig, outfile_sig)
write_tsv(dmr_cor, outfile_cor)
write_tsv(as.data.frame(pval_insertion_dmr), outfile_insertions)
write_tsv(as.data.frame(pval_deletions_dmr), outfile_deletions)
