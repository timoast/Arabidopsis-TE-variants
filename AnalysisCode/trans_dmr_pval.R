library(readr)
library(dplyr)

# Load data
# Will do insertions and deletions separately then combine and plot with different colors
# Will create separate plots for C-DMRs and CG-DMRs
cdmr <- read_tsv("../ProcessedData/c_dmr_allC.tsv.gz", col_names = T, na = c("None")) %>%
  select(-stop)

cgdmr <- read_tsv("../ProcessedData/cg_dmrs_allC.tsv.gz", col_names = T, na = c("None")) %>%
  select(-stop)

insertions <- read_tsv("../ProcessedData/TE_insertions_matrix.tsv.gz", col_names=T) %>%
  rowwise() %>%
  filter(min(accession_poly_count, accession_nopoly_count) > 7) %>%
  select(-(stop:accession_nopoly_count))

deletions <- read_tsv("../ProcessedData/TE_deletions_matrix.tsv.gz", col_names=T) %>%
  rowwise() %>%
  filter(min(accession_poly_count, accession_nopoly_count) > 7) %>%
  select(-(stop:accession_nopoly_count))

# filter out DMRs with lots of NAs (>15%)
cdmr <- cdmr[rowSums(is.na(cdmr[,2:ncol(cdmr)])) < ((ncol(cdmr) -1) * 0.15), ]
cgdmr <- cgdmr[rowSums(is.na(cgdmr[,2:ncol(cgdmr)])) < ((ncol(cgdmr) -1) * 0.15), ]

# save the position information. The order will stay the same, so can later match p values to chr positions
cdmr_positions <- cdmr[,1:2]
cgdmr_positions <- cgdmr[,1:2]
insertion_positions <- insertions[,1:2]
deletion_positions <- deletions[,1:2]

# Remove columns without mC data
colnames(insertions) <- gsub("-", "_", colnames(insertions))
colnames(deletions) <- gsub("-", "_", colnames(deletions))
havename <- intersect(names(cdmr),names(insertions))

insertions <- insertions[,colnames(insertions) %in% havename]
deletions <- deletions[,colnames(deletions) %in% havename]

cdmr <- cdmr[,colnames(cdmr) %in% havename]
cgdmr <- cgdmr[,colnames(cgdmr) %in% havename]

# order the accessions the same way in both TE and DMR sets
insertions <- insertions[,order(names(insertions))]
deletions <- deletions[,order(names(deletions))]

cdmr <- cdmr[,order(names(cdmr))]
cgdmr <- cdmr[,order(names(cgdmr))]

# Perform correlations
# Insertions
cor_insertions_cdmr <- cor(t(insertions), t(cdmr), use="pairwise.complete.obs")
cor_insertions_cgdmr <- cor(t(insertions), t(cgdmr), use="pairwise.complete.obs")
# Deletions
cor_deletions_cdmr <- cor(t(deletions), t(cdmr), use="pairwise.complete.obs")
cor_deletions_cgdmr <- cor(t(deletions), t(cgdmr), use="pairwise.complete.obs")

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

# insertions
perm_insertions_cdmr <- permute_labels(cdmr, insertions, 500)
perm_insertions_cgdmr <- permute_labels(cgdmr, insertions, 500)

# deletions
perm_deletions_cdmr <- permute_labels(cdmr, deletions, 500)
perm_deletions_cgdmr <- permute_labels(cgdmr, deletions, 500)

# Calculate p-value for each correlation based on above bootstrap
lookup_perm <- function(row, column, perm) {
  d <- c()
  for(x in 1:length(perm)) {
    d <- c(d, perm_matricies[[x]][row, column])
  }
  return(d)
}

estimate_pval <- function(real, perm) {
  pvals <- real[,]
  n <- length(perm)
  for(column in 1:ncol(real)) {
    for(row in 1:nrow(real)) {
      d <- lookup_perm(row, column, perm_matricies)
      p <- sum(abs(d) > abs(real[row, column])) / n
      pvals[row, column] <- p
    }
  }
  return(pvals)
}

# insertions
pval_insertion_cdmr <- estimate_pval(cor_insertions_cdmr, perm_insertions_cdmr)
pval_insertion_cgdmr <- estimate_pval(cor_insertions_cgdmr, perm_insertions_cgdmr)

# deletions
pval_deletions_cdmr <- estimate_pval(cor_deletions_cdmr, perm_deletions_cdmr)
pval_deletions_cgdmr <- estimate_pval(cor_deletions_cgdmr, perm_deletions_cgdmr)

# Get row and column number of each significant correlation, look up positions of TE and DMR
# Gives dataframe with coordinates of correlated DMR-TE pairs
get_points <- function(te_positions, dmr_positions, dat, cutoff) {
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
        if(dat[row, column] < cutoff) {
          df[x,] <- c(dmr_positions[column,1],
                      dmr_positions[column,2],
                      te_positions[row,1],
                      te_positions[row,2])
          x <- x + 1
        }
      }
    }
  }
  return(df)
}

# Insertions
cdmr_insertions <- get_points(insertion_positions, cdmr_positions, pval_insertion_cdmr, 0.01)
cgdmr_insertions <- get_points(insertion_positions, cgdmr_positions, pval_insertion_cgdmr, 0.01)

# Deletions
cdmr_deletions <- get_points(deletion_positions, cdmr_deletions, pval_deletions_cdmr, 0.01)
cgdmr_deletions <- get_points(deletion_positions, cgdmr_deletions, pval_deletions_cgdmr, 0.01)

# Join insertions and deletions, adding column with insertion/deletion
cdmr_insertions <- mutate(cdmr_insertions, class = "Insertion")
cgdmr_insertions <- mutate(cgdmr_insertions, class = "Insertion")
cdmr_deletions <- mutate(cdmr_deletions, class = "Deletion")
cgdmr_deletions <- mutate(cgdmr_deletions, class = "Deletion")

cdmr_sig <- rbind(cdmr_insertions, cdmr_deletions)
cgdmr_sig <- rbind(cgdmr_insertions, cgdmr_deletions)

# Now plot as scatterplot



