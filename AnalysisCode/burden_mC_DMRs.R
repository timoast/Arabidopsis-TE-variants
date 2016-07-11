library(dplyr)
library(readr)
library(ggplot2)

set.seed(1)

# load data
mC <- read.table("../ProcessedData/c_dmr_allC.tsv.gz", header = T, na.strings = c("None"))
mCG <- read.table("../ProcessedData/cg_dmrs_allC.tsv.gz", header = T, na.strings = c("None"))

# exclude bud tissue methylomes (drop column if contains the word "bud")
mC <- mC[,!grepl("bud", colnames(mC))]
mCG <- mCG[,!grepl("bud", colnames(mCG))]

# Load the DMR data and filter where the DMR is within 1 kb of a TE variant, and the variant is rare
TE_DMR <- read_tsv("../ProcessedData/TE_C_DMR_distances.bed.gz", col_names = F) %>% 
  filter(X17 < 1000) %>% 
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15) %>%
  filter(FrequencyClassification == "Rare")

TE_CG_DMR <- read_tsv("../ProcessedData/TE_CG_DMR_distances.bed.gz", col_names = F) %>%
  filter(X17 < 1000) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15) %>%
  filter(FrequencyClassification == "Rare")

# convert DNA methylation values to ranks and exclude those with no coverage (NA)
# to maintain the same range of ranks for each TE variant, convert to percentiles excluding the NAs
# ie, we divide each rank by the number of non-NA ranks & multiply by 100,
# so the highest becomes 100 and lowest 1.
# We can them compare between ranks, and 100 will always be the highest in the population etc.

# CDMR
numNAsCDMR <- rowSums(is.na(mC[,4:ncol(mC)]))
ranks_mc <- mC[0,4:ncol(mC)]
for(i in 1:nrow(mC)) {
  l <- ncol(mC) - 3 - numNAsCDMR[i]  # -3 because there are 3 coordinate columns that are not ranks
  d <- rank(mC[i,4:ncol(mC)], ties.method = "first", na.last = "keep")  # keep NAs when calculating ranks
  ranks_mc[i,] <- d / l * 100
}
ranks_mc <- cbind(chr = mC$chr, start = mC$start, stop = mC$stop, ranks_mc)

#CGDMR
numNAsCGDMR <- rowSums(is.na(mCG[,4:ncol(mCG)]))
ranks_mCG <- mCG[0,4:ncol(mCG)]
for(i in 1:nrow(mCG)) {
  l <- ncol(mCG) - 3 - numNAsCGDMR[i]
  d <- rank(mCG[i,4:ncol(mCG)], ties.method = "first", na.last = "keep")
  ranks_mCG[i,] <- d / l * 100
}
ranks_mCG <- cbind(chr = mCG$chr, start = mCG$start, stop = mCG$stop, ranks_mCG)

# save mC ranks
write_tsv(ranks_mc, "../ProcessedData/c_dmr_mc_ranks.tsv")
system("gzip ../ProcessedData/c_dmr_mc_ranks.tsv")

write_tsv(ranks_mCG, "../ProcessedData/cg_dmr_mc_ranks.tsv")
system("gzip ../ProcessedData/cg_dmr_mc_ranks.tsv")

# Function to look up mC rank for the accession containing each rare TE variant
gatherData <- function(tepav, ranks) {
  data <- vector()
  for(i in 1:nrow(tepav)){
    row <- tepav[i,]
    indel <- ifelse(row[["AbsenceClassification"]] == "No insertion", "pos_accessions", "neg_accessions")
    accessions <- unlist(strsplit(gsub("-", "_", as.character(row[[indel]])), ","))
    accessions <- intersect(accessions, colnames(ranks))
    if(length(accessions) > 0){
      for(acc in accessions){
        dmr <- filter(ranks, chr == row[["dmr_chr"]], start == row[["dmr_start"]], stop == row[["dmr_stop"]])
        if(tryCatch(nrow(dmr) > 1)) {stop("non-unique DMR coordinates")}
        data <- c(data, as.numeric(dmr[1,acc]))
      }
    }
  }
  return(data)
}

# same function, but choose an accession rank at random rather than the one containing the rare variant
gatherDataRand <- function(tepav, ranks) {
  data <- vector()
  all_names <- colnames(ranks)[4:ncol(ranks)]
  for(i in 1:nrow(tepav)){
    row <- tepav[i,]
    indel <- ifelse(row[["AbsenceClassification"]] == "No insertion", "pos_accessions", "neg_accessions")
    accessions <- unlist(strsplit(gsub("-", "_", as.character(row[[indel]])), ","))
    accessions <- intersect(accessions, colnames(ranks))
    if(length(accessions) > 0){
      rand <- sample(all_names, length(accessions), replace = FALSE)
      for(acc in rand){
        dmr <- filter(ranks, chr == row[["dmr_chr"]], start == row[["dmr_start"]], stop == row[["dmr_stop"]])
        if(tryCatch(nrow(dmr) > 1)) {stop("non-unique DMR coordinates")}
        data <- c(data, as.numeric(dmr[1,acc]))
      }
    }
  }
  return(data)
}

# Get the TE deletion variants and TE insertion varaints separately
# Count the mC ranks for each TE variant
# Repeat for a random selection of accession names
# C-DMRs
deletions <- filter(TE_DMR, AbsenceClassification == "True deletion")
insertions <- filter(TE_DMR, AbsenceClassification == "No insertion")

deletion.ranks <- gatherData(deletions, ranks_mc)
insertion.ranks <- gatherData(insertions, ranks_mc)

# randomized
deletion.rand <- gatherDataRand(deletions, ranks_mc)
insertion.rand <- gatherDataRand(insertions, ranks_mc)

# CG-DMRs
cg.deletions <- filter(TE_CG_DMR, AbsenceClassification == "True deletion")
cg.insertions <- filter(TE_CG_DMR, AbsenceClassification == "No insertion")

cg.deletion.ranks <- gatherData(cg.deletions, ranks_mCG)
cg.insertion.ranks <- gatherData(cg.insertions, ranks_mCG)

# randomized
cg.deletion.rand <- gatherDataRand(cg.deletions, ranks_mCG)
cg.insertion.rand <- gatherDataRand(cg.insertions, ranks_mCG)

# plots
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

dot_plot <- function(real, random, title, brk = 72) {
  
  real.hist <- hist(real, breaks = brk, plot = F)
  rand.hist <- hist(random, breaks = brk, plot = F)
  
  y_range <- c(min(real.hist$counts) - 10, max(real.hist$counts) + 10)
  if(y_range[1] < 0) { y_range[1] = 0 }
  
  model.real <- lm(real.hist$counts ~ poly(real.hist$breaks[2:length(real.hist$breaks)], 2, raw=TRUE))
  rsq.real <- signif(summary(model.real)$adj.r.squared, 3)
  p.real <- signif(lmp(model.real), 3)
  
  model.rand <- lm(rand.hist$counts ~ poly(rand.hist$breaks[2:length(rand.hist$breaks)], 2, raw=TRUE))
  rsq.rand <- signif(summary(model.rand)$adj.r.squared, 3)
  p.rand <- signif(lmp(model.rand), 3)
  
  plot(real.hist$counts, pch=19, col = "blue", las = 1, cex = 0.6,
       main = paste(title, " real"), ylim = y_range, xlab = "DNA methylation level rank", ylab = "Rare variants")
  lines(model.real$fitted.values)
  legend("top", paste("Rsq = ", rsq.real, ", p = ", p.real), bty = "n")
  
  plot(rand.hist$counts, pch=19, col = "blue", las = 1, cex = 0.6,
       main = paste(title, " random"), ylim = y_range, xlab = "DNA methylation level rank", ylab = "Rare variants")
  lines(model.rand$fitted.values)
  legend("top", paste("Rsq = ", rsq.rand, ", p = ", p.rand), bty = "n")
}

pdf("../Plots/burden_c_dmr_methylation_dot_plots.pdf", width = 6, height = 6, useDingbats=FALSE)
par(mfrow=c(2,2))
dot_plot(deletion.ranks, deletion.rand, "Deletions")
dot_plot(insertion.ranks, insertion.rand, "Insertions")
dev.off()

pdf("../Plots/burden_cg_dmr_methylation_dot_plots.pdf", width = 6, height = 6, useDingbats=FALSE)
par(mfrow=c(2,2))
dot_plot(cg.deletion.ranks, cg.deletion.rand, "Deletions")
dot_plot(cg.insertion.ranks, cg.insertion.rand, "Insertions")
dev.off()

