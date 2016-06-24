library(dplyr)
library(readr)
library(ggplot2)

set.seed(1)

# load data
mC <- read.table("../ProcessedData/c_dmr_allC.tsv.gz", header = T)
mCG <- read.table("../ProcessedData/cg_dmrs_allC.tsv.gz", header = T)

TE_DMR <- read_tsv("../ProcessedData/TE_C_DMR_distances.bed.gz", col_names = F) %>%
  filter(X16 < 1000) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15) %>%
  filter(FrequencyClassification == "Rare")

TE_CG_DMR <- read_tsv("../ProcessedData/TE_CG_DMR_distances.bed.gz", col_names = F) %>%
  filter(X16 < 1000) %>%
  rename(dmr_chr = X1, dmr_start = X2, dmr_stop = X3,
         pos_accessions = X8, neg_accessions = X9,
         AbsenceClassification = X13, FrequencyClassification = X15) %>%
  filter(FrequencyClassification == "Rare")

# convert DNA methylation values to ranks
ranks_mc <- mC[0,4:ncol(mC)]
for(i in 1:nrow(mC)) {
  ranks_mc[i,] <- rank(mC[i,4:ncol(mC)], ties.method = "first")
}
ranks_mc <- cbind(chr = mC$chr, start = mC$start, stop = mC$stop, ranks_mc)

ranks_mCG <- mCG[0,4:ncol(mCG)]
for(i in 1:nrow(mCG)) {
  ranks_mCG[i,] <- rank(mCG[i,4:ncol(mCG)], ties.method = "first")
}
ranks_mCG <- cbind(chr = mCG$chr, start = mCG$start, stop = mCG$stop, ranks_mCG)

# save mC ranks
write_tsv(ranks_mc, "../ProcessedData/c_dmr_mc_ranks.tsv")
system("gzip ../ProcessedData/c_dmr_mc_ranks.tsv")

write_tsv(ranks_mc, "../ProcessedData/cg_dmr_mc_ranks.tsv")
system("gzip ../ProcessedData/cg_dmr_mc_ranks.tsv")

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

