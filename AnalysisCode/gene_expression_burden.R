library(dplyr)

set.seed(1)

# this data only has accession names with RNA data available
expression <- read.table("../RawData/gene_expression.tsv.gz", sep = ""), header = T, row.names = 1)
upstream <- read.table('../ProcessedData/upstream_insertion.tsv', sep = ""), header = F)
exon <- read.table('../ProcessedData/exon_insertion.tsv', sep = ""), header = F)
intron <- read.table("../ProcessedData/intron_insertion.tsv", sep = ""), header = F)
utr5 <- read.table("../ProcessedData/utr5_insertion.tsv", sep = ""), header = F)
utr3 <- read.table("../ProcessedData/utr3_insertion.tsv", sep = ""), header = F)
downstream <- read.table("../ProcessedData/downstream_insertion.tsv", sep = ""), header = F)

upstream <- upstream[!duplicated(upstream),]
exon <- exon[!duplicated(exon),]
intron <- intron[!duplicated(intron),]
utr5 <- utr5[!duplicated(utr5),]
utr3 <- utr3[!duplicated(utr3),]
downstream <- downstream[!duplicated(downstream),]

header <- c("gene", 'TE', 'pos_accessions', 'neg_accessions')
colnames(upstream) <- header
colnames(exon) <- header
colnames(intron) <- header
colnames(utr3) <- header
colnames(utr5) <- header
colnames(downstream) <- header

upstream[] <- lapply(upstream, as.character)
exon[] <- lapply(exon, as.character)
utr3[] <- lapply(utr3, as.character)
utr5[] <- lapply(utr5, as.character)
intron[] <- lapply(intron, as.character)
downstream[] <- lapply(downstream, as.character)

find.len <- function(d) {
  return(length(unlist(strsplit(d, ","))))
}

upstream %>%
  rowwise() %>%
  filter(find.len(pos_accessions) <= 4) -> upstream.five

exon %>%
  rowwise() %>%
  filter(find.len(pos_accessions) <= 4) -> exon.five

intron %>%
  rowwise() %>%
  filter(find.len(pos_accessions) <= 4) -> intron.five

utr3 %>%
  rowwise() %>%
  filter(find.len(pos_accessions) <= 4) -> utr3.five

utr5 %>%
  rowwise() %>%
  filter(find.len(pos_accessions) <= 4) -> utr5.five

downstream %>%
  rowwise() %>%
  filter(find.len(pos_accessions) <= 4) -> downstream.five

## Convert expression values to ranks
ranks <- expression[0,]

for(i in 1:nrow(expression)) {
  ranks[i,] <- rank(expression[i,], ties.method = "first")
}

rownames(ranks) <- rownames(expression)


gatherData <- function(insertions, ranks) {
  data <- vector()
  for(i in 1:nrow(insertions)){
    row <- insertions[i,]
    gene <- as.character(row[[1]])
    accessions <- unlist(strsplit(gsub("-", "_", as.character(row[[3]])), ","))
    accessions <- intersect(accessions, colnames(ranks))
    if(length(accessions) > 0){
      for(acc in accessions){
        data <- c(data, as.numeric(ranks[gene,acc]))
      }
    }
  }
  return(data)
}

permuteLabels <- function(ranks, insertions) {
  data <- vector()
  all_names <- colnames(ranks)
  for(i in 1:nrow(insertions)){
    row <- insertions[i,]
    gene <- as.character(row[[1]])
    accessions <- unlist(strsplit(gsub("-", "_", as.character(row[[3]])), ","))
    accessions <- intersect(accessions, colnames(ranks))
    if(length(accessions) > 0){
      rand <- sample(all_names, length(accessions), replace = FALSE)
      for(acc in rand){
        data <- c(data, as.numeric(ranks[gene,acc]))
      }
    }
  }
  return(data)
}

exon.ranks <- gatherData(exon.five, ranks)
upstream.ranks <- gatherData(upstream.five, ranks)
downstream.ranks <- gatherData(downstream.five, ranks)
intron.ranks <- gatherData(intron.five, ranks)
utr5.ranks <- gatherData(utr5.five, ranks)
utr3.ranks <- gatherData(utr3.five, ranks)

exon.rand <- permuteLabels(ranks, exon.five)
upstream.rand <- permuteLabels(ranks, upstream.five)
downstream.rand <- permuteLabels(ranks, downstream.five)
utr5.rand <- permuteLabels(ranks, utr5.five)
utr3.rand <- permuteLabels(ranks, utr3.five)
intron.rand <- permuteLabels(ranks, intron.five)

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
  
  plot(real.hist$counts, pch=19, col = "blue", las=1, cex = 0.6,
       main = paste(title, " real"), ylim = y_range, xlab = "Expression bin", ylab = "Rare variants")
  lines(model.real$fitted.values)
  legend("top", paste("Rsq = ", rsq.real, ", p = ", p.real), bty = "n")
  
  plot(rand.hist$counts, pch=19, col = "blue", las = 1, cex = 0.6,
       main = paste(title, " random"), ylim = y_range, xlab = "Expression bin", ylab = "Rare variants")
  lines(model.rand$fitted.values)
  legend("top", paste("Rsq = ", rsq.rand, ", p = ", p.rand), bty = "n")
  
}

## Plots
pdf("../Plots/burden_gene_expression_dot_plots.pdf", width = 4, height = 4, useDingbats=FALSE)
dot_plot(exon.ranks, exon.rand, "Exon")
dot_plot(upstream.ranks, upstream.rand, "Upstream")
dot_plot(downstream.ranks, downstream.rand, "Downstream")
dot_plot(intron.ranks, intron.rand, "Intron")
dot_plot(utr3.ranks, utr3.rand, "3' UTR")
dot_plot(utr5.ranks, utr5.rand, "5' UTR")
dev.off()
