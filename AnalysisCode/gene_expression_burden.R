library(dplyr)
library(readr)

set.seed(1)

expression <- read.table("../RawData/gene_expression.tsv.gz", header = T, row.names = 1)

gene_list <- read_tsv("../ProcessedData/gene_list.tsv", col_names = c("gene", "classification")) %>%
        filter(classification != "transposable_element_gene")

expression <- expression[rownames(expression) %in% gene_list$gene,]

upstream <- read_tsv("../ProcessedData/GeneFeatures/gene_upstream_regions_intersections.bed.gz", col_names = F) %>%
  select(X5:X12, X16) %>%
  rename(pos_accessions = X5, neg_accessions = X6, gene = X16)

downstream <- read_tsv("../ProcessedData/GeneFeatures/gene_downstream_regions_intersections.bed.gz", col_names = F) %>%
  select(X5:X12, X16) %>%
  rename(pos_accessions = X5, neg_accessions = X6, gene = X16)

exon <- read_tsv("../ProcessedData/GeneFeatures/exon_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6)

intron <- read_tsv("../ProcessedData/GeneFeatures/intron_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6)

utr5 <- read_tsv("../ProcessedData/GeneFeatures/utr5_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6)

utr3 <- read_tsv("../ProcessedData/GeneFeatures/utr3_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6)

upstream <- upstream[!duplicated(upstream),]	
exon <- exon[!duplicated(exon),]	
intron <- intron[!duplicated(intron),]	
utr5 <- utr5[!duplicated(utr5),]	
utr3 <- utr3[!duplicated(utr3),]	
downstream <- downstream[!duplicated(downstream),]	

find.len <- function(d) {
  return(length(unlist(strsplit(d, ","))))
}

upstream %>%
  rowwise() %>%
  filter(X12 == "Rare") -> upstream.rare

exon %>%
  rowwise() %>%
  filter(X12 == "Rare") -> exon.rare

intron %>%
  rowwise() %>%
  filter(X12 == "Rare") -> intron.rare

utr3 %>%
  rowwise() %>%
  filter(X12 == "Rare") -> utr3.rare

utr5 %>%
  rowwise() %>%
  filter(X12 == "Rare") -> utr5.rare

downstream %>%
  rowwise() %>%
  filter(X12 == "Rare") -> downstream.rare

## Convert expression values to ranks
ranks <- expression[0,]

for(i in 1:nrow(expression)) {
  ranks[i,] <- rank(expression[i,], ties.method = "first")
}
rownames(ranks) <- rownames(expression)

# save ranks
write_tsv(cbind(gene=rownames(ranks), ranks), "../ProcessedData/expression_ranks.tsv")
system("gzip ../ProcessedData/expression_ranks.tsv")

gatherData <- function(tepav, indel, ranks) {  # set indel as either pos_accessions or neg_accessions
  data <- vector()
  for(i in 1:nrow(tepav)){
    row <- tepav[i,]
    gene <- as.character(row[["gene"]])
    accessions <- unlist(strsplit(gsub("-", "_", as.character(row[[indel]])), ","))
    accessions <- intersect(accessions, colnames(ranks))
    if(length(accessions) > 0){
      for(acc in accessions){
        data <- c(data, as.numeric(ranks[gene, acc]))
      }
    }
  }
  return(data)
}

permuteLabels <- function(ranks, tepav, indel) {
  data <- vector()
  all_names <- colnames(ranks)
  for(i in 1:nrow(tepav)){
    row <- tepav[i,]
    gene <- as.character(row[["gene"]])
    accessions <- unlist(strsplit(gsub("-", "_", as.character(row[[indel]])), ","))
    accessions <- intersect(accessions, colnames(ranks))
    if(length(accessions) > 0){
      rand <- sample(all_names, length(accessions), replace = FALSE)
      for(acc in rand){
        data <- c(data, as.numeric(ranks[gene, acc]))
      }
    }
  }
  return(data)
}

exon.ranks.insertion <- gatherData(filter(exon.rare, X10 == "No insertion"), "pos_accessions", ranks)
exon.ranks.deletion <- gatherData(filter(exon.rare, X10 == "True deletion"), "neg_accessions", ranks)

upstream.ranks.insertion <- gatherData(filter(upstream.rare, X10 == "No insertion"), "pos_accessions", ranks)
upstream.ranks.deletion <- gatherData(filter(upstream.rare, X10 == "True deletion"), "neg_accessions", ranks)

downstream.ranks.insertion <- gatherData(filter(downstream.rare, X10 == "No insertion"), "pos_accessions", ranks)
downstream.ranks.deletion <- gatherData(filter(downstream.rare, X10 == "True deletion"), "neg_accessions", ranks)

intron.ranks.insertion <- gatherData(filter(intron.rare, X10 == "No insertion"), "pos_accessions", ranks)
intron.ranks.deletion <- gatherData(filter(intron.rare, X10 == "True deletion"), "neg_accessions", ranks)

utr5.ranks.insertion <- gatherData(filter(utr5.rare, X10 == "No insertion"), "pos_accessions", ranks)
utr5.ranks.deletion <- gatherData(filter(utr5.rare, X10 == "True deletion"), "neg_accessions", ranks)

utr3.ranks.insertion <- gatherData(filter(utr3.rare, X10 == "No insertion"), "pos_accessions", ranks)
utr3.ranks.deletion <- gatherData(filter(utr3.rare, X10 == "True deletion"), "neg_accessions", ranks)

exon.rand <- permuteLabels(ranks, filter(exon.rare, X10 == "No insertion"), "pos_accessions")
upstream.rand <- permuteLabels(ranks, filter(upstream.rare, X10 == "No insertion"), "pos_accessions")
downstream.rand <- permuteLabels(ranks, filter(downstream.rare, X10 == "No insertion"), "pos_accessions")
utr5.rand <- permuteLabels(ranks, filter(utr5.rare, X10 == "No insertion"), "pos_accessions")
utr3.rand <- permuteLabels(ranks, filter(utr3.rare, X10 == "No insertion"), "pos_accessions")
intron.rand <- permuteLabels(ranks, filter(intron.rare, X10 == "No insertion"), "pos_accessions")

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

dot_plot <- function(insertions, deletions, random, title, brk = 72) {
  
  insertion.hist <- hist(insertions, breaks = brk, plot = F)
  deletion.hist <- hist(deletions, breaks = brk, plot = F)
  combined.hist <- hist(c(insertions, deletions), breaks = brk, plot = F)
  rand.hist <- hist(random, breaks = brk, plot = F)
  
  y_range <- c(min(insertion.hist$counts) - 10, max(insertion.hist$counts) + 10)
  if(y_range[1] < 0) { y_range[1] = 0 }
  
  model.insertion <- lm(insertion.hist$counts ~ poly(insertion.hist$breaks[2:length(insertion.hist$breaks)], 2, raw=TRUE))
  rsq.insertion <- signif(summary(model.insertion)$adj.r.squared, 3)
  p.insertion <- signif(lmp(model.insertion), 3)
  
  model.deletion <- lm(deletion.hist$counts ~ poly(deletion.hist$breaks[2:length(deletion.hist$breaks)], 2, raw=TRUE))
  rsq.deletion <- signif(summary(model.deletion)$adj.r.squared, 3)
  p.deletion <- signif(lmp(model.deletion), 3)
  
  model.combined <- lm(combined.hist$counts ~ poly(combined.hist$breaks[2:length(combined.hist$breaks)], 2, raw=TRUE))
  rsq.combined <- signif(summary(model.combined)$adj.r.squared, 3)
  p.combined <- signif(lmp(model.combined), 3)
  
  model.rand <- lm(rand.hist$counts ~ poly(rand.hist$breaks[2:length(rand.hist$breaks)], 2, raw=TRUE))
  rsq.rand <- signif(summary(model.rand)$adj.r.squared, 3)
  p.rand <- signif(lmp(model.rand), 3)
  
  plot(insertion.hist$counts, pch=19, col = "blue", las=1, cex = 0.6,
       main = paste(title, " TE insertions"), ylim = y_range, xlab = "Expression bin", ylab = "Rare variants")
  lines(model.insertion$fitted.values)
  legend("top", paste("Rsq = ", rsq.insertion, ", p = ", p.insertion), bty = "n")

  plot(deletion.hist$counts, pch=19, col = "blue", las=1, cex = 0.6,
       main = paste(title, " TE deletion"), xlab = "Expression bin", ylab = "Rare variants")
  lines(model.deletion$fitted.values)
  legend("top", paste("Rsq = ", rsq.deletion, ", p = ", p.deletion), bty = "n")

  plot(combined.hist$counts, pch=19, col = "blue", las=1, cex = 0.6,
       main = paste(title, " combined \n insertions and deletion"), xlab = "Expression bin", ylab = "Rare variants")
  lines(model.combined$fitted.values)
  legend("top", paste("Rsq = ", rsq.combined, ", p = ", p.combined), bty = "n")
  
  plot(rand.hist$counts, pch=19, col = "blue", las = 1, cex = 0.6,
       main = paste(title, " random"), ylim = y_range, xlab = "Expression bin", ylab = "Rare variants")
  lines(model.rand$fitted.values)
  legend("top", paste("Rsq = ", rsq.rand, ", p = ", p.rand), bty = "n")
  
}

## Plots
pdf("../Plots/burden_gene_expression_dot_plots.pdf", width = 10, height = 15, useDingbats=FALSE)
par(mfrow=c(6,4))
dot_plot(upstream.ranks.insertion, upstream.ranks.deletion, upstream.rand, "Upstream")
dot_plot(utr5.ranks.insertion, utr5.ranks.deletion, utr5.rand, "5' UTR")
dot_plot(exon.ranks.insertion, exon.ranks.deletion, exon.rand, "Exon")
dot_plot(intron.ranks.insertion, intron.ranks.deletion, intron.rand, "Intron")
dot_plot(utr3.ranks.insertion, utr3.ranks.deletion, utr3.rand, "3' UTR")
dot_plot(downstream.ranks.insertion, downstream.ranks.deletion, downstream.rand, "Downstream")
dev.off()
