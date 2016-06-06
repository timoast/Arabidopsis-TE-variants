library(dplyr)	
library(ggplot2)	
library(qvalue)	

data_path = "/Users/timstuart/Dropbox/PhD/Projects/Population_epigenetics/RNA/"	

# this data only has accession names with RNA data available	
expression <- read.table(paste(data_path, "expression.tsv", sep = ""), header = T, row.names = 1)	
upstream <- read.table(paste(data_path, 'upstream_insertion.tsv', sep = ""), header = F)	
exon <- read.table(paste(data_path, 'exon_insertion.tsv', sep = ""), header = F)	
intron <- read.table(paste(data_path, "intron_insertion.tsv", sep = ""), header = F)	
utr5 <- read.table(paste(data_path, "utr5_insertion.tsv", sep = ""), header = F)	
utr3 <- read.table(paste(data_path, "utr3_insertion.tsv", sep = ""), header = F)	
downstream <- read.table(paste(data_path, "downstream_insertion.tsv", sep = ""), header = F)	

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
	
# filter where insertion in at least 5 accessions (MAF >3%)	
upstream %>%	
  rowwise() %>%	
  filter(find.len(pos_accessions) > 4) -> upstream.five	
	
exon %>%	
  rowwise() %>%	
  filter(find.len(pos_accessions) > 4) -> exon.five	
	
intron %>%	
  rowwise() %>%	
  filter(find.len(pos_accessions) > 4) -> intron.five	
	
utr3 %>%	
  rowwise() %>%	
  filter(find.len(pos_accessions) > 4) -> utr3.five	
	
utr5 %>%	
  rowwise() %>%	
  filter(find.len(pos_accessions) > 4) -> utr5.five	
	
downstream %>%	
  rowwise() %>%	
  filter(find.len(pos_accessions) > 4) -> downstream.five	
 	
lookupExp <- function(accessions, gene, exp) {	
  x <- vector()	
  for(acc in accessions) {	
    x <- c(x, as.numeric(exp[gene,acc]))	
  }	
  return(x)	
}	
	
test.sig.cor <- function(genes, expression){	
  data <- matrix(nrow = nrow(genes), ncol = 4)	
  # gene, te, acc names with, acc names without, p-val, cor	
  for(i in 1:nrow(genes)){	
    row <- genes[i,]	
    gene <- as.character(row[[1]])	
    te <- as.character(row[[2]])	
    ins.acc <- unlist(strsplit(gsub("-", "_", row[[3]]), ","))	
    no.ins.acc <- unlist(strsplit(gsub("-", "_", row[[4]]), ","))	
    ins.acc <- intersect(ins.acc, colnames(expression))	
    no.ins.acc <- intersect(no.ins.acc, colnames(expression))	
    ins.data <- lookupExp(ins.acc, gene, expression)	
    ins.data <- ins.data[is.finite(ins.data)]	
    no.ins.data <- lookupExp(no.ins.acc, gene, expression)	
    no.ins.data <- no.ins.data[is.finite(no.ins.data)]	
    if(length(ins.data) > 5 && length(no.ins.data) > 5) {	
      pval <- wilcox.test(ins.data, no.ins.data)[[3]]	
      r <- cor(c(ins.data, no.ins.data),	
          c(rep(1, length(ins.data)),	
            rep(2, length(no.ins.data))),	
          use="pairwise.complete.obs")	
      lfc <- log2(mean(ins.data)/mean(no.ins.data))	
      data[i,] <- c(pval, r, lfc, gene)	
    }	
  }	
  return(data)	
}	

upstream.data <- test.sig.cor(upstream.five, expression)	
exon.data <- test.sig.cor(exon.five, expression)	
intron.data <- test.sig.cor(intron.five, expression)	
utr3.data <- test.sig.cor(utr3.five, expression)	
utr5.data <- test.sig.cor(utr5.five, expression)	
downstream.data <- test.sig.cor(downstream.five, expression)	
	
# make one dataframe	
exon.df <- as.data.frame(exon.data) %>% mutate(feature = "Exon", sample_size = nrow(exon.data))	
intron.df <- as.data.frame(intron.data) %>% mutate(feature = "Intron", sample_size = nrow(intron.data))	
utr3.df <- as.data.frame(utr3.data) %>% mutate(feature = "3' UTR", sample_size = nrow(utr3.data))	
utr5.df <- as.data.frame(utr5.data) %>% mutate(feature = "5' UTR", sample_size = nrow(utr5.data))	
upstream.df <- as.data.frame(upstream.data) %>% mutate(feature = "Upstream", sample_size = nrow(upstream.data))	
downstream.df <- as.data.frame(downstream.data) %>% mutate(feature = "Downstream", sample_size = nrow(downstream.data))	
	
all.data <- bind_rows(exon.df, intron.df, utr3.df, utr5.df, upstream.df, downstream.df)	
colnames(all.data) <- c("pvalue", "cor", "log2FoldChange", "gene", "feature", "sample_size")	
all.data <- transform(all.data, pvalue = as.numeric(pvalue))	
all.data <- transform(all.data, log2FoldChange = as.numeric(log2FoldChange))	
 	
q <- qvalue(all.data$pvalue, lambda = seq(0, 0.6, 0.05), smooth.df = 4)	
qvals <- as.data.frame(q$qvalues)	
qval.data <- bind_cols(all.data, qvalue = qvals)	
qval.data <- rename(qval.data, qvalue = `q$qvalues`)	
 	
sig_threshold <- 0.01	
fc <- log2(2)	
qval.data %>%	
  rowwise() %>%	
  filter(!is.na(pvalue)) %>%	
  mutate(sig = ifelse(qvalue<sig_threshold & log2FoldChange > fc, "Up",	
                      ifelse(qvalue<sig_threshold & log2FoldChange < -fc, "Down", "NS")))-> qval.data	
	
qval.data$feature <- factor(qval.data$feature, levels = c("Upstream", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream"))	
	
ggplot(qval.data, aes(log2FoldChange, -log10(qvalue)), label = gene) +	
  geom_point(aes(color=sig, size=0.1)) + facet_wrap(~feature, nrow=1) + theme_bw(base_size = 8)+	
  geom_hline(yintercept = -log10(sig_threshold), size = 0.1) +	
  geom_vline(xintercept = -fc, size = 0.1) +	
  geom_vline(xintercept = fc, size = 0.1) +	
  scale_size_identity() +	
  scale_color_manual(values = c("blue", "grey", "red")) +	
  theme(legend.position="none") +	
  xlim(-6, 6)	+
  ggsave("volcano_plots.pdf", width = 20, height = 4, units = "cm", useDingbats=F)	
	
qval.data %>%
  filter(sig != "NS") %>%
  write.table(., file="sig_genes.tsv", quote = F, sep = "\t", row.names = F)