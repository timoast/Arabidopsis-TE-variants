library(dplyr)
library(ggplot2)
library(qvalue)
library(readr)

expression <- read_tsv("../RawData/gene_expression.tsv.gz", col_names = TRUE)

gene_list <- read_tsv("../ProcessedData/gene_list.tsv", col_names = c("gene", "classification")) %>%
        filter(classification != "transposable_element_gene")

expression <- inner_join(expression, gene_list)

upstream <- read_tsv("../ProcessedData/GeneFeatures/gene_upstream_regions_intersections.bed.gz", col_names = F) %>%
  select(X5:X12, X16) %>%
  rename(pos_accessions = X5, neg_accessions = X6, AbsenceClassification = X10, gene = X16)

downstream <- read_tsv("../ProcessedData/GeneFeatures/gene_downstream_regions_intersections.bed.gz", col_names = F) %>%
  select(X5:X12, X16) %>%
  rename(pos_accessions = X5, neg_accessions = X6, AbsenceClassification = X10, gene = X16)

exon <- read_tsv("../ProcessedData/GeneFeatures/exon_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6, AbsenceClassification = X10)

intron <- read_tsv("../ProcessedData/GeneFeatures/intron_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6, AbsenceClassification = X10)

utr5 <- read_tsv("../ProcessedData/GeneFeatures/utr5_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6, AbsenceClassification = X10)

utr3 <- read_tsv("../ProcessedData/GeneFeatures/utr3_intersections.bed.gz", col_names = F) %>%
  mutate(gene = substr(X16, 8, 16)) %>%
  select(X5:X12, gene) %>%
  rename(pos_accessions = X5, neg_accessions = X6, AbsenceClassification = X10)

upstream <- upstream[!duplicated(upstream),]	
exon <- exon[!duplicated(exon),]	
intron <- intron[!duplicated(intron),]	
utr5 <- utr5[!duplicated(utr5),]	
utr3 <- utr3[!duplicated(utr3),]	
downstream <- downstream[!duplicated(downstream),]	
	
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
    x <- c(x, as.numeric(exp[exp$gene == gene, acc]))	
  }	
  return(x)	
}	
	
test.sig.cor <- function(genes, expression){	
  data <- matrix(nrow = nrow(genes), ncol = 5)	
  # gene, te, acc names with, acc names without, p-val, cor	
  for(i in 1:nrow(genes)){	
    row <- genes[i,]
    gene <- as.character(row[["gene"]])
    absence <- row[["AbsenceClassification"]]
    ins.acc <- unlist(strsplit(gsub("-", "_", row[["pos_accessions"]]), ","))	
    no.ins.acc <- unlist(strsplit(gsub("-", "_", row[["neg_accessions"]]), ","))	
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
      data[i,] <- c(pval, r, lfc, gene, absence)	
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
colnames(all.data) <- c("pvalue", "cor", "log2FoldChange", "gene", "AbsenceClassification", "feature", "sample_size")	
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
  ggsave("../Plots/volcano_plots.pdf", width = 20, height = 4, units = "cm", useDingbats=F)	

ggplot(filter(qval.data, AbsenceClassification == "True deletion"), aes(log2FoldChange, -log10(qvalue)), label = gene) +	
  geom_point(aes(color=sig, size=0.1)) + facet_wrap(~feature, nrow=1) + theme_bw(base_size = 8)+	
  geom_hline(yintercept = -log10(sig_threshold), size = 0.1) +	
  geom_vline(xintercept = -fc, size = 0.1) +	
  geom_vline(xintercept = fc, size = 0.1) +	
  scale_size_identity() +	
  scale_color_manual(values = c("blue", "grey", "red")) +	
  theme(legend.position="none") +	
  xlim(-6, 6) +
  ggsave("../Plots/volcano_plots_deletions.pdf", width = 20, height = 4, units = "cm", useDingbats=F)	


ggplot(filter(qval.data, AbsenceClassification == "No insertion"), aes(log2FoldChange, -log10(qvalue)), label = gene) +	
  geom_point(aes(color=sig, size=0.1)) + facet_wrap(~feature, nrow=1) + theme_bw(base_size = 8)+	
  geom_hline(yintercept = -log10(sig_threshold), size = 0.1) +	
  geom_vline(xintercept = -fc, size = 0.1) +	
  geom_vline(xintercept = fc, size = 0.1) +	
  scale_size_identity() +	
  scale_color_manual(values = c("blue", "grey", "red")) +	
  theme(legend.position="none") +	
  xlim(-6, 6) +
  ggsave("../Plots/volcano_plots_insertions.pdf", width = 20, height = 4, units = "cm", useDingbats=F)	


ggplot(filter(qval.data, is.na(AbsenceClassification)), aes(log2FoldChange, -log10(qvalue)), label = gene) +	
  geom_point(aes(color=sig, size=0.1)) + facet_wrap(~feature, nrow=1) + theme_bw(base_size = 8)+	
  geom_hline(yintercept = -log10(sig_threshold), size = 0.1) +	
  geom_vline(xintercept = -fc, size = 0.1) +	
  geom_vline(xintercept = fc, size = 0.1) +	
  scale_size_identity() +	
  scale_color_manual(values = c("blue", "grey", "red")) +	
  theme(legend.position="none") +	
  xlim(-6, 6) +
  ggsave("../Plots/volcano_plots_na_absence_call.pdf", width = 20, height = 4, units = "cm", useDingbats=F)	

qval.data %>%
  filter(sig != "NS") %>%
  write.table(., file="../ProcessedData/sig_genes.tsv", quote = F, sep = "\t", row.names = F)

print("Number of genes in each feature:")
qval.data %>%
  group_by(feature) %>%
  summarize(count = n())

print("Number of significantly differentially expressed genes:")
sig <- filter(qval.data, sig != "NS")
length(unique(sig$gene))
