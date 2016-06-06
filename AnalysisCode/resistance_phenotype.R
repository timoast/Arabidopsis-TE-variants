library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

download.file("http://journals.plos.org/plosgenetics/article/asset?unique&id=info:doi/10.1371/journal.pgen.0010060.sd002",
              destfile = "../RawData/resistance_dataset.csv") 
data <- read_csv("../RawData/resistance_dataset.csv", col_names = FALSE)
tepav <- read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = TRUE)

te.var <- tepav %>%
  filter(TE == "AT3TE29460")

pos.acc <- te.var[2,5][[1]]
pos <- strsplit(pos.acc, ",")[[1]]

names <- tepav[1,5:6]
accessions <- strsplit(names[[2]], ",")
all_schmitz <- gsub("_", "-", accessions[[1]])

data %>%
  filter(X1 %in% all_schmitz) %>%
  mutate(TE = ifelse(X1 %in% pos, TRUE, FALSE)) -> filtered_data

filtered_data <- transform(filtered_data, X6 = as.numeric(X6))
filtered_data <- transform(filtered_data, X7 = as.numeric(X7))
filtered_data <- transform(filtered_data, X8 = as.numeric(X8))
filtered_data <- transform(filtered_data, X9 = as.numeric(X9))

filtered_data %>%
  group_by(TE) %>%
  summarise(avrPpH3 = mean(X6, na.rm = T) * 100, avrRpm1 = mean(X7, na.rm = T) * 100,
            avrRpt2 = mean(X8, na.rm = T) * 100, avrB = mean(X9, na.rm = T) * 100, sample_size = n()) -> stats

reshaped <- head(melt(stats), -2)

ggplot(reshaped, aes(variable, value, TE)) +
  geom_bar(stat = "identity", aes(fill = TE), position = "dodge", color = "black", size = 0.2) + theme_bw() +
  ylim(0, 100) + ylab("Accessions resistant (%)") + xlab("avr genes") +
  scale_fill_brewer("Blues") + ggtitle("P. syringae resistance") + 
  theme(legend.position = "null") +
  ggsave(height=10,width=10, units = "cm", path="../Plots", filename="resistance.pdf", useDingbats=FALSE)