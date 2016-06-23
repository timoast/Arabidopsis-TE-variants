library(dplyr)
library(readr)

random <- read_tsv(file = "../ProcessedData/kee_random_intersections.bed", col_names = F)
kee <- read_tsv(file = "../ProcessedData/KEE_intersection_count.tsv",
                col_names = c("chr", "start", "stop", "KEE", "TE variants"))

pdf("../Plots/KEE_region_intersections.pdf", height = 4, width = 6)
hist(random$X4, breaks = 80, col = "lightblue",
     xlab = "TE insertions",
     main = "Frequency of TE insertions\nin random 300 kb regions",
     xlim = c(0, 300))
abline(v = kee$`TE variants`)
dev.off()

kee %>%
  rowwise() %>%
  mutate(pval = table(random$X4 > `TE variants`)["TRUE"][[1]] / 10000) %>%
  write_tsv(., "../ProcessedData/kee_pvals.tsv")