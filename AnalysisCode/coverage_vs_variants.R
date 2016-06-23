library(dplyr)
library(readr)
library(tidyr)

cov <- read_tsv("../ProcessedData/coverages.tsv", col_names = c("Accession", "Coverage"))
ins <- read_tsv("../ProcessedData/insertion_counts.tsv", col_names = c("Accession", "Insertions"))
del <- read_tsv("../ProcessedData/deletion_counts.tsv", col_names = c("Accession", "Absences"))

dat <- left_join(cov, ins) %>% left_join(del) %>%
  gather(indel, Count, Insertions, Absences)

model.del <- lm(Count ~ Coverage, data=filter(dat, indel == "Absences"))
rsq.del <- signif(summary(model.del)$adj.r.squared, 3)

model.ins <- lm(Count ~ Coverage, data=filter(dat, indel == "Insertions"))
rsq.ins <- signif(summary(model.ins)$adj.r.squared, 3)

pdf("../Plots/coverage_variant_count.pdf", height=4, width = 6.5, useDingbats = F)
par(mfrow=c(1,2))
plot(Count ~ Coverage, data=filter(dat, indel == "Absences"), ylim=c(0, 2000), cex=0.8, las=2)
abline(model.del)
title("TE absences")
legend("top", paste("Rsq = ", rsq.del), bty="n")

plot(Count ~ Coverage, data=filter(dat, indel == "Insertions"), ylim=c(0, 800), cex=0.8, las=2)
abline(model.ins)
title("TE insertions")
legend("top", paste("Rsq = ", rsq.ins), bty = "n")
dev.off()