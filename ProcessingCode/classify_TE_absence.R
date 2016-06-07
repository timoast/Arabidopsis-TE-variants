library(dplyr)

poly <- read_tsv("../ProcessedData/poly_te_genotyped.tsv", col_names = F)

polarized <- poly %>%
  rowwise() %>%
  mutate(TE_present = length(unlist(strsplit(X5, ","))),
         TE_absent = length(unlist(strsplit(X6, ","))),
         Minor_allele = ifelse(TE_present < TE_absent, "presence", "absence"),
         Absence_classification = if(TE_absent < 43) {
           "True deletion"
         } else if(TE_absent > 172) {
           "No insertion"
         } else {
           "NA"
         },
         MAF = min(TE_present, TE_absent) / 216) %>%
  ungroup()

colnames(polarized) <- c("chromosome", "start", "end", "TE", "Accessions_TE_present", "Accessions_TE_absent", "Count_TE_present", "Count_TE_absent", "Minor_allele", "Absence_classification", "MAF")
write.table(polarized, file = "../RawData/tepav_regenerated.tsv",
            quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
