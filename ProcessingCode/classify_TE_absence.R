library(dplyr)

poly <- read_tsv("../ProcessedData/poly_te_genotyped.tsv.gz", col_names = F)

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
         MAF = min(TE_present, TE_absent) / 216,
         Frequency_classification = ifelse(MAF < 0.03, "Rare", "Common")) %>%
  ungroup()

colnames(polarized) <- c("chromosome", "start", "end", "TE", "Accessions_TE_present", "Accessions_TE_absent", "Count_TE_present", "Count_TE_absent", "Minor_allele", "Absence_classification", "MAF", "Frequency_classification")
# add the LD data
polarized <- mutate(polarized, id = paste(chromosome, start, TE))

ld <- read_tsv("../RawData/snp_association.tsv.gz") %>%
  rowwise() %>%
  mutate(id = paste(chrom, pos, TEID), LD = keys[[flag]]) %>%
  select(id, LD)

data <- left_join(polarized, ld, by = "id") %>%  select(-id)

write_tsv(data,  "../RawData/tepav_regenerated.tsv")
system("gzip ../RawData/tepav_regenerated.tsv")