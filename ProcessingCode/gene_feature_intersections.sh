# This script will intersect the TE variants with all genomic features, including
# 2 kb upstream region
# 5' UTR
# Exon
# Intron
# 3' UTR
# 2 kb downstream region
# Pseudogenes
# Other Col-0 TEs
# DNase I hypersensitivity sites

gunzip -c ../RawData/TEPID_TEPAV.tsv.gz | sed 1d > tepid_tepav.tsv

# need to remove redundancy in gene models first
gunzip -c ../ProcessedData/GeneFeatures/exons.gff.gz \
| bedtools intersect -a - -b tepid_tepav.tsv -wb \
> ../ProcessedData/GeneFeatures/exon_intersections.gff

rm tepid_tepav.tsv