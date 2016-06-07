# This file will create coordinate files for gene upstream and downstream regions
# and gene introns, and place the output in ProcessedData

# Upstream
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") print $0}' - \
| awk 'BEGIN {FS=OFS="\t"} {$3 = "upstream"; $5 = $4; $4 = $4 - 2000; print $0}' - \
| gzip - \
> ../ProcessedData/GeneFeatures/gene_upstream_regions.gff.gz

# Downstream
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") print $0}' - \
| awk 'BEGIN {FS=OFS="\t"} {$3 = "downstream"; $4 = $5; $5 = $5 + 2000; print $0}' - \
| gzip - \
> ../ProcessedData/GeneFeatures/gene_downstream_regions.gff.gz

# Exons
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "exon") print $0}' - \
| sort -k1,1 -k4,4n - \
> ../ProcessedData/GeneFeatures/exons.gff

# Introns
python add_feature_between.py ../ProcessedData/GeneFeatures/exons.gff \
| gzip - \
> ../ProcessedData/GeneFeatures/introns.gff.gz

# 5' UTR
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "five_prime_UTR") print $0}' - \
| gzip - \
> ../ProcessedData/GeneFeatures/utr5.gff.gz

# 3' UTR
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "three_prime_UTR") print $0}' - \
| gzip - \
> ../ProcessedData/GeneFeatures/utr3.gff.gz

# Pseudogene
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "pseudogene") print $0}' - \
| gzip - \
> ../ProcessedData/GeneFeatures/pseudogene.gff.gz

gff2bed < ../ProcessedData/GeneFeatures/exons.gff \
| bedtools merge -c 10 -o distinct -i - \
| gzip - \
> ../ProcessedData/GeneFeatures/exons.bed.gz
rm ../ProcessedData/GeneFeatures/exons.gff