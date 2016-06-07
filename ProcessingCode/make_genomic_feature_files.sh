# This file will create coordinate files for gene upstream and downstream regions
# and gene introns, and place the output in ProcessedData

# Upstream
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") print $0}' - \
| awk 'BEGIN {FS=OFS="\t"} {$3 = "upstream"; $5 = $4; $4 = $4 - 2000; print $0}' - \
| gzip - \
> ../ProcessedData/gene_upstream_regions.gff.gz

# Downstream
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") print $0}' - \
| awk 'BEGIN {FS=OFS="\t"} {$3 = "downstream"; $4 = $5; $5 = $5 + 2000; print $0}' - \
| gzip - \
> ../ProcessedData/gene_downstream_regions.gff.gz

# Introns
wget ftp://ftp.ensemblgenomes.org//pub/release-31/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.31.gff3.gz

gt gff3 -addintrons Arabidopsis_thaliana.TAIR10.31.gff3.gz \
| grep "intron" - \
| awk 'BEGIN {FS=OFS="\t"} {$1 = "chr"$1; print $0}' - \
| bedtools intersect -a - -b ../RawData/TAIR10_GFF3_genes.gff.gz -wb \
| awk 'BEGIN {FS=OFS="\t"} {if ($12 == "gene") print $10, $11, $3, $4, $5, $15, $16, $17, $18}' - \
| gzip - \
> ../ProcessedData/introns.gff.gz