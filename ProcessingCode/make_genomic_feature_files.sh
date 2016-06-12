# This file will create coordinate files for gene upstream and downstream regions
# and gene introns, and place the output in ProcessedData
# It will then intersect the TE variants with all genomic features, including
# 2 kb upstream region
# 5' UTR
# Exon
# Intron
# 3' UTR
# 2 kb downstream region
# Pseudogenes
# Other Col-0 TEs
# DNase I hypersensitivity sites

# Information of the allele frequency or insertion classification (insertion or deletion)
# can be found in the ProcessedData/GeneFeatures/*_intersections.bed.gz files

# temporary decompressed TEPAV file
gunzip -c ../RawData/TEPID_TEPAV.tsv.gz | sed 1d > tepid_tepav.tsv

touch ../ProcessedData/gene_feature_counts.csv

# Upstream
echo Upstream
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") print $0}' - \
| awk 'BEGIN {FS=OFS="\t"} {$3 = "upstream"; $5 = $4; $4 = $4 - 2000; if ($4 < 0) $4 = 1; print $0}' - \
| gff2bed \
> ../ProcessedData/GeneFeatures/gene_upstream_regions.bed

bedtools intersect -a tepid_tepav.tsv -b ../ProcessedData/GeneFeatures/gene_upstream_regions.bed -wb \
> ../ProcessedData/GeneFeatures/gene_upstream_regions_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/gene_upstream_regions_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/gene_upstream_regions.bed \
>> ../ProcessedData/gene_feature_counts.csv

# Downstream
echo Downstream
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") print $0}' - \
| awk 'BEGIN {FS=OFS="\t"} {$3 = "downstream"; $4 = $5; $5 = $5 + 2000; print $0}' - \
| gff2bed \
> ../ProcessedData/GeneFeatures/gene_downstream_regions.bed

bedtools intersect -a tepid_tepav.tsv -b ../ProcessedData/GeneFeatures/gene_downstream_regions.bed -wb \
> ../ProcessedData/GeneFeatures/gene_downstream_regions_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/gene_downstream_regions_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/gene_downstream_regions.bed \
>> ../ProcessedData/gene_feature_counts.csv

# Exons
echo Exons
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "exon") print $0}' - \
| sort -k1,1 -k4,4n - \
> ../ProcessedData/GeneFeatures/exons.gff

# Introns
echo Introns
python add_feature_between.py ../ProcessedData/GeneFeatures/exons.gff \
| gff2bed \
| bedtools merge -c 10 -o distinct -i - \
> ../ProcessedData/GeneFeatures/introns.bed

bedtools intersect -wb -a tepid_tepav.tsv -b ../ProcessedData/GeneFeatures/introns.bed \
> ../ProcessedData/GeneFeatures/intron_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/intron_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/introns.bed \
>> ../ProcessedData/gene_feature_counts.csv

# Convert exons to bed
gff2bed < ../ProcessedData/GeneFeatures/exons.gff \
| bedtools merge -c 10 -o distinct -i - \
> ../ProcessedData/GeneFeatures/exons.bed
rm ../ProcessedData/GeneFeatures/exons.gff

# intersect exons with TEPAV
bedtools intersect -a tepid_tepav.tsv -b ../ProcessedData/GeneFeatures/exons.bed -wb \
> ../ProcessedData/GeneFeatures/exon_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/exon_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/exons.bed \
>> ../ProcessedData/gene_feature_counts.csv

# 5' UTR
echo UTR5
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "five_prime_UTR") print $0}' - \
| gff2bed \
| bedtools merge -c 10 -o distinct -i - \
> ../ProcessedData/GeneFeatures/utr5.bed

bedtools intersect -wb -a tepid_tepav.tsv -b ../ProcessedData/GeneFeatures/utr5.bed \
> ../ProcessedData/GeneFeatures/utr5_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/utr5_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/utr5.bed \
>> ../ProcessedData/gene_feature_counts.csv

# 3' UTR
echo UTR3
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "three_prime_UTR") print $0}' - \
| gff2bed \
| bedtools merge -c 10 -o distinct -i - \
> ../ProcessedData/GeneFeatures/utr3.bed

bedtools intersect -wb -a tepid_tepav.tsv -b ../ProcessedData/GeneFeatures/utr3.bed \
> ../ProcessedData/GeneFeatures/utr3_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/utr3_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/utr3.bed \
>> ../ProcessedData/gene_feature_counts.csv

# Pseudogene
echo Pseudogene
gunzip -c ../RawData/TAIR10_GFF3_genes.gff.gz \
| awk '{if ($3 == "pseudogene") print $0}' - \
| gff2bed \
> ../ProcessedData/GeneFeatures/pseudogene.bed

bedtools merge -c 10 -o distinct -i ../ProcessedData/GeneFeatures/pseudogene.bed \
| bedtools intersect -a tepid_tepav.tsv -b - -wb \
> ../ProcessedData/GeneFeatures/pseudogene_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/pseudogene_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
wc -l ../ProcessedData/GeneFeatures/pseudogene.bed \
>> ../ProcessedData/gene_feature_counts.csv

# Other TEs
echo TEs
gunzip -c ../RawData/TAIR9_TE.bed.gz \
| bedtools merge -c 5 -o distinct -i - \
| bedtools intersect -a tepid_tepav.tsv -b - -wb \
> ../ProcessedData/GeneFeatures/te_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/te_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
gunzip -c ../RawData/TAIR9_TE.bed.gz \
| wc -l \
>> ../ProcessedData/gene_feature_counts.csv

# DHS sites
gunzip -c ../RawData/Sullivan_DHS_PE_peaks_control.bed.gz \
| bedtools intersect -a tepid_tepav.tsv -b - -wb \
> ../ProcessedData/GeneFeatures/dhs_intersections.bed

# record number of intersections
wc -l ../ProcessedData/GeneFeatures/dhs_intersections.bed \
>> ../ProcessedData/gene_feature_counts.csv

# record total number
gunzip -c ../RawData/Sullivan_DHS_PE_peaks_control.bed.gz \
| wc -l \
>> ../ProcessedData/gene_feature_counts.csv

rm tepid_tepav.tsv
gzip ../ProcessedData/GeneFeatures/*