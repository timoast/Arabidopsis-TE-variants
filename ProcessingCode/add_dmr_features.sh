# Alter DMR files with the genomic features that they intersect

bedtools intersect -a ../RawData/c_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz -wb \
| cut -f1,2,3,7 \
> c_dmr_intersections.bed

bedtools intersect -a ../RawData/c_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz -v \
| awk 'BEGIN {FS=OFS="\t"} {print $0, "intergenic"}' \
| cat c_dmr_intersections.bed - \
| sort -k1,1 -k2,2n \
| gzip \
> ../ProcessedData/c_dmr_intersections.tsv.gz

rm c_dmr_intersections.bed

bedtools intersect -a ../RawData/cg_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz -wb \
| cut -f1,2,3,7 \
> cg_dmr_intersections.bed

bedtools intersect -a ../RawData/cg_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz -v \
| awk 'BEGIN {FS=OFS="\t"} {print $0, "intergenic"}' \
| cat cg_dmr_intersections.bed - \
| sort -k1,1 -k2,2n \
| gzip \
> ../ProcessedData/cg_dmr_intersections.tsv.gz

rm cg_dmr_intersections.bed