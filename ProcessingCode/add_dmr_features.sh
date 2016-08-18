# Alter DMR files with the genomic features that they intersect
# Add distance from each DMR to the closest TE

bedtools intersect -wa -wb -a ../RawData/c_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz \
    | cut -f1,2,3,7 \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 4 -o distinct \
	  > c_dmr_intersections.bed

bedtools intersect -wa -a ../RawData/c_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz -v \
    | awk 'BEGIN {FS=OFS="\t"} {print $0, "intergenic"}' \
    | cat c_dmr_intersections.bed - \
    | sort -k1,1 -k2,2n \
    | gzip \
	  > ../ProcessedData/c_dmr_intersections.tsv.gz

rm c_dmr_intersections.bed

bedtools intersect -wa -wb -a ../RawData/cg_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz \
    | cut -f1,2,3,7 \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 4 -o distinct \
	  > cg_dmr_intersections.bed

bedtools intersect -wa -a ../RawData/cg_dmrs.tsv.gz -b ../ProcessedData/GeneFeatures/genomic_features.bed.gz -v \
    | awk 'BEGIN {FS=OFS="\t"} {print $0, "intergenic"}' \
    | cat cg_dmr_intersections.bed - \
    | sort -k1,1 -k2,2n \
    | gzip \
	  > ../ProcessedData/cg_dmr_intersections.tsv.gz

rm cg_dmr_intersections.bed

# distance from each C-DMR to closest TE
bedtools closest -d -t first -a ../RawData/c_dmrs.tsv.gz \
         -b ../RawData/TAIR9_TE.bed.gz \
    | cut -f 1,2,3,8,11 \
    | gzip \
	  > ../ProcessedData/c_dmrs_closest_te.tsv.gz 

# distance from each C-DMR to closest TE variant
gzip -dc ../RawData/TEPID_TEPAV.tsv.gz \
    | sed 1d \
    | bedtools closest -d -t first -a ../RawData/c_dmrs.tsv.gz \
	       -b - \
    | cut -f 1,2,3,17 \
    | gzip \
	  > ../ProcessedData/c_dmr_tepav_distance.tsv.gz
