#! /bin/sh

# to get mean length and number of each TE variant type (in R):
# library(readr)
# library(dplyr)
# 
# read_tsv("../RawData/TEPID_TEPAV.tsv.gz", col_names = T) %>%
#   mutate(l = end-start) %>%
#   group_by(Absence_classification) %>%
#   summarise(ml = mean(l), c = n())

all=23095
insertions=15077
deletions=5856
inter=2162

all_size=501
ins_size=134
del_size=1170
inter_size=1249

for (( i=1; i <= 10000; i++ )); do
    ### All
    # choose random coordinates
    python choose_random_coords.py $all $all_size \
    | sort -k1,1 -k2,2n - \
    > ../ProcessedData/random_coords.bed

    # c-dmrs
    bedtools closest -a ../ProcessedData/random_coords.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_c_dmr_rep_10k.txt

    #cg-dmrs
    bedtools closest -a ../ProcessedData/random_coords.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_cg_dmr_rep_10k.txt

    ### insertions
    # choose random coordinates
    python choose_random_coords.py $insertions $ins_size \
    | sort -k1,1 -k2,2n - \
    > ../ProcessedData/random_coords_insertions.bed

    # c-dmrs
    bedtools closest -a ../ProcessedData/random_coords_insertions.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_c_dmr_rep_10k_insertions.txt

    #cg-dmrs
    bedtools closest -a ../ProcessedData/random_coords_insertions.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_cg_dmr_rep_10k_insertions.txt

    ### deletions
    # choose random coordinates
    python choose_random_coords.py $deletions $del_size \
    | sort -k1,1 -k2,2n - \
    > ../ProcessedData/random_coords_deletions.bed

    # c-dmrs
    bedtools closest -a ../ProcessedData/random_coords_deletions.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_c_dmr_rep_10k_deletions.txt

    #cg-dmrs
    bedtools closest -a ../ProcessedData/random_coords_deletions.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_cg_dmr_rep_10k_deletions.txt

    ### intermediate frequency
    # choose random coordinates
    python choose_random_coords.py $inter $inter_size \
    | sort -k1,1 -k2,2n - \
    > ../ProcessedData/random_coords_intermediate.bed

    # c-dmrs
    bedtools closest -a ../ProcessedData/random_coords_intermediate.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_c_dmr_rep_10k_intermediate.txt

    #cg-dmrs
    bedtools closest -a ../ProcessedData/random_coords_intermediate.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
    | awk 'BEGIN {FS=OFS="\t"} {if ($7 < 1000) print $0}' - \
    | wc -l \
    >> ../ProcessedData/random_selections_cg_dmr_rep_10k_intermediate.txt
done

bedtools closest -a ../ProcessedData/random_coords.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_c_dmr.bed.gz

bedtools closest -a ../ProcessedData/random_coords.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_cg_dmr.bed.gz

bedtools closest -a ../ProcessedData/random_coords_insertions.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_c_dmr_insertions.bed.gz

bedtools closest -a ../ProcessedData/random_coords_insertions.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_cg_dmr_insertions.bed.gz

bedtools closest -a ../ProcessedData/random_coords_deletions.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_c_dmr_deletions.bed.gz

bedtools closest -a ../ProcessedData/random_coords_deletions.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_cg_dmr_deletions.bed.gz

bedtools closest -a ../ProcessedData/random_coords_intermediate.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_c_dmr_intermediate.bed.gz

bedtools closest -a ../ProcessedData/random_coords_intermediate.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_cg_dmr_intermediate.bed.gz