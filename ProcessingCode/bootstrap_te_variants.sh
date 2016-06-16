tepav=23095
av_size=501


for (( i=1; i <= 10000; i++ )); do
     # choose random coordinates
    python choose_random_coords.py $tepav $av_size \
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
done

bedtools closest -a ../ProcessedData/random_coords.bed -b ../RawData/c_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_c_dmr.bed.gz

bedtools closest -a ../ProcessedData/random_coords.bed -b ../RawData/cg_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/random_coords_closest_cg_dmr.bed.gz