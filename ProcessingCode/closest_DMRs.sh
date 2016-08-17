gzip -dc ../RawData/TEPID_TEPAV.tsv.gz | sed 1d \
| bedtools closest -d -t first -b - -a ../ProcessedData/c_dmr_intersections.tsv.gz -wb \
| gzip - \
> ../ProcessedData/TE_C_DMR_distances.bed.gz

gzip -dc ../RawData/TEPID_TEPAV.tsv.gz | sed 1d \
| bedtools closest -d -t first -b - -a ../ProcessedData/cg_dmr_intersections.tsv.gz -wb \
| gzip - \
> ../ProcessedData/TE_CG_DMR_distances.bed.gz
