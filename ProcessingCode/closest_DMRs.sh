gunzip -c ../RawData/TEPID_TEPAV.tsv.gz | sed 1d \
| bedtools closest -b - -a ../RawData/c_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/TE_C_DMR_distances.bed.gz

gunzip -c ../RawData/TEPID_TEPAV.tsv.gz | sed 1d \
| bedtools closest -b - -a ../RawData/cg_dmrs.tsv.gz -d -wb \
| gzip - \
> ../ProcessedData/TE_CG_DMR_distances.bed.gz