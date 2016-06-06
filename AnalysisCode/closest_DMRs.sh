gunzip -c ../RawData/TEPID_TEPAV.tsv.gz | sed 1d \
| bedtools closest -b - -a ../RawData/c_dmrs.tsv.gz -d -wb \
> ../ProcessedData/TE_C_DMR_distances.bed

gunzip -c ../RawData/TEPID_TEPAV.tsv.gz | sed 1d \
| bedtools closest -b - -a ../RawData/cg_dmrs.tsv.gz -d -wb \
> ../ProcessedData/TE_CG_DMR_distances.bed