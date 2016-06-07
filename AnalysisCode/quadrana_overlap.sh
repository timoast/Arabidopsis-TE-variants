# finds number of TE variants reported by Quadrana et al. are within 1 kb of our TE variants.

gunzip -c ../RawData/TEPID_TEPAV.tsv.gz \
| sed 1d \
| bedtools closest -d -b - -a ../RawData/Quadrana_TEPAV.tsv.gz \
| awk '$16 < 1000 {print $0}' \
| wc -l