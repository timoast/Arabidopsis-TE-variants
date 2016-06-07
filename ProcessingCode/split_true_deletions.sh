# This script will split TE variants into two categories
# 1. TE absences that are due to the lack of an insertion in the ancestral genome
# 2. TE absences that are due to a true deletion of the ancestral TE
# Some TE variants will have an NA call, where the classification is ambiguous
# these will not be in either file

gunzip -c ../RawData/TEPID_TEPAV.tsv.gz \
| sed 1d \
| awk 'BEGIN {FS=OFS="\t"} {if ($10 == "True deletion") print $0}' - \
| gzip - \
> ../ProcessedData/TEPAV_true_deletions.tsv.gz

gunzip -c ../RawData/TEPID_TEPAV.tsv.gz \
| sed 1d \
| awk 'BEGIN {FS=OFS="\t"} {if ($10 == "No insertion") print $0}' - \
| gzip - \
> ../ProcessedData/TEPAV_no_insertion.tsv.gz