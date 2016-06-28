#! /bin/sh

# finds coverages of all samples (from tepid-discover log file)
# run from directory containing all sample subdirectories

head */tepid_discover_log* \
| python get_cov.py - \
> ../ProcessedData/coverages.tsv

wc -l */insertions* \
| python get_counts.py - \
| sed '$ d' \
> ../ProcessedData/insertion_counts.tsv

wc -l */deletions* \
| python get_counts.py - \
| sed '$ d' \
> ../ProcessedData/deletion_counts.tsv

wc -l */second_pass_insertion* \
| python get_counts.py -\
| sed '$ d' \
> ../ProcessedData/refinement_stats.tsv