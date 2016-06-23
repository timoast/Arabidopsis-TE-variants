#! /bin/sh

gzip -dc ../RawData/TEPID_TEPAV.tsv.gz \
| sed 1d \
| bedtools intersect -a ../RawData/KEE_regions_150kb.bed -b - -c \
> ../ProcessedData/KEE_intersection_count.tsv

python choose_random_coords.py 10000 300000 \
| sort -k1,1 -k2,2n - \
> ../ProcessedData/kee_random_10k.bed

gzip -dc ../RawData/TEPID_TEPAV.tsv.gz \
| sed 1d \
| bedtools intersect -a ../ProcessedData/kee_random_10k.bed -b - -c \
> ../ProcessedData/kee_random_intersections.bed