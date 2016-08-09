#! /bin/bash

# make a temp folder to hold the intermediate output
mkdir ./temp

# filter where distance between TEs is more than 5 kb (this might be a low cutoff)
awk '$7 > 5000' ../ProcessedData/long_range_ld_te.tsv > ./temp/long_range_ld_5kb.tsv

# switch order of columns and cat, so we can intersect with the significant pval coordinates
awk 'BEGIN {FS=OFS="\t"} {print $4,$5,$6,$1,$2,$3,$7}' ./temp/long_range_ld_5kb.tsv \
    | cat - ./temp/long_range_ld_5kb.tsv \
    > ./temp/both_combinations_long_range_ld.tsv

# first need to change column order and add a second coordinate to the pvals file
awk 'BEGIN {FS=OFS="\t"} {$6 = $4 + 10; print $3,$4,$6,$1,$2,$5}' ../ProcessedData/cdmr_sig_pvals.tsv > ./temp/reordered_pvals.tsv

# intersect with the significant pvalues
bedtools intersect -a ./temp/reordered_pvals.tsv -b ./temp/both_combinations_long_range_ld.tsv -wb > ./temp/intersections_pval.tsv

# now add column with distance between DMR and TE, NA if different chromosome
# filter where distance is less than 2 kb

# columns are coords of DMR-correlated TE insertion, insertion/deletion classification
# chr start of correlated DMR, coords of TE intersecting the inserted TE, coords of TE found
# in same accessions as the correlated TE, distance between the two TEs
# found in the same accessions (NA if chromosomes are different)

# if chromosomes different insert NA, else calculate the distance between correlated DMR and TE, print those where distance is less than 2000
awk 'BEGIN {FS=OFS="\t"} {if ($4 != $10) $14 = "NA"; else $14 = sqrt(($5 - $11)^2); if ($14 < 2000) print $0}' ./temp/intersections_pval.tsv \
    > ../ProcessedData/ghost_trans_associations.tsv

rm -r ./temp
