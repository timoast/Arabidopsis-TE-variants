#! /bin/bash

gzip -dc ../RawData/TAIR10_GFF3_genes.gff.gz \
    | cut -f 9 \
    | grep Note \
    | awk -F'[;=]' 'BEGIN{OFS="\t"} {print $2, $4}' \
	  > ../ProcessedData/gene_list.tsv
