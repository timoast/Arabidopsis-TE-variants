#! /usr/bin/env python

""" Example GFF
chr1    TAIR9   exon    3631    3913    .   +   .   Parent=AT1G01010.1
chr1    TAIR9   exon    3996    4276    .   +   .   Parent=AT1G01010.1
chr1    TAIR9   exon    4486    4605    .   +   .   Parent=AT1G01010.1
chr1    TAIR9   exon    4706    5095    .   +   .   Parent=AT1G01010.1
chr1    TAIR9   exon    5174    5326    .   +   .   Parent=AT1G01010.1
chr1    TAIR9   exon    5439    5899    .   +   .   Parent=AT1G01010.1
chr1    TAIR9   exon    8571    8737    .   -   .   Parent=AT1G01020.1
chr1    TAIR9   exon    8417    8464    .   -   .   Parent=AT1G01020.1
chr1    TAIR9   exon    8236    8325    .   -   .   Parent=AT1G01020.1
chr1    TAIR9   exon    7942    7987    .   -   .   Parent=AT1G01020.1
"""

"""
Usage:
python add_feature_between.py infile.gff > output.gff
Read gff file and coordinates of intervening sequences and label with given feature name
Input must be position sorted
sort -k1,1 -k4,4n file > out
"""

import pandas as pd
import sys


def add_introns(f):
    with open(f, 'r') as f:
        table = pd.read_csv(f, sep='\t', header=-1)
    l = len(table)
    for x in range(1,l):
        if table.iloc[x][8] == table.iloc[x-1][8]:
            start = table.iloc[x-1][4]
            end = table.iloc[x][3]
            strand = table.iloc[x][6]
            chrom = table.iloc[x][0]
            data = [
                chrom,
                table.iloc[x][1],
                'intron',
                start,
                end,
                table.iloc[x][5],
                table.iloc[x][6],
                table.iloc[x][7],
                table.iloc[x][8]]
            print("\t".join(map(str,data)))


add_introns(sys.argv[1])
