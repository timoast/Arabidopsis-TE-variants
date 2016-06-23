#! /usr/bin/env python

import sys

def count_variants(inf):
    for line in inf:
        line = line.split()
        count = line[0]
        accession = line[1].split("/")[0]
        print(accession + "\t" + count)

count_variants(sys.stdin)