#! /usr/bin/env python

import sys


def getCoverages(inf):
    for line in inf:
        if line.startswith("Sample"):
            nm = line.split()[1].strip("\n")
        elif line.startswith("coverage ="):
            cov = line.split()[2].split("x")[0]
            print(nm + "\t" + cov)
            cov = None
            nm = None
        else:
            pass

getCoverages(sys.stdin)