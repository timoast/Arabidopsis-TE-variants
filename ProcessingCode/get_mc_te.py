#! /usr/bin/env python

from __future__ import division
import MySQLdb
from argparse import ArgumentParser
import gzip

"""
Get DNA methylation levels around regions in a file of input coordinates
Input files should have chromosome (matching chromosome format in MySQL table), start, stop, feature name (eg name of TE)
This file can be gzipped
Will try to guess if there is a header present in the file and skip it, but it's best not to have a header
You can specify the number and size of bins to get mC data for, and the mC context
Default is 20 bins of 100 bp upstream and downstream (feature +/- 2 kb), and 20 bins of equal size within the feature itself
Default is to get all mC in all contexts
The output will be printed to stdout, so you can pipe to another program or to a file
"""

parser = ArgumentParser(description='gather mC data in n bins from mySQL database between given intervals, +/- flanking regions')
parser.add_argument('-a', '--host', help='host', required=True)
parser.add_argument('-u', '--user', help='mySQL user', required=True)
parser.add_argument('-p', '--password', help='mySQL password', required=True)
parser.add_argument('-d', '--database', help='mySQL database name', required=True)
parser.add_argument('-t', '--table', help='mySQL table name', required=True)
parser.add_argument('-i', '--infile', help='input coordinate file', required=True)
parser.add_argument('-c', '--context', help='mC context', required=False, default=False)
parser.add_argument('-b', '--binsize', help='Bin size', required=False, default=100, type=int)
parser.add_argument('-n', '--numberbins', help='Number of bins', required=False, default=20, type=int)
options = parser.parse_args()


def get_data(chrom, start, stop, options, cursor):
    """
    fetch mC data from table in database in between given coordinates plus upstream and downstream regions
    input set of coordinates, MySQL database, table and mC context
    average methylation in n bins within the given coordinates
    average methylation in upstream and downstream regions
    """
    window = options.binsize * options.numberbins  # size of flanking region to query each side, default 2 kb
    upstream = start - window
    levels = query_region(upstream, start, stop, options, chrom, cursor)
    return levels


def query_region(upstream, start, stop, options, chrom, cursor):
    total_bins = options.numberbins*3
    data = []
    # need to do three parts: upstream, body, downstream
    for x in xrange(total_bins):    # between 0 and 60
        if x < options.numberbins:  # between 0 and 20
            # upstream region
            mc_start = upstream + (x*options.binsize)           # starts 0 offset (0 x binsize)
            mc_end = upstream + ((x+1)* options.binsize)        # stops at upstream + 1xbinsize
        elif options.numberbins < x < options.numberbins*2:  # between 20 and 40
            # body
            # binsize is different for body (non-standard length)
            a = int(x-(options.numberbins))       # start 0: x = 20, x - numberbins (20) = 0
            body_binsize = (stop - start) / options.numberbins  # will be the length of feature evenly divided into 20 bins
            mc_start = start + (a*body_binsize)   # starts with 0 offset, then increase by 20 each step
            mc_end = start + ((a+1)*body_binsize) # one binsize larger
        else:  # between 40 and 60
            # downstream region
            a = int(x-(options.numberbins*2))        # starts 0: 40 - numberbins (20) * 2 = 0
            mc_start = stop + (a*options.binsize)    # starts with stop == stop, increase by 1 binsize each step
            mc_end = stop + ((a+1)* options.binsize) # one binsize larger
        if options.context is not False:
            query = "select sum(mc), sum(h) from {db} where context = '{c}' and assembly = '{chrom}' and (position between {mc_start} and {mc_end});".format(
                db=options.database + '.' + options.table,
                chrom=chrom,
                c=options.context,
                mc_start=mc_start,
                mc_end=mc_end)
        else:
            query = "select sum(mc), sum(h) from {db} where assembly = '{chrom}' and (position between {mc_start} and {mc_end});".format(
                db=options.database + '.' + options.table,
                chrom=chrom,
                mc_start=mc_start,
                mc_end=mc_end)
        cursor.execute(query)
        results = cursor.fetchall()
        for row in results:
            mc = row[0]
            h = row[1]
        if h == 0 or h is None:
            level = 0.
        else:
            level = float(mc / h)
            if level < 0.:
                level = 0.
        data.append(level)
    return data


def setup_sql(options):
    link = MySQLdb.connect(options.host, options.user, options.password)
    cursor = link.cursor()
    return cursor, link


def process_all(options):
    cursor, link = setup_sql(options)
    if options.infile.endswith(".gz"):
        coordinate_file = gzip.open(options.infile, "rb")
    else:
        coordinate_file = open(options.infile, "r")
    bins = ['bin_'+str(x) for x in range(int(options.numberbins*3))]
    print('coords'+"\t"+"\t".join(bins)+"\n")  # header
    window = options.binsize * options.numberbins
    for line in coordinate_file:
        line = line.rsplit()
        if line[0] == "chromosome":
            pass # header
        else:
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            feature_name = line[3]
            if ((start-window) <= 0) or chrom == "Mt" or chrom == "Pt":
                pass
            else:
                mc_values = get_data(chrom, start, stop, options, cursor)
                rounded = [round(x, 4) for x in mc_values]
                print(feature_name+"\t"+"\t".join(map(str, rounded))+"\n")
    cursor.close()
    link.close()
    coordinate_file.close()

process_all(options)
