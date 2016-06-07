#! /usr/bin/env python

from __future__ import division
import MySQLdb
import numpy as np
from argparse import ArgumentParser
import gzip


parser = ArgumentParser(description='gather mC data from mySQL database between given intervals')
parser.add_argument('-a', '--host', help='host', required=True)
parser.add_argument('-u', '--user', help='mySQL user', required=True)
parser.add_argument('-p', '--password', help='mySQL password', required=True)
parser.add_argument('-d', '--database', help='mySQL database name', required=True)
parser.add_argument('-i', '--infile', help='input coordinate file', required=True)
parser.add_argument('-o', '--outfile', help='output file', required=True)
parser.add_argument('-c', '--context', help='mC context', required=False, default=False)
parser.add_argument('-b', '--binsize', help='Bin size', required=False, default=200, type=int)
parser.add_argument('-n', '--numberbins', help='Number of bins', required=False, default=20, type=int)
parser.add_argument('-m', '--all_mc_tables', help='all mC table names', required=True)
options = parser.parse_args()


def get_data(chrom, start, stop, options, cursor, tables):
    """
    fetch mC data from table in database in upstream and downstream regions from deletion sites
    input set of deletion coordinates with list of accessions containing deletion
    sum methylation in bins upstream and downstream of deletion coordinates
    return aggregate data for accessions containing and not containing TE deletion
    """
    data = []
    window = options.binsize * options.numberbins  # default 4 kb
    upstream = int(start - (window/2))
    x = 0
    levels = np.zeros(options.numberbins)
    for table in tables:
        x += 1
        levels += query_region(upstream, stop, options, table, chrom, cursor)
    levels = levels / x
    return levels.tolist()


def query_region(upstream, stop, options, table, chrom, cursor):
    data = np.zeros(options.numberbins)
    # need to do two parts, upstream and downstream
    for x in xrange(options.numberbins):
        if x < options.numberbins/2:
            # upstream region
            mc_start = upstream + (x*options.binsize)
            mc_end = upstream + ((x+1)* options.binsize)
        else:
            # downstream region
            a = int(x-(options.numberbins/2))
            mc_start = stop + (a*options.binsize)
            mc_end = stop + ((a+1)* options.binsize)
        if options.context is not False:
            query = "select sum(mc), sum(h) from {db} where class = '{c}' and assembly = {chrom} and (position between {mc_start} and {mc_end});".format(
                db=options.database + '.' + table,
                chrom=chrom,
                c=options.context,
                mc_start=mc_start,
                mc_end=mc_end)
        else:
            query = "select sum(mc), sum(h) from {db} where assembly = {chrom} and (position between {mc_start} and {mc_end});".format(
                db=options.database + '.' + table,
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
        data[x] += level
    return data


def setup_sql(options):
    link = MySQLdb.connect(options.host, options.user, options.password)
    cursor = link.cursor()
    return cursor, link


def group_tables(pos_accessions, neg_accessions, options):
    table_prefix = 'mC_calls_'
    prefixed_pos_accessions = [table_prefix+i for i in pos_accessions]
    prefixed_neg_accessions = [table_prefix+i for i in neg_accessions]
    all_mc_tables = [line.strip() for line in open(options.all_mc_tables, 'r')]
    no_bud_tables = []
    for i in all_mc_tables:
        if i.endswith('bud'):
            pass
        else:
            no_bud_tables.append(i)
    pos_tables = [i for i in prefixed_pos_accessions if i in no_bud_tables]
    neg_tables = [i for i in prefixed_neg_accessions if i in no_bud_tables]
    return pos_tables, neg_tables


def process_all(options):
    cursor, link = setup_sql(options)
    if options.infile.endswith(".gz"):
        insertions = gzip.open(options.infile, "rb")
    else:
        insertions = open(options.infile, "r")
    with open(options.outfile, 'w+') as outf:
        ins_bins = ['ins'+str(x) for x in range(int(options.numberbins))]
        no_ins_bins = ['no_ins'+str(x) for x in range(int(options.numberbins))]
        outf.write('coords'+"\t"+"\t".join(ins_bins)+"\t"+"\t".join(no_ins_bins)+"\n")  # header
        window = options.binsize * options.numberbins
        for line in insertions:
            line = line.rsplit()
            if line[0] == "chromosome":
                pass # header
            else:
                chrom = line[0].strip('chr')
                start = int(line[1])
                stop = int(line[2])
                coords = 'chr'+chrom+','+str(start)+','+str(stop)
                pos_accessions = line[4].replace('-', '_').split(',')
                neg_accessions = line[5].replace('-', '_').split(',')
                pos_tables, neg_tables = group_tables(pos_accessions, neg_accessions, options)
                if ((start-window) <= 0) or len(pos_tables) == 0 or len(neg_tables) == 0 or chrom == "Mt" or chrom == "Pt":
                    pass
                else:
                    mc_values_accessions = get_data(chrom, start, stop, options, cursor, pos_tables)
                    mc_values_non_accessions = get_data(chrom, start, stop, options, cursor, neg_tables)
                    all_data = mc_values_accessions + mc_values_non_accessions
                    outf.write(coords+"\t"+"\t".join(map(str,all_data))+"\n")
    cursor.close()
    link.close()
    insertions.close()

process_all(options)
