#! /usr/bin/python

from __future__ import division
import MySQLdb
from argparse import ArgumentParser
import gzip

parser = ArgumentParser(description='Get mC within +/- 200 bp flanking regions of input coordinates for all samples')
parser.add_argument('-x', '--host', help='host', required=True)
parser.add_argument('-u', '--user', help='mySQL user', required=True)
parser.add_argument('-p', '--password', help='mySQL password', required=True)
parser.add_argument('-d', '--database', help='mySQL database name', required=True)
parser.add_argument('-f', '--infile', help='input coordinate file', required=True)
options = parser.parse_args()

class mcDatabase:

    def __init__(self):
        self.link = None

    def setupSQL(self, host, user, password, database):
        """ setup connection with mysql """
        link = MySQLdb.connect(host, user, password)
        self.link = link
        cursor = link.cursor()
        cursor.execute("use {};".format(database))
        cursor.execute("show tables;")
        results = cursor.fetchall()
        self.tables = [row[0] for row in results]
        cursor.close()

    def gatherDataFromMySQL(self, line, tables):
        """ collect mC data from each table in database between coordinates in line """
        cursor = self.link.cursor()
        data = []
        upstream = int(line[1]) - 200
        downstream = int(line[2]) + 200
        for i in tables:
            query = "select sum(mc), sum(h) from {db} where assembly = {ch} and (position between {u} and {s} or position between {e} and {d})".format(
                db=str(i),
                ch=line[0].strip('chr'),
                u=upstream,
                s=line[1],
                e=line[2],
                d=downstream)
            cursor.execute(query)
            results = cursor.fetchall()
            mc = processResults(results)
            data.append(mc)
        cursor.close()
        return data


def processResults(results):
    for row in results:
        mc = row[0]
        h = row[1]
    if h == 0 or h is None:
        return None
    elif h > 0:
        level = float(mc / h)
        return round(level, 4)
    else:
        raise Exception('h values in table incorrect')


def filter_tables(inp):
    d = []
    for i in inp:
        if i.endswith("_bud"):
            pass
        else:
            d.append(i)
    return d


mc = mcDatabase()
mc.setupSQL(options.host, options.user, options.password, options.database)
mc.tables.remove("snps")
tables = filter_tables(mc.tables)
print 'chr\tstart\tstop\t'+'\t'.join(tables)  # header

if options.infile.endswith(".gz"):
    inf = gzip.open(options.infile, "rb")
else:
    inf = open(options.infile, "r")
for line in inf:
    line = line.rsplit()
    chrom = line[0]
    start = line[1]
    stop = line[2]
    coords = [chrom, start, stop]
    data = mc.gatherDataFromMySQL(line, tables)
    data = (str(x) for x in data)
    print '\t'.join(coords) + '\t' + '\t'.join(data)

