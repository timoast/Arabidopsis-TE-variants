#! /usr/bin/python

from __future__ import division
import MySQLdb
from argparse import ArgumentParser
import gzip

"""
usage: get_mc_DMRs_context.py [-h] -x HOST -u USER -p PASSWORD -d DATABASE -f
                              FILE -c CONTEXT

Get mC within coordinates for multiple samples

optional arguments:
  -h, --help            show this help message and exit
  -x HOST, --host HOST  host
  -u USER, --user USER  mySQL user
  -p PASSWORD, --password PASSWORD
                        mySQL password
  -d DATABASE, --database DATABASE
                        mySQL database name
  -f FILE, --file FILE  input coordinate file
  -c CONTEXT, --context CONTEXT
                        DNA methylation context
"""

parser = ArgumentParser(description='Get mC within coordinates for multiple samples')
parser.add_argument('-x', '--host', help='host', required=True)
parser.add_argument('-u', '--user', help='mySQL user', required=True)
parser.add_argument('-p', '--password', help='mySQL password', required=True)
parser.add_argument('-d', '--database', help='mySQL database name', required=True)
parser.add_argument('-f', '--file', help='input coordinate file', required=True)
parser.add_argument('-c', '--context', help='DNA methylation context', required=True)
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

    def gatherDataFromMySQL(self, line, context):
        """ collect mC data from each table in database between coordinates in line """
        cursor = self.link.cursor()
        data = []
        for i in self.tables:
            query = "select sum(mc), sum(h) from {db} where class in ('{context}') and assembly = {ch} and position between {s} and {e}".format(
                db=str(i),
                ch=line[0].strip('chr'),
                context=context,
                s=line[1],
                e=line[2])
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

mc = mcDatabase()
mc.setupSQL(options.host, options.user, options.password, options.database)
mc.tables.remove("snps")
names =  [x.split('mC_calls_')[1] for x in mc.tables]
print 'chr\tstart\tstop\t'+'\t'.join(names)  # header
if options.file.endswith(".gz"):
    infile = gzip.open(options.file, "rb")
else:
    infile = open(options.file, "r")
for line in infile:
    line = line.rsplit()
    data = mc.gatherDataFromMySQL(line, str(options.context))
    data = (str(x) for x in data)
    print '\t'.join(line) + '\t' + '\t'.join(data)
infile.close()