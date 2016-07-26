#! /usr/bin/env python

# run from directory containing mC_calls files

import os
import MySQLdb
from argparse import ArgumentParser


"""
MySQL commands:
create database population_epigenetics;
use population_epigenetics;
create table mC_calls_Aa_0 (assembly int(1), position int(9), strand char(1), class char(3), mc int(3), h int(3), uc int(3));
insert into mC_calls_Aa_0 (assembly, position, strand, class, mc, h, uc) values (1, 10, '+', 'CHH', 4, 5, 1);
...
alter table mC_calls_Aa_0 add index (position);
"""

def create_database(options):
    link = MySQLdb.connect(options.host, options.user, options.password)
    cursor = link.cursor()
    query = "create database {db};".format(db=options.database)
    cursor.execute(query)
    cursor.close()
    link.close()

def add_table(options, fname):
    """Make new table and add data"""
    with open(fname, 'r') as infile:
        link = MySQLdb.connect(options.host, options.user, options.password)
        cursor = link.cursor()
        accession_name = fname.replace('calls_', '.').split('.')
        accession_name = accession_name[1]
        cursor.execute("use {db};".format(db=options.database))
        build = """create table mC_calls_{name} (assembly int(1), position int(9), strand char(1), class char(3), mc int(3), h int(3), uc int(3));""".format(name=accession_name)
        cursor.execute(build)
        for line in infile:
            line = line.rsplit()
            chrom = line[0]
            if chrom == 'chrom':
                pass  # header
            else:
                chrom = int(chrom)
                pos = int(line[1])
                strand = line[2]
                mc_class = line[3]
                if mc_class[:2] == 'CG':
                    mc_class = 'CG'
                else:
                    mc_class = mc_class[0] + mc_class[1:].replace('A', 'H').replace('T', 'H').replace('C', 'H')
                mc = int(line[4])
                h = int(line[5])
                uc = h - mc
                add_data = """insert into mC_calls_{name} (assembly, position, strand, class, mc, h, uc) values ({chrom}, {pos}, '{strand}', '{mc_class}', {mc}, {h}, {uc});""".format(name=accession_name,
                                                                                                                                                                                       chrom=chrom,
                                                                                                                                                                                       pos=pos,
                                                                                                                                                                                       strand=strand,
                                                                                                                                                                                       mc_class=mc_class,
                                                                                                                                                                                       mc=mc,
                                                                                                                                                                                       h=h,
                                                                                                                                                                                       uc=uc)
                cursor.execute(add_data)
        cursor.close()
        link.close()


def indexTables(options):
    link = MySQLdb.connect(options.host, options.user, options.password)
    cursor = link.cursor()
    acc_data = {}
    col_data = {}
    cursor.execute("use {};".format(options.database))
    cursor.execute("show tables")
    tables_result = cursor.fetchall()
    names = []
    for row in tables_result:
        name = row[0]
        names.append(name)
    for name in names:
        query = "alter table {n} add index (position)".format(n=name)
        cursor.execute(query)
    cursor.close()
    link.close()


parser = ArgumentParser(description='Build MySQL tables from file')
parser.add_argument('--host', help='MySQL host address', required=True)
parser.add_argument('-u', '--user', help='MySQL user name', required=True)
parser.add_argument('-p', '--password', help='MySQL password', required=True)
parser.add_argument('-d', '--database', help='MySQL database name', required=True)
options = parser.parse_args()

create_database(options)

for files in os.listdir('.'):
    if files.startswith("GSM"):
        add_table(options, file)

indexTables(options)

