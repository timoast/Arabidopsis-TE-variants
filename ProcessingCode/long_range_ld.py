#! /usr/bin/env python

from __future__ import division
from argparse import ArgumentParser
import gzip


parser = ArgumentParser(description="Check for long-range LD in set of TE variants")
parser.add_argument('-t', '--transposons', help='file containing the TE variants', required=True)
parser.add_argument('-c', '--cutoff', help='fraction of matching accession names to be called similar', required=True)
parser.add_argument('-o', '--outfile', help='name of output file', required=True)
options = parser.parse_args()


def read_file(options):
    """ read TE variants into list """
    if options.transposons.endswith(".gz"):
        infile = gzip.open(options.transposons, "rb")
    else:
        infile = open(options.transposons, "r")
    d = []
    for line in infile:
        line = line.rsplit()
        if line[0] == "chromosome":
            pass  # header
        else:
            pos_acc = line[4]
            neg_acc = line[5]
            ident = ",".join([line[0], line[1], line[2]])
            d.append([pos_acc, neg_acc, ident])
    infile.close()
    return(d)


def compare_lists(list1, list2):
    """take two lists, return fraction of matches between the two lists"""
    number_match = len(set(list1) & set(list2))
    l1 = len(list1)
    l2 = len(list2)
    fraction = number_match / max(l1, l2)
    return(fraction)


def calc_distance(c1, c2):
    """ Find the distance between coordinates if chromosome is the same
        otherwise just return NA """
    if c1[0] == c2[0]:  # chromosome is the same
        d = abs(int(c1[1]) - int(c2[1]))
    else:
        d = "NA"
    return(d)
        

def get_similar(data, cutoff):
    """iterate over dictionary items, return items with similarity above cutoff"""
    similar_variants = []  # holds coordinates of similar TE variants
    # convert to only half the symmetrical matrix so we don't do 2x too many comparisons and get duplicate results
    for i in xrange(0, len(data)):
        for j in xrange(i, len(data)):
            acc_pos_1, acc_neg_1, coords_1 = [x.split(",") for x in data[i]]
            acc_pos_2, acc_neg_2, coords_2 = [x.split(",") for x in data[j]]
            if (coords_1 == coords_2) or (len(acc_pos_1) < 5) or (len(acc_neg_1) < 5):  # same TEPAV or only one accession
                pass
            else:
                # look at minor allele
                if len(acc_pos_1) < len(acc_neg_1):
                    fraction = compare_lists(acc_pos_1, acc_pos_2)
                else:
                    fraction = compare_lists(acc_neg_1, acc_neg_2)
                if fraction >= float(cutoff):
                    d = calc_distance(coords_1, coords_2)
                    all_info = coords_1 + coords_2 + [str(d)]
                    similar_variants.append(all_info)
                else:
                    pass
    return(similar_variants)


def main(options):
    data = read_file(options)
    highly_similar_pairs = get_similar(data, options.cutoff)
    with open(options.outfile, "w+") as outfile:
        for i in highly_similar_pairs:
            outfile.write("\t".join(i) + "\n")

                      
main(options)
