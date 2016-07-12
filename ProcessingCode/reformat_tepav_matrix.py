#! /usr/bin/env python

def convert_matrix(infile, names):
    """
    convert TEPAV genotyped file to presence / absence matrix
    """
    if infile.endswith(".gz"):
        inf = gzip.open(infile, "rb")
    else:
        inf = open(infile, "r")
    for line in inf:
        line = line.rsplit()
        if line[0] == "chromosome":
            pass  # header
        else:
            chrom = line[0]
            start = line[1]
            stop = line[2]
            TE = line[3]
            n_te = str(len(TE.split(",")))
            pos = line[4].split(",")
            neg = line[5].split(",")
            te_id = "\t".join([chrom, start, stop])
            status = get_status(pos, neg, names)
            column_ordered = []
            for i in names:
                column_ordered.append(status[i])
            noNA = filter(lambda x: x != "NA", status.values())            
            noNA = map(int, noNA)
            pos_count = sum(noNA)
            l = len(noNA)
            neg_count = l - pos_count
            print(te_id + "\t" + TE + "\t" + n_te + "\t"  + str(pos_count) + "\t" + str(neg_count) + "\t" + "\t".join(column_ordered))
    inf.close()


def create_names_dict(infile):
    """
    read accession names into list
    """
    return [name.strip("\n") for name in open(infile, "r")]


def get_status(pos, neg, names):
    """
    create dictionary of TEPAV status for each accession (1/0/NA)
    """
    status = {}
    for i in names:
        if i in pos:
            status[i] = "1"
        elif i in neg:
            status[i] = "0"
        else:
            status[i] = "NA"
    return status


if __name__ == "__main__":
    from argparse import ArgumentParser
    import gzip

    parser = ArgumentParser(description='Convert TEPAV genotyped file to presence/absence matrix')
    parser.add_argument('-n', '--names', help='All sample names newline separated', required=True)
    parser.add_argument('-i', '--infile', help="Input TEPAV file", required=True)
    options = parser.parse_args()

    names = create_names_dict(options.names)
    print "chrom" +"\t" + "pos" + "\t" + "stop" + "\t" + "TEID" + "\t" + "number_of_te" + "\t" + "accession_poly_count" + "\t" + "accession_nopoly_count" + "\t" + "\t".join(names)
    convert_matrix(options.infile, names)