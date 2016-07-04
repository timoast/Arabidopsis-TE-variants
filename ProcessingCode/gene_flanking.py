#! /usr/bin/env python


def flanking(options):
    if options.file.endswith(".gz"):
        infile = gzip.open(options.file, "rb")
    else:
        infile = open(options.file, "r")
    for line in infile:
        line = line.rsplit()
        if line[2] == "gene":
            chrom = line[0]
            strand = line[6]
            info = line[8]
            gene = info.split(';')[0].split("=")[1]
            if (options.region == "upstream" and strand == "+") or (options.region == "downstream" and strand == "-"):
                start = str(max(0, int(line[3])-2000))
                stop = line[3]
            elif (options.region == "upstream" and strand == "-") or (options.region == "downstream" and strand == "+"):
                start = line[4]
                stop = str(int(line[4]) + 2000)
            else:
                raise Exception("Error: incorrect strand information in gff file")
            print(chrom+"\t"+start+"\t"+stop+"\t"+gene+"\t"+strand+"\t"+"TAIR10"+"\t"+options.region+"\t"+info)
        else:
            pass
    infile.close()


if __name__ == "__main__":
    import gzip
    from argparse import ArgumentParser


    parser = ArgumentParser(description='Get gene upstream or downstream regions from gff file')
    parser.add_argument('-r', '--region', help='upstream or downstream', required=True)
    parser.add_argument('-f', '--file', help='gff file', required=True)
    options = parser.parse_args()
    flanking(options)