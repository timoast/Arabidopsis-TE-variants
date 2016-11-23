library(RSVSim)

# load genome
# first download genome from TAIR
tair10 = readDNAStringSet("../../RawData/tair10.fa.gz")

# load data and read into correct format
# Must run the python code first to make these files
del_tair_r1 <- read.table("tair_del_r1.bed", header=TRUE)
del_tair_r2 <- read.table("tair_del_r2.bed", header=TRUE)
del_tair_r3 <- read.table("tair_del_r3.bed", header=TRUE)
del_tair_r4 <- read.table("tair_del_r4.bed", header=TRUE)
del_tair_r5 <- read.table("tair_del_r5.bed", header=TRUE)

ins_tair_r1 <- read.table("tair_ins_r1.bed", header=TRUE)
ins_tair_r2 <- read.table("tair_ins_r2.bed", header=TRUE)
ins_tair_r3 <- read.table("tair_ins_r3.bed", header=TRUE)
ins_tair_r4 <- read.table("tair_ins_r4.bed", header=TRUE)
ins_tair_r5 <- read.table("tair_ins_r5.bed", header=TRUE)

ath_deletions_r1 <- with(del_tair_r1, GRanges(chr, IRanges(start, stop)))
ath_deletions_r2 <- with(del_tair_r2, GRanges(chr, IRanges(start, stop)))
ath_deletions_r3 <- with(del_tair_r3, GRanges(chr, IRanges(start, stop)))
ath_deletions_r4 <- with(del_tair_r4, GRanges(chr, IRanges(start, stop)))
ath_deletions_r5 <- with(del_tair_r5, GRanges(chr, IRanges(start, stop)))

ath_insertions_r1 <- with(ins_tair_r1, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r2 <- with(ins_tair_r2, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r3 <- with(ins_tair_r3, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r4 <- with(ins_tair_r4, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r5 <- with(ins_tair_r5, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))

# simulate variation
simulateSV(output="./ath/r1",
           genome=tair10,
           regionsDels=ath_deletions_r1,
           regionsIns=ath_insertions_r1,
           random=FALSE)

simulateSV(output="./ath/r2",
           genome=tair10,
           regionsDels=ath_deletions_r2,
           regionsIns=ath_insertions_r2,
           random=FALSE)

simulateSV(output="./ath/r3",
           genome=tair10,
           regionsDels=ath_deletions_r3,
           regionsIns=ath_insertions_r3,
           random=FALSE)

simulateSV(output="./ath/r4",
           genome=tair10,
           regionsDels=ath_deletions_r4,
           regionsIns=ath_insertions_r4,
           random=FALSE)

simulateSV(output="./ath/r5",
           genome=tair10,
           regionsDels=ath_deletions_r5,
           regionsIns=ath_insertions_r5,
           random=FALSE)