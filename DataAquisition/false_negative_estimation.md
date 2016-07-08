## False negative estimation for TE absences and TE insertions

Using coordinates of Col-0 specific TE insertions from Quadrana et al. 2016, eLife

### False negative estimation for TE absence calls.  

This is very straightforward, simply check how many of the Col-0 specific insertions are detected as absent in Ler using the Ler short read data.  

```bash
bedtools intersect -a quadrana_col_te.bed -b deletions_ler.bed | wc -l
## 173
```

173 out of 193 TE absences in Ler are detected with TEPID (89.6%)

### False negative estimation for TE insertion calls.  

This is more complicated, we will  map Col-0 reads to the Ler assembly, and should be able to pick detect the Col-0 specific TE insertions as insertions relative to the Ler genome using Col-0 short reads.  

I will use 2x90 PE reads from Jiang et al. 2014, Genome Research. I pooled together a couple of Col-0 control replicates to get the final dataset.  

```bash
cat Col-0_Control_G10-1_1.fastq Col-0_Control_G10-2_1.fastq | pigz - > col_1.fq.gz
cat Col-0_Control_G10-1_2.fastq Col-0_Control_G10-2_2.fastq | pigz - > col_2.fq.gz
```

Build the Ler genome index for bowtie2 and yaha

```bash
bowtie2-build polished_assembly.fasta bt2
yaha -g polished_assembly.fasta
```

Map using tepid-map  

```bash
tepid-map -x bt2 -y polished_assembly.X15_01_65525S -p 20 -s 200 -n ler \
-1 col_1.fq.gz -2 col_2.fq.gz
```

Now find TE insertions. Since we don't have a TE annotation for Ler, we just use a dummy annotation and look at the intermediate files that TE calls are made from.  

```bash
touch dummy_te.bed
tepid-discover -n ler -c ler.bam -s ler.split.bam -t dummy_te.bed -p 10 -k
```

Now check what the read count thresholds are for calling a TE insertion, given a TE annotation 

```bash
grep coverage tepid_discover_log_ler.txt
## coverage = 0x
```

As the coverage is estimated by TEPID using the first chromosome, this is one of the short scaffolds and gives an inaccurate estimate. Instead, we will estimate using one of the larger scaffolds, using the same code as is used by TEPID.  


```python
#! /usr/bin/env python

"""
calc_cov.py
"""

import pysam
import sys

def calc_cov(bam_name, start, stop):
    """
    calculates average coverage
    """
    bam = pysam.AlignmentFile(bam_name)
    x = 0
    l = 0
    for read in bam.pileup("scf7180000001373|quiver", start, stop):
        x += read.n
        l += 1
    if l > 0 and x > 0:
        bam.close()
        return int(x/l)
    else:
        bam.close()
        return 0

cov = calc_cov(sys.argv[1], 100000, 120000)
print(cov)
```

```bash
python calc_cov.py ler.bam
## 58
```

This will set the thresholds for insertion calls to >11 split reads. So we look for >11 split reads at the TE insertion sites, in the intermediate file.  

```bash
bedtools intersect -c -a quadrana_ler_te.bed -b disc_split_ins.temp \
| awk '$4 > 11' | wc -l
## 173
```

173 out of the 193 Col-0 insertions found in the 10 largest scaffolds (89.6%) have enough
read evidence to be called as an insertion by TEPID, had a high-quality TE annotation been
given. This amounts to a false negative rate of 10%.
